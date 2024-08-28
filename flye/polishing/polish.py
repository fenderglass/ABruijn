#(c) 2016 by Authors
#This file is a part of ABruijn program.
#Released under the BSD license (see LICENSE file)

"""
Runs polishing binary in parallel and concatentes output
"""

import logging
import subprocess, multiprocessing
import os
from collections import defaultdict
import gzip

from flye.polishing.alignment import (make_alignment, get_contigs_info,
                                      merge_chunks, split_into_chunks)
from flye.utils.sam_parser import SynchronizedSamReader
from flye.polishing.bubbles import make_bubbles
import flye.utils.fasta_parser as fp
from flye.utils.utils import which
import flye.config.py_cfg as cfg
from flye.six import iteritems


POLISH_BIN = "flye-modules"

logger = logging.getLogger()


class PolishException(Exception):
    pass


def check_binaries():
    if not which(POLISH_BIN):
        raise PolishException("polishing binary was not found. "
                              "Did you run 'make'?")
    try:
        devnull = open(os.devnull, "w")
        subprocess.check_call([POLISH_BIN, "polisher", "-h"], stderr=devnull)
    except subprocess.CalledProcessError as e:
        if e.returncode == -9:
            logger.error("Looks like the system ran out of memory")
        raise PolishException(str(e))
    except OSError as e:
        raise PolishException(str(e))


def polish(contig_seqs, read_seqs, work_dir, num_iters, num_threads, read_platform,
           read_type, output_progress):
    """
    High-level polisher interface
    """
    logger_state = logger.disabled
    if not output_progress:
        logger.disabled = True

    subs_matrix = os.path.join(cfg.vals["pkg_root"],
                               cfg.vals["err_modes"][read_platform]["subs_matrix"])
    hopo_matrix = os.path.join(cfg.vals["pkg_root"],
                               cfg.vals["err_modes"][read_platform]["hopo_matrix"])
    use_hopo = cfg.vals["err_modes"][read_platform]["hopo_enabled"]
    use_hopo = use_hopo and (read_type == "raw")
    stats_file = os.path.join(work_dir, "contigs_stats.txt")
    bed_coverage = os.path.join(work_dir, "base_coverage.bed.gz")

    bam_input = read_seqs[0].endswith("bam")

    prev_assembly = contig_seqs
    contig_lengths = None
    coverage_stats = None
    for i in range(num_iters):
        logger.info("Polishing genome (%d/%d)", i + 1, num_iters)

        ####
        if not bam_input:
            logger.info("Running minimap2")
            alignment_file = os.path.join(work_dir, "minimap_{0}.bam".format(i + 1))
            make_alignment(prev_assembly, read_seqs, num_threads,
                           read_platform, read_type, alignment_file)
        else:
            logger.info("Polishing with provided bam")
            alignment_file = read_seqs[0]

        #####
        logger.info("Separating alignment into bubbles")
        contigs_info = get_contigs_info(prev_assembly)
        bubbles_file = os.path.join(work_dir,
                                    "bubbles_{0}.fasta".format(i + 1))
        coverage_stats, mean_aln_error = \
            make_bubbles(alignment_file, contigs_info, prev_assembly,
                         read_platform, num_threads,
                         bubbles_file)

        logger.info("Alignment error rate: %f", mean_aln_error)
        consensus_out = os.path.join(work_dir, "consensus_{0}.fasta".format(i + 1))
        polished_file = os.path.join(work_dir, "polished_{0}.fasta".format(i + 1))
        if os.path.getsize(bubbles_file) == 0:
            logger.info("No reads were aligned during polishing")
            if not output_progress:
                logger.disabled = logger_state
            open(stats_file, "w").write("#seq_name\tlength\tcoverage\n")
            open(polished_file, "w")
            gzip.open(bed_coverage, "wt")
            return polished_file, stats_file

        #####
        logger.info("Correcting bubbles")
        _run_polish_bin(bubbles_file, subs_matrix, hopo_matrix,
                        consensus_out, num_threads, output_progress, use_hopo)
        polished_fasta, polished_lengths, bubble_coverages = _compose_sequence(consensus_out)
        fp.write_fasta_dict(polished_fasta, polished_file)

        #Cleanup
        os.remove(bubbles_file)
        os.remove(consensus_out)
        if not bam_input:
            os.remove(alignment_file)

        contig_lengths = polished_lengths
        prev_assembly = polished_file

    with open(stats_file, "w") as f:
        f.write("#seq_name\tlength\tcoverage\n")
        for ctg_id in contig_lengths:
            if ctg_id not in coverage_stats:
                coverage_stats[ctg_id] = 0
            f.write("{0}\t{1}\t{2}\n".format(ctg_id,
                    contig_lengths[ctg_id], coverage_stats[ctg_id]))

    with gzip.open(bed_coverage, "wt") as f:
        f.write("#seq_name\tstart\tend\tcoverage\n")
        for ctg_id in bubble_coverages:
            cur_pos = 0
            for (bub_len, bub_cov) in bubble_coverages[ctg_id]:
                f.write("{0}\t{1}\t{2}\t{3}\n".format(ctg_id, cur_pos, cur_pos + bub_len, bub_cov))
                cur_pos += bub_len

    if not output_progress:
        logger.disabled = logger_state

    return prev_assembly, stats_file


def generate_polished_edges(edges_file, gfa_file, polished_contigs, work_dir,
                            platform, read_type, polished_stats, num_threads):
    """
    Generate polished graph edges sequences by extracting them from
    polished contigs
    """
    logger.debug("Generating polished GFA")

    edges_new_coverage = {}
    with open(polished_stats, "r") as f:
        for line in f:
            if line.startswith("#"):
                continue
            ctg, _len, coverage = line.strip().split()
            ctg_id = ctg.split("_")[1]
            edges_new_coverage[ctg_id] = int(coverage)

    alignment_file = os.path.join(work_dir, "edges_aln.bam")
    polished_dict = fp.read_sequence_dict(polished_contigs)
    make_alignment(polished_contigs, [edges_file], num_threads,
                   platform, read_type, alignment_file)
    aln_reader = SynchronizedSamReader(alignment_file,
                                       polished_dict, multiprocessing.Manager(),
                                       cfg.vals["max_read_coverage"])
    aln_by_edge = defaultdict(list)

    #getting one best alignment for each contig
    #for ctg in polished_dict:
    #    ctg_aln = aln_reader.get_alignments(ctg)
    for aln in aln_reader.get_all_alignments():
        aln_by_edge[aln.qry_id].append(aln)
    #logger.debug("Bam parsing done")

    MIN_CONTAINMENT = 0.9
    updated_seqs = 0
    edges_dict = fp.read_sequence_dict(edges_file)
    for edge in edges_dict:
        if edge in aln_by_edge:
            aln_by_edge[edge].sort(key=lambda a: a.qry_end - a.qry_start, reverse=True)
            main_aln = aln_by_edge[edge][0]
            map_start = main_aln.trg_start
            map_end = main_aln.trg_end
            for aln in aln_by_edge[edge]:
                if aln.trg_id == main_aln.trg_id and aln.trg_sign == main_aln.trg_sign:
                    map_start = min(map_start, aln.trg_start)
                    map_end = max(map_end, aln.trg_end)

            new_seq = polished_dict[main_aln.trg_id][map_start : map_end]
            if main_aln.qry_sign == "-":
                new_seq = fp.reverse_complement(new_seq)

            #print(edge, main_aln.qry_len, len(new_seq), main_aln.qry_start, main_aln.qry_end)
            if len(new_seq) / aln.qry_len > MIN_CONTAINMENT:
                edges_dict[edge] = new_seq
                updated_seqs += 1

    #writes fasta file with polished egdes
    #edges_polished = os.path.join(work_dir, "polished_edges.fasta")
    #fp.write_fasta_dict(edges_dict, edges_polished)

    #writes gfa file with polished edges
    with open(os.path.join(work_dir, "polished_edges.gfa"), "w") as gfa_polished, \
         open(gfa_file, "r") as gfa_in:
        for line in gfa_in:
            if line.startswith("S"):
                seq_id = line.split()[1]
                coverage_tag = line.split()[3]
                seq_num = seq_id.split("_")[1]
                if seq_num in edges_new_coverage:
                    #logger.info("from {0} to {1}".format(coverage_tag, edges_new_coverage[seq_num]))
                    coverage_tag = "dp:i:{0}".format(edges_new_coverage[seq_num])
                gfa_polished.write("S\t{0}\t{1}\t{2}\n"
                                    .format(seq_id, edges_dict[seq_id], coverage_tag))
            else:
                gfa_polished.write(line)

    logger.debug("%d sequences remained unpolished",
                 len(edges_dict) - updated_seqs)
    os.remove(alignment_file)


def filter_by_coverage(args, stats_in, contigs_in, stats_out, contigs_out):
    """
    Filters out contigs with low coverage
    """
    SUBASM_MIN_COVERAGE = 1
    HARD_MIN_COVERAGE = cfg.vals["hard_minimum_coverage"]
    RELATIVE_MIN_COVERAGE = cfg.vals["relative_minimum_coverage"]

    ctg_stats = {}
    sum_cov = 0
    sum_length = 0

    with open(stats_in, "r") as f:
        for line in f:
            if line.startswith("#"): continue
            tokens = line.split("\t")
            ctg_id, ctg_len, ctg_cov = tokens[0], int(tokens[1]), int(tokens[2])
            ctg_stats[ctg_id] = (ctg_len, ctg_cov)
            sum_cov += ctg_cov * ctg_len
            sum_length += ctg_len

    mean_coverage = int(sum_cov / sum_length)
    coverage_threshold = None
    if args.read_type == "subasm":
        coverage_threshold = SUBASM_MIN_COVERAGE
    elif args.meta:
        coverage_threshold = HARD_MIN_COVERAGE
    else:
        coverage_threshold = int(round(mean_coverage /
                                       RELATIVE_MIN_COVERAGE))
        coverage_threshold = max(HARD_MIN_COVERAGE, coverage_threshold)
    logger.debug("Mean contig coverage: %d, selected threshold: %d",
                 mean_coverage, coverage_threshold)

    filtered_num = 0
    filtered_seq = 0
    good_fasta = {}
    for hdr, seq in fp.stream_sequence(contigs_in):
        if ctg_stats[hdr][1] >= coverage_threshold:
            good_fasta[hdr] = seq
        else:
            filtered_num += 1
            filtered_seq += ctg_stats[hdr][0]
    logger.debug("Filtered %d contigs of total length %d",
                 filtered_num, filtered_seq)

    fp.write_fasta_dict(good_fasta, contigs_out)
    with open(stats_out, "w") as f:
        f.write("#seq_name\tlength\tcoverage\n")
        for ctg_id in good_fasta:
            f.write("{0}\t{1}\t{2}\n".format(ctg_id,
                    ctg_stats[ctg_id][0], ctg_stats[ctg_id][1]))


def _run_polish_bin(bubbles_in, subs_matrix, hopo_matrix,
                    consensus_out, num_threads, output_progress, use_hopo):
    """
    Invokes polishing binary
    """
    cmdline = [POLISH_BIN, "polisher", "--bubbles", bubbles_in, "--subs-mat", subs_matrix,
               "--hopo-mat", hopo_matrix, "--out", consensus_out,
               "--threads", str(num_threads)]
    if not output_progress:
        cmdline.append("--quiet")

    if use_hopo:
        cmdline.append("--enable-hopo")

    try:
        subprocess.check_call(cmdline)
    except subprocess.CalledProcessError as e:
        if e.returncode == -9:
            logger.error("Looks like the system ran out of memory")
        raise PolishException(str(e))
    except OSError as e:
        raise PolishException(str(e))


def _compose_sequence(consensus_file):
    """
    Concatenates bubbles consensuses into genome
    """
    consensuses = defaultdict(list)
    coverage = defaultdict(list)
    with open(consensus_file, "r") as f:
        header = True
        for line in f:
            if header:
                tokens = line.strip().split(" ")
                if len(tokens) != 4:
                    raise Exception("Bubble format error")

                ctg_id = tokens[0][1:]
                ctg_pos = int(tokens[1])
                coverage = int(tokens[2])
                ctg_sub_pos = int(tokens[3])
            else:
                consensuses[ctg_id].append((ctg_pos, ctg_sub_pos, coverage, line.strip()))
            header = not header

    polished_fasta = {}
    polished_stats = {}
    polished_coverages = {}
    for ctg_id, seqs in iteritems(consensuses):
        seqs.sort(key=lambda p: (p[0], p[1]))
        sorted_seqs = [p[3] for p in seqs]
        bubble_coverages = [(len(p[3]), p[2]) for p in seqs]
        concat_seq = "".join(sorted_seqs)
        #mean_coverage = sum(coverage[ctg_id]) / len(coverage[ctg_id])
        polished_fasta[ctg_id] = concat_seq
        polished_stats[ctg_id] = len(concat_seq)
        polished_coverages[ctg_id] = bubble_coverages

    return polished_fasta, polished_stats, polished_coverages
