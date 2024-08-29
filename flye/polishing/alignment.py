#(c) 2016 by Authors
#This file is a part of ABruijn program.
#Released under the BSD license (see LICENSE file)

"""
Runs Minimap2 and parses its output
"""

from __future__ import absolute_import
from __future__ import division
import os
from collections import namedtuple
import subprocess
import logging
import datetime
from copy import copy

import flye.utils.fasta_parser as fp
from flye.utils.utils import which, get_median
from flye.utils.sam_parser import AlignmentException
from flye.six import iteritems
from flye.six.moves import range


logger = logging.getLogger()
MINIMAP_BIN = "flye-minimap2"
SAMTOOLS_BIN = "flye-samtools"

ContigInfo = namedtuple("ContigInfo", ["id", "length", "type"])


def check_binaries():
    if not which(MINIMAP_BIN):
        raise AlignmentException("minimap2 is not installed")
    if not which(SAMTOOLS_BIN):
        raise AlignmentException("samtools is not installed")
    if not which("sort"):
        raise AlignmentException("UNIX sort utility is not available")


def make_alignment(reference_file, reads_file, num_proc,
                   platform, read_type, out_alignment):
    mode = platform + "-" + read_type
    _run_minimap(reference_file, reads_file, num_proc, mode, out_alignment)


def get_contigs_info(contigs_file):
    contigs_info = {}
    contigs_fasta = fp.read_sequence_dict(contigs_file)
    for ctg_id, ctg_seq in iteritems(contigs_fasta):
        contig_type = ctg_id.split("_")[0]
        contigs_info[ctg_id] = ContigInfo(ctg_id, len(ctg_seq),
                                          contig_type)

    return contigs_info


def shift_gaps(seq_trg, seq_qry):
    """
    Shifts all ambigious query gaps to the right
    """
    lst_trg, lst_qry = list("$" + seq_trg + "$"), list("$" + seq_qry + "$")
    is_gap = False
    gap_start = 0
    for i in range(len(lst_trg)):
        if is_gap and lst_qry[i] != "-":
            is_gap = False
            swap_left = gap_start - 1
            swap_right = i - 1

            while (swap_left > 0 and swap_right >= gap_start and
                   lst_qry[swap_left] == lst_trg[swap_right]):
                lst_qry[swap_left], lst_qry[swap_right] = \
                            lst_qry[swap_right], lst_qry[swap_left]
                swap_left -= 1
                swap_right -= 1

        if not is_gap and lst_qry[i] == "-":
            is_gap = True
            gap_start = i

    return "".join(lst_qry[1 : -1])


def get_uniform_alignments(alignments):
    """
    Leaves top alignments for each position within contig
    assuming uniform coverage distribution
    """
    if not alignments:
        return []

    WINDOW = 100
    MIN_COV = 20
    GOOD_RATE = 0.66
    MIN_QV = 20

    def is_reliable(aln):
        return not aln.is_secondary and aln.map_qv >= MIN_QV
        #return not aln.is_secondary and not aln.is_supplementary and aln.map_qv >= MIN_QV

    seq_len = alignments[0].trg_len
    ctg_id = alignments[0].trg_id

    #split contig into windows, get median read coverage over all windows and
    #determine the quality threshold cutoffs for each window
    wnd_primary_cov = [0 for _ in range(seq_len // WINDOW + 1)]
    wnd_all_cov = [0 for _ in range(seq_len // WINDOW + 1)]

    for aln in alignments:
        for i in range(aln.trg_start // WINDOW, aln.trg_end // WINDOW + 1):
            if is_reliable(aln):
                wnd_primary_cov[i] += 1
            wnd_all_cov[i] += 1

    cov_threshold = max(int(get_median(wnd_primary_cov)), MIN_COV)
    orig_primary_cov = copy(wnd_primary_cov)

    selected_alignments = []
    original_sequence = 0
    primary_sequence = 0
    secondary_sequence = 0
    primary_aln = 0
    secondary_aln = 0

    def _aln_score(aln):
        wnd_good = 0
        wnd_bad = 0
        for i in range(aln.trg_start // WINDOW, aln.trg_end // WINDOW + 1):
            if wnd_primary_cov[i] < cov_threshold:
                wnd_good += 1
            else:
                wnd_bad += 1
        return wnd_good, wnd_bad

    sec_aln_scores = {}
    for aln in alignments:
        original_sequence += aln.trg_end - aln.trg_start

        #always keep primary alignments, regardless of local coverage
        if is_reliable(aln):
            primary_sequence += aln.trg_end - aln.trg_start
            primary_aln += 1
            selected_alignments.append(aln)

        #if alignment is secondary, count how many windows it helps to improve
        else:
            wnd_good, wnd_bad = _aln_score(aln)
            sec_aln_scores[aln.qry_id] = (wnd_good, wnd_bad, aln)

    #now, greedily add secondaty alignments, until they add useful coverage
    _score_fun = lambda x: (sec_aln_scores[x][0] - 2 * sec_aln_scores[x][1],
                            sec_aln_scores[x][2].trg_end - sec_aln_scores[x][2].trg_start)
    sorted_sec_aln = [x for x in sorted(sec_aln_scores, reverse=True, key=_score_fun)]
    for aln_id in sorted_sec_aln:
        aln = sec_aln_scores[aln_id][2]
        #recompute scores
        wnd_good, wnd_bad = _aln_score(aln)
        to_take = wnd_good / (wnd_good + wnd_bad) > GOOD_RATE
        if to_take:
            selected_alignments.append(aln)
            secondary_aln += 1
            secondary_sequence += aln.trg_end - aln.trg_start
            for i in range(aln.trg_start // WINDOW, aln.trg_end // WINDOW + 1):
                wnd_primary_cov[i] += 1

        #logger.debug("\tSec score: {} {} {} {}".format(aln_id, wnd_good, wnd_bad, to_take))

    #logger.debug("Seq: {0} pri_cov: {1} all_cov: {2}".format(ctg_id, get_median(orig_primary_cov),
    #                                                         get_median(wnd_all_cov)) + "\n" +
    #            "\tOriginal seq: {0}, reads: {1}".format(original_sequence, len(alignments)) + "\n" +
    #            "\tPrimary seq: {0}, reads: {1}".format(primary_sequence, primary_aln) + "\n" +
    #            "\tSecondary seq: {0}, reads: {1}".format(secondary_sequence, secondary_aln) + "\n" +
    #            "\tSelected size: {0}, median coverage: {1}".format(len(selected_alignments), get_median(wnd_primary_cov)))

    median_cov = get_median(wnd_primary_cov)
    #mean_cov = sum(wnd_primary_cov) / (len(wnd_primary_cov) + 1)

    return selected_alignments, median_cov


def split_into_chunks(fasta_in, chunk_size, fasta_out):
    out_dict = {}
    for header, seq in fp.stream_sequence(fasta_in):
        #print len(seq)
        for i in range(0, max(len(seq) // chunk_size, 1)):
            chunk_hdr = "{0}$chunk_{1}".format(header, i)
            start = i * chunk_size
            end = (i + 1) * chunk_size
            if len(seq) - end < chunk_size:
                end = len(seq)

            #print(start, end)
            out_dict[chunk_hdr] = seq[start : end]

    fp.write_fasta_dict(out_dict, fasta_out)


def merge_chunks(fasta_in, fold_function=lambda l: "".join(l)):
    """
    Merges sequence chunks. Chunk names are in format `orig_name$chunk_id`.
    Each chunk is as dictionary entry. Value type is arbitrary and
    one can supply a custom fold function
    """
    def name_split(h):
        orig_hdr, chunk_id = h.rsplit("$", 1)
        return orig_hdr, int(chunk_id.rsplit("_", 1)[1])

    out_dict = {}
    cur_seq = []
    cur_contig = None
    for hdr in sorted(fasta_in, key=name_split):
        orig_name, dummy_chunk_id = name_split(hdr)
        if orig_name != cur_contig:
            if cur_contig != None:
                out_dict[cur_contig] = fold_function(cur_seq)
            cur_seq = []
            cur_contig = orig_name
        cur_seq.append(fasta_in[hdr])

    if cur_seq:
        out_dict[cur_contig] = fold_function(cur_seq)

    return out_dict


def _run_minimap(reference_file, reads_files, num_proc, reads_type, out_file):
    #SAM_HEADER = "\'@PG|@HD|@SQ|@RG|@CO\'"
    work_dir = os.path.dirname(out_file)
    stderr_file = os.path.join(work_dir, "minimap.stderr")
    SORT_THREADS = "4"
    SORT_MEM = "4G" if os.path.getsize(reference_file) > 100 * 1024 * 1024 else "1G"
    BATCH = "5G" if os.path.getsize(reference_file) > 100 * 1024 * 1024 else "1G"

    mode = None
    extra_args = []
    if reads_type in ["pacbio-raw", "pacbio-corrected", "pacbio-subasm"]:
        mode = "map-pb"
    if reads_type == "pacbio-hifi":
        mode = "map-hifi"
    elif reads_type in ["nano-raw", "nano-corrected"]:
        mode = "map-ont"
    elif reads_type == "nano-nano_hq":
        mode = "lr:hq"
        extra_args = ["-k", "17"]

    cmdline = [MINIMAP_BIN, "'" + reference_file + "'"]
    cmdline.extend(["'" + read_file + "'" for read_file in reads_files])
    cmdline.extend(["-x", mode, "-t", str(num_proc)])
    cmdline.extend(extra_args)

    #Produces gzipped SAM sorted by reference name. Since it's not sorted by
    #read name anymore, it's important that all reads have SEQ.
    #is sam_output not set, produces PAF alignment
    #a = SAM output, p = min primary-to-seconday score
    #N = max secondary alignments
    #--sam-hit-only = don't output unmapped reads
    #--secondary-seq = custom option to output SEQ for seqcondary alignment with hard clipping
    #-L: move CIGAR strings for ultra-long reads to the separate tag
    #-Q don't output fastq quality
    tmp_prefix = os.path.join(os.path.dirname(out_file),
                              "sort_" + datetime.datetime.now().strftime("%y%m%d_%H%M%S"))
    cmdline.extend(["-a", "-p", "0.5", "-N", "10", "--sam-hit-only", "-L", "-K", BATCH,
                    "-z", "1000", "-Q", "--secondary-seq", "-I", "64G"])
    cmdline.extend(["|", SAMTOOLS_BIN, "view", "-T", "'" + reference_file + "'", "-u", "-"])
    cmdline.extend(["|", SAMTOOLS_BIN, "sort", "-T", "'" + tmp_prefix + "'", "-O", "bam",
                    "-@", SORT_THREADS, "-l", "1", "-m", SORT_MEM])
    cmdline.extend(["-o", "'" + out_file + "'"])

    #logger.debug("Running: " + " ".join(cmdline))
    try:
        devnull = open(os.devnull, "wb")
        subprocess.check_call(["/bin/bash", "-c",
                              "set -eo pipefail; " + " ".join(cmdline)],
                              stderr=open(stderr_file, "w"),
                              stdout=open(os.devnull, "w"))
        subprocess.check_call(SAMTOOLS_BIN + " index -@ 4 " + "'" + out_file + "'", shell=True)
        #os.remove(stderr_file)

    except (subprocess.CalledProcessError, OSError) as e:
        logger.error("Error running minimap2, terminating. See the alignment error log "
                     " for details: " + stderr_file)
        logger.error("Cmd: " + " ".join(cmdline))
        raise AlignmentException(str(e))
