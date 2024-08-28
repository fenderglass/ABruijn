#(c) 2016 by Authors
#This file is a part of ABruijn program.
#Released under the BSD license (see LICENSE file)

"""
Quick and dirty alignment consensus
"""

import logging
from collections import defaultdict
from flye.six import itervalues

import multiprocessing
import traceback

from flye.polishing.alignment import shift_gaps, get_uniform_alignments
from flye.utils.sam_parser import SynchronizedSamReader, SynchonizedChunkManager
import flye.config.py_cfg as cfg
import flye.utils.fasta_parser as fp
from flye.utils.utils import process_in_parallel

logger = logging.getLogger()


class Profile(object):
    __slots__ = ("insertions", "matches", "nucl")

    def __init__(self):
        self.insertions = defaultdict(str)
        self.matches = defaultdict(int)
        self.nucl = "-"


def _thread_worker(aln_reader, chunk_feeder, platform, results_queue, error_queue):
    try:
        while True:
            ctg_region = chunk_feeder.get_chunk()
            if ctg_region is None:
                break
            ctg_aln = aln_reader.get_alignments(ctg_region.ctg_id, ctg_region.start,
                                                ctg_region.end)
            ctg_id = ctg_region.ctg_id
            if len(ctg_aln) == 0:
                continue

            ctg_aln = aln_reader.trim_and_transpose(ctg_aln, ctg_region.start, ctg_region.end)
            ctg_aln, _mean_cov = get_uniform_alignments(ctg_aln)

            profile, aln_errors = _contig_profile(ctg_aln, platform)
            sequence = _flatten_profile(profile)
            results_queue.put((ctg_id, ctg_region.start, sequence, aln_errors))

    except Exception as e:
        logger.error("Thread exception")
        logger.error(traceback.format_exc())
        error_queue.put(e)


def get_consensus(alignment_path, contigs_path, contigs_info, num_proc,
                  platform):
    """
    Main function
    """

    CHUNK_SIZE = 1000000
    contigs_fasta = fp.read_sequence_dict(contigs_path)
    mp_manager = multiprocessing.Manager()
    aln_reader = SynchronizedSamReader(alignment_path, contigs_fasta, mp_manager,
                                       max_coverage=cfg.vals["max_read_coverage"],
                                       use_secondary=True)
    chunk_feeder = SynchonizedChunkManager(contigs_fasta, mp_manager, CHUNK_SIZE)

    #manager = multiprocessing.Manager()
    results_queue = mp_manager.Queue()
    error_queue = mp_manager.Queue()

    process_in_parallel(_thread_worker, (aln_reader, chunk_feeder,
                            platform, results_queue, error_queue), num_proc)

    if not error_queue.empty():
        raise error_queue.get()

    chunk_consensus = defaultdict(list)
    total_aln_errors = []
    while not results_queue.empty():
        ctg_id, region_start, ctg_seq, aln_errors = results_queue.get()
        total_aln_errors.extend(aln_errors)
        if len(ctg_seq) > 0:
            chunk_consensus[ctg_id].append((region_start, ctg_seq))

    out_fasta = {}
    for ctg in chunk_consensus:
        sorted_chunks = [x[1] for x in sorted(chunk_consensus[ctg], key=lambda p: p[0])]
        out_fasta[ctg] = "".join(sorted_chunks)

    mean_aln_error = sum(total_aln_errors) / (len(total_aln_errors) + 1)
    logger.info("Alignment error rate: %f", mean_aln_error)

    return out_fasta


def _contig_profile(alignment, platform):
    """
    Computes alignment profile
    """

    if not alignment:
        return []

    genome_len = alignment[0].trg_len

    aln_errors = []
    profile = [Profile() for _ in range(genome_len)]
    #max_aln_err = cfg.vals["err_modes"][platform]["max_aln_error"]
    for aln in alignment:
        #if aln.err_rate > max_aln_err: continue
        aln_errors.append(aln.err_rate)

        #after gap shifting it is possible that
        #two gaps are aligned against each other
        qry_seq = shift_gaps(aln.trg_seq, aln.qry_seq)
        trg_seq = shift_gaps(qry_seq, aln.trg_seq)

        trg_pos = aln.trg_start
        for trg_nuc, qry_nuc in zip(trg_seq, qry_seq):
            if trg_nuc == "-":
                trg_pos -= 1
            if trg_pos >= genome_len:
                trg_pos -= genome_len

            #total += 1
            prof_elem = profile[trg_pos]
            if trg_nuc == "-" and qry_nuc != "-":
                prof_elem.insertions[aln.qry_id] += qry_nuc
            else:
                prof_elem.nucl = trg_nuc
                prof_elem.matches[qry_nuc] += 1

            trg_pos += 1

    #print "len", genome_len, "median coverage", cov_threshold
    #print "total bases: ", total, "discarded bases: ", discarded
    #print "filtered", float(discarded) / total
    #print ""

    return profile, aln_errors


def _flatten_profile(profile):
    growing_seq = []
    ins_group = defaultdict(int)

    for elem in profile:
        pos_matches = elem.matches
        pos_insertions = elem.insertions
        pos_nucl = elem.nucl

        ins_group.clear()
        for ins_str in itervalues(pos_insertions):
            ins_group[ins_str] += 1

        match_and_del_num = sum(itervalues(pos_matches))
        del_num = pos_matches["-"]
        num_ins = len(pos_insertions)

        max_match = pos_nucl
        if len(pos_matches):
            max_match = max(sorted(pos_matches), key=pos_matches.get)
        max_insert = None
        if ins_group:
            max_insert = max(sorted(ins_group), key=ins_group.get)

        is_deletion = max_match == "-" or del_num > match_and_del_num // 3
        is_insertion = max_insert and max_insert != "-" and num_ins > match_and_del_num // 3

        if not is_deletion:
            growing_seq.append(max_match)
        if is_insertion:
            growing_seq.append(max_insert)

    return "".join(growing_seq)
