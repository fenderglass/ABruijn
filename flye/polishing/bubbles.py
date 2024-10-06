#(c) 2016 by Authors
#This file is a part of ABruijn program.
#Released under the BSD license (see LICENSE file)

"""
Separates alignment into small bubbles for further correction
"""

import logging
from bisect import bisect
from collections import defaultdict
from queue import Queue

import multiprocessing
import traceback

import flye.utils.fasta_parser as fp
import flye.config.py_cfg as cfg
from flye.polishing.alignment import shift_gaps, get_uniform_alignments
from flye.utils.sam_parser import SynchronizedSamReader, SynchonizedChunkManager
from flye.utils.utils import process_in_parallel, get_median


logger = logging.getLogger()


class ProfileInfo(object):
    __slots__ = ("nucl", "insertions", "propagated_ins", "num_deletions",
                 "num_missmatch", "coverage")

    def __init__(self):
        self.nucl = ""
        #self.num_inserts = 0
        self.propagated_ins = 0
        self.insertions = defaultdict(str)
        self.num_deletions = 0
        self.num_missmatch = 0
        self.coverage = 0


class Bubble(object):
    __slots__ = ("contig_id", "position", "sub_position", "branches", "consensus")

    def __init__(self, contig_id, position):
        self.contig_id = contig_id
        self.position = position
        self.sub_position = 0
        self.branches = []
        self.consensus = ""


def _thread_worker(aln_reader, chunk_feeder, contigs_info, err_mode,
                   results_queue, error_queue, bubbles_file,
                   bubbles_file_lock):
    """
    Will run in parallel
    """
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
            ref_seq = aln_reader.get_region_sequence(ctg_region.ctg_id, ctg_region.start,
                                                     ctg_region.end)

            #since we are working with contig chunks, tranform alignment coorinates
            ctg_aln = aln_reader.trim_and_transpose(ctg_aln, ctg_region.start, ctg_region.end)
            ctg_aln, mean_cov = get_uniform_alignments(ctg_aln)

            profile, aln_errors = _compute_profile(ctg_aln, ref_seq)
            partition, num_long_bubbles = _get_partition(profile, err_mode)
            ctg_bubbles = _get_bubble_seqs(ctg_aln, profile, partition, ctg_id)

            ##
            coverage_cap = 0.9 * cfg.vals["max_read_coverage"]
            if mean_cov > coverage_cap:
                mean_cov = aln_reader.get_median_depth(ctg_region.ctg_id, ctg_region.start,
                                                       ctg_region.end)
            ##

            ctg_bubbles, num_empty = _postprocess_bubbles(ctg_bubbles)
            ctg_bubbles, num_long_branch = _split_long_bubbles(ctg_bubbles)

            #transform coordinates back
            for b in ctg_bubbles:
                b.position += ctg_region.start

            if bubbles_file_lock:
                bubbles_file_lock.acquire()

            _output_bubbles(ctg_bubbles, open(bubbles_file, "a"))
            results_queue.put((ctg_id, len(ctg_bubbles), num_long_bubbles,
                               num_empty, num_long_branch, aln_errors,
                               mean_cov))
            
            if bubbles_file_lock:
                bubbles_file_lock.release()

            del profile
            del ctg_bubbles

    except Exception as e:
        logger.error("Thread exception")
        logger.error(traceback.format_exc())
        error_queue.put(e)


def make_bubbles(alignment_path, contigs_info, contigs_path,
                 err_mode, num_proc, bubbles_out):
    """
    The main function: takes an alignment and returns bubbles
    """
    CHUNK_SIZE = 1000000

    contigs_fasta = fp.read_sequence_dict(contigs_path)
    manager = None if num_proc == 1 else multiprocessing.Manager()
    aln_reader = SynchronizedSamReader(alignment_path, contigs_fasta, manager,
                                       cfg.vals["max_read_coverage"], use_secondary=True)
    chunk_feeder = SynchonizedChunkManager(contigs_fasta, manager, chunk_size=CHUNK_SIZE)

    if manager:
        results_queue = manager.Queue()
        error_queue = manager.Queue()
        bubbles_out_lock = multiprocessing.Lock()

        process_in_parallel(_thread_worker, (aln_reader, chunk_feeder, contigs_info, err_mode,
                         results_queue, error_queue, bubbles_out, bubbles_out_lock), num_proc)
    else:
        results_queue = Queue()
        error_queue = Queue()
        bubbles_out_lock = None

        _thread_worker(aln_reader, chunk_feeder, contigs_info, err_mode,
                results_queue, error_queue, bubbles_out, bubbles_out_lock)
        
    if not error_queue.empty():
        raise error_queue.get()

    #logging
    total_bubbles = 0
    total_long_bubbles = 0
    total_long_branches = 0
    total_empty = 0
    total_aln_errors = []
    coverage_stats = defaultdict(list)

    while not results_queue.empty():
        (ctg_id, num_bubbles, num_long_bubbles,
            num_empty, num_long_branch,
            aln_errors, mean_coverage) = results_queue.get()
        total_long_bubbles += num_long_bubbles
        total_long_branches += num_long_branch
        total_empty += num_empty
        total_aln_errors.extend(aln_errors)
        total_bubbles += num_bubbles
        coverage_stats[ctg_id].append(mean_coverage)

    for ctg in coverage_stats:
        coverage_stats[ctg] = int(sum(coverage_stats[ctg]) / len(coverage_stats[ctg]))

    mean_aln_error = sum(total_aln_errors) / (len(total_aln_errors) + 1)
    logger.debug("Generated %d bubbles", total_bubbles)
    logger.debug("Split %d long bubbles", total_long_bubbles)
    logger.debug("Skipped %d empty bubbles", total_empty)
    logger.debug("Skipped %d bubbles with long branches", total_long_branches)
    ###

    return coverage_stats, mean_aln_error


def _output_bubbles(bubbles, out_stream):
    """
    Outputs list of bubbles into file
    """
    for bubble in bubbles:
        if len(bubble.branches) == 0:
            raise Exception("No branches in a bubble")
        out_stream.write(">{0} {1} {2} {3}\n".format(bubble.contig_id,
                                                     bubble.position,
                                                     len(bubble.branches),
                                                     bubble.sub_position))
        out_stream.write(bubble.consensus + "\n")
        for branch_id, branch in enumerate(bubble.branches):
            out_stream.write(">{0}\n".format(branch_id))
            out_stream.write(branch + "\n")

    out_stream.flush()


def _split_long_bubbles(bubbles):
    MAX_BUBBLE = cfg.vals["max_bubble_length"]
    #MAX_BUBBLE = 50
    #MAX_BRANCH = MAX_BUBBLE * 1.5
    new_bubbles = []
    long_branches = 0

    for bubble in bubbles:
        median_branch = sorted(bubble.branches, key=len)[len(bubble.branches) // 2]
        num_chunks = len(median_branch) // MAX_BUBBLE
        #if len(median_branch) > MAX_BRANCH:
        if num_chunks > 1:
            #logger.debug("Splitting: pos:{0} len:{1}".format(bubble.position, len(median_branch)))
            long_branches += 1

            for part_num in range(num_chunks):
                new_branches = []
                for b in bubble.branches:
                    chunk_len = len(b) // num_chunks
                    start = part_num * chunk_len
                    end = (part_num + 1) * chunk_len if part_num != num_chunks - 1 else len(b)
                    new_branches.append(b[start:end])

                new_bubbles.append(Bubble(bubble.contig_id, bubble.position))
                new_bubbles[-1].consensus = new_branches[0]
                new_bubbles[-1].branches = new_branches
                new_bubbles[-1].sub_position = part_num

        else:
            new_bubbles.append(bubble)

    return new_bubbles, long_branches


def _postprocess_bubbles(bubbles):
    MAX_BUBBLE = cfg.vals["max_bubble_length"]
    MAX_BRANCHES = cfg.vals["max_bubble_branches"]

    new_bubbles = []
    empty_bubbles = 0
    for bubble in bubbles:
        if len(bubble.branches) == 0:
            empty_bubbles += 1
            continue

        median_branch = sorted(bubble.branches, key=len)[len(bubble.branches) // 2]
        if len(median_branch) == 0:
            #logger.debug("Median branch with zero length: {0}".format(bubble.position))
            empty_bubbles += 1
            continue

        #only take branches that are not significantly differ in length from the median
        new_branches = []
        for branch in bubble.branches:
            incons_rate = abs(len(branch) - len(median_branch)) / len(median_branch)
            if incons_rate < 0.5 and len(branch) > 0:
                new_branches.append(branch)

        #checking again (since we might have tossed some branches)
        if len(new_branches) == 0:
            empty_bubbles += 1
            continue

        #if bubble consensus has very different length from all the branchs, replace
        #consensus with the median branch instead
        if abs(len(median_branch) - len(bubble.consensus)) > len(median_branch) // 2:
            bubble.consensus = median_branch

        #finally, keep only MAX_BRANCHES
        if len(new_branches) > MAX_BRANCHES:
            left = len(new_branches) // 2 - MAX_BRANCHES // 2
            right = left + MAX_BRANCHES
            new_branches = sorted(new_branches, key=len)[left:right]

        new_bubbles.append(Bubble(bubble.contig_id, bubble.position))
        new_bubbles[-1].consensus = bubble.consensus
        new_bubbles[-1].branches = new_branches

    return new_bubbles, empty_bubbles


def _is_solid_kmer(profile, position, err_mode):
    """
    Checks if the kmer at given position is solid
    """
    MISSMATCH_RATE = cfg.vals["err_modes"][err_mode]["solid_missmatch"]
    INS_RATE = cfg.vals["err_modes"][err_mode]["solid_indel"]
    SOLID_LEN = cfg.vals["solid_kmer_length"]

    for i in range(position, position + SOLID_LEN):
        if profile[i].coverage == 0:
            return False
        local_missmatch = (profile[i].num_missmatch +
                           profile[i].num_deletions) / profile[i].coverage
        #local_ins = len(profile[i].insertions) / profile[i].coverage
        local_ins = profile[i].propagated_ins / profile[i].coverage
        if local_missmatch > MISSMATCH_RATE or local_ins > INS_RATE:
            return False
    return True


def _is_simple_kmer(profile, position):
    """
    Checks if the kmer with center at the given position is simple
    """
    SIMPLE_LEN = cfg.vals["simple_kmer_length"]

    extended_len = SIMPLE_LEN * 2
    nucl_str = [p.nucl for p in profile[position - extended_len // 2 :
                                        position + extended_len // 2]]

    #single nucleotide homopolymers
    for i in range(extended_len // 2 - SIMPLE_LEN // 2,
                   extended_len // 2 + SIMPLE_LEN // 2 - 1):
        if nucl_str[i] == nucl_str[i + 1]:
            return False

    #dinucleotide homopolymers
    for shift in [0, 1]:
        for i in range(SIMPLE_LEN - shift - 1):
            pos = extended_len // 2 - SIMPLE_LEN + shift + i * 2
            if (nucl_str[pos : pos + 2] == nucl_str[pos + 2 : pos + 4]):
                return False

    #trinucleotide homopolymers
    #for shift in [0, 1, 2]:
    #    for i in xrange(SIMPLE_LEN - shift - 1):
    #        pos = shift + i * 3
    #        if (nucl_str[pos : pos + 3] == nucl_str[pos + 3 : pos + 6]):
    #            #logger.debug("tri" + "".join(nucl_str))
    #            return False

    return True


def _compute_profile(alignment, ref_sequence):
    """
    Computes alignment profile
    """
    if len(alignment) == 0:
        raise Exception("No alignmemnts!")
    genome_len = alignment[0].trg_len

    #max_aln_err = cfg.vals["err_modes"][platform]["max_aln_error"]
    min_aln_len = min(cfg.vals["min_polish_aln_len"], genome_len // 2)
    aln_errors = []
    #filtered = 0
    profile = [ProfileInfo() for _ in range(genome_len)]

    for i in range(genome_len):
        profile[i].nucl = ref_sequence[i]

    for aln in alignment:
        #if aln.err_rate > max_aln_err or len(aln.qry_seq) < min_aln_len:
        if len(aln.qry_seq) < min_aln_len:
            #filtered += 1
            continue

        aln_errors.append(aln.err_rate)

        qry_seq = shift_gaps(aln.trg_seq, aln.qry_seq)
        trg_seq = shift_gaps(qry_seq, aln.trg_seq)

        trg_pos = aln.trg_start
        for trg_nuc, qry_nuc in zip(trg_seq, qry_seq):
            if trg_nuc == "-":
                trg_pos -= 1
            #if trg_pos >= genome_len:
            #    trg_pos -= genome_len

            prof_elem = profile[trg_pos]
            if trg_nuc == "-":
                prof_elem.insertions[aln.qry_id] += qry_nuc
                #prof_elem.num_inserts += 1
            else:
                #prof_elem.nucl = trg_nuc
                prof_elem.coverage += 1

                if qry_nuc == "-":
                    prof_elem.num_deletions += 1
                elif trg_nuc != qry_nuc:
                    prof_elem.num_missmatch += 1

            trg_pos += 1

    for i in range(genome_len):
        for ins_read, ins_str in profile[i].insertions.items():
            profile[i].propagated_ins += 1
            span = len(ins_str)
            for j in range(max(0, i - span), i):
                profile[j].propagated_ins += 1
            for j in range(i + 1, min(i + span + 1, genome_len)):
                profile[j].propagated_ins += 1
            

    #logger.debug("Filtered: {0} out of {1}".format(filtered, len(alignment)))
    return profile, aln_errors


def _get_partition(profile, err_mode):
    """
    Partitions genome into sub-alignments at solid regions / simple kmers
    """
    #logger.debug("Partitioning genome")
    SOLID_LEN = cfg.vals["solid_kmer_length"]
    SIMPLE_LEN = cfg.vals["simple_kmer_length"]
    MAX_BUBBLE = cfg.vals["max_bubble_length"]

    solid_flags = [False for _ in range(len(profile))]
    prof_pos = 0
    while prof_pos < len(profile) - SOLID_LEN:
        if _is_solid_kmer(profile, prof_pos, err_mode):
            for i in range(prof_pos, prof_pos + SOLID_LEN):
                solid_flags[i] = True
            prof_pos += SOLID_LEN
        else:
            prof_pos += 1

    partition = []
    prev_partition = SOLID_LEN

    long_bubbles = 0
    prof_pos = SOLID_LEN
    while prof_pos < len(profile) - SOLID_LEN:
        cur_partition = prof_pos + SIMPLE_LEN // 2
        landmark = (all(solid_flags[prof_pos : prof_pos + SIMPLE_LEN]) and
                    _is_simple_kmer(profile, cur_partition))

        if prof_pos - prev_partition > MAX_BUBBLE:
            long_bubbles += 1

        if landmark or prof_pos - prev_partition > MAX_BUBBLE:
            partition.append(cur_partition)
            prev_partition = cur_partition
            prof_pos += SOLID_LEN
        else:
            prof_pos += 1

    #logger.debug("Partitioned into {0} segments".format(len(partition) + 1))
    #logger.debug("Long bubbles: {0}".format(long_bubbles))

    return partition, long_bubbles


def _get_bubble_seqs(alignment, profile, partition, contig_id):
    """
    Given genome landmarks, forms bubble sequences
    """
    if not partition or not alignment:
        return []

    ctg_len = alignment[0].trg_len

    bubbles = []
    ext_partition = [0] + partition + [ctg_len]
    for p_left, p_right in zip(ext_partition[:-1], ext_partition[1:]):
        bubbles.append(Bubble(contig_id, p_left))
        consensus = [p.nucl for p in profile[p_left : p_right]]
        bubbles[-1].consensus = "".join(consensus)

    for aln in alignment:
        bubble_id = bisect(ext_partition, aln.trg_start) - 1
        next_bubble_start = ext_partition[bubble_id + 1]
        chromosome_end = aln.trg_end >= ext_partition[-1]

        incomplete_segment = aln.trg_start > ext_partition[bubble_id]
        trg_pos = aln.trg_start
        branch_start = 0
        for i, trg_nuc in enumerate(aln.trg_seq):
            if trg_nuc == "-":
                continue
            #if trg_pos >= contig_info.length:
                #trg_pos -= contig_info.length

            if trg_pos >= next_bubble_start:
                if not incomplete_segment:
                    branch_seq = fp.to_acgt(aln.qry_seq[branch_start : i].replace("-", ""))
                    bubbles[bubble_id].branches.append(branch_seq)

                incomplete_segment = False
                bubble_id = bisect(ext_partition, trg_pos) - 1
                next_bubble_start = ext_partition[bubble_id + 1]
                branch_start = i

            trg_pos += 1

        if chromosome_end:
            branch_seq = fp.to_acgt(aln.qry_seq[branch_start:].replace("-", ""))
            bubbles[-1].branches.append(branch_seq)

    return bubbles
