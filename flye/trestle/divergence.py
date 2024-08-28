#(c) 2016-2018 by Authors
#This file is a part of Flye program.
#Released under the BSD license (see LICENSE file)

"""
Created on Wed Jan  4 03:50:31 2017

@author: jeffrey_yuan
"""

import logging
from collections import defaultdict

import multiprocessing
import os.path

from flye.polishing.alignment import shift_gaps
from flye.utils.sam_parser import SynchronizedSamReader, SynchonizedChunkManager
import flye.utils.fasta_parser as fp
from flye.utils.utils import process_in_parallel
import flye.config.py_cfg as config

logger = logging.getLogger()

class Profile(object):
    __slots__ = ("insertions", "matches", "nucl")

    def __init__(self):
        self.insertions = defaultdict(int)
        self.matches = defaultdict(int)
        self.nucl = "-"


def _thread_worker(aln_reader, chunk_feeder, contigs_info, platform, results_queue,
                   error_queue):
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

            profile, aln_errors = _contig_profile(ctg_aln, platform,
                                                  contigs_info[ctg_id].length)
            #sequence = _flatten_profile(profile)
            results_queue.put((ctg_id, profile, aln_errors))

    except Exception as e:
        error_queue.put(e)


def _contig_profile(alignment, platform, genome_len):
    """
    Computes alignment profile
    """
    #max_aln_err = config.vals["err_modes"][platform]["max_aln_error"]
    aln_errors = []
    profile = [Profile() for _ in range(genome_len)]
    for aln in alignment:
        #if aln.err_rate > max_aln_err: continue
        aln_errors.append(aln.err_rate)

        qry_seq = shift_gaps(aln.trg_seq, aln.qry_seq)
        trg_seq = shift_gaps(qry_seq, aln.trg_seq)
        #qry_seq = aln.qry_seq
        #trg_seq = aln.trg_seq

        trg_pos = aln.trg_start
        for trg_nuc, qry_nuc in zip(trg_seq, qry_seq):
            if trg_nuc == "-":
                trg_pos -= 1
            if trg_pos >= genome_len:
                trg_pos -= genome_len

            prof_elem = profile[trg_pos]
            if trg_nuc == "-":
                prof_elem.insertions[qry_nuc] += 1
            else:
                prof_elem.nucl = trg_nuc
                prof_elem.matches[qry_nuc] += 1

            trg_pos += 1

    return profile, aln_errors


def _count_freqs(elem):
    matches = elem.matches
    insertions = elem.insertions
    nucl = elem.nucl
    del_key = "-"

    coverage = sum(matches.values())

    mat_ct = 0
    if nucl in matches:
        mat_ct = matches[nucl]

    subs = {key:matches[key] for key in matches
                                if (key != nucl and key != del_key)}
    max_sub_ct = 0
    max_sub_base = ""
    if subs:
        max_sub_base = max(subs, key=subs.get)
        max_sub_ct = subs[max_sub_base]

    del_ct = 0
    if del_key in matches:
        del_ct = matches[del_key]

    max_ins_key = "^"
    max_ins_ct = 0
    if insertions:
        max_ins_base = max(insertions, key=insertions.get)
        max_ins_ct = insertions[max_ins_base]
        max_ins_key = "^{0}".format(max_ins_base)

    return {'cov':coverage, 'mat_base':nucl, 'mat_ct':mat_ct,
                            'sub_base':max_sub_base, 'sub_ct':max_sub_ct,
                            'del_base':del_key, 'del_ct':del_ct,
                            'ins_base':max_ins_key, 'ins_ct':max_ins_ct}


def _call_position(ind, counts, pos, sub_thresh, del_thresh, ins_thresh):
    over_thresh = False
    if counts['cov']:
        if counts['sub_ct'] / float(counts['cov']) >= sub_thresh:
            pos['sub'].append(ind)
            over_thresh = True
        if counts['del_ct'] / float(counts['cov']) >= del_thresh:
            pos['del'].append(ind)
            over_thresh = True
        if counts['ins_ct'] / float(counts['cov']) >= ins_thresh:
            pos['ins'].append(ind)
            over_thresh = True

        if over_thresh:
            pos['total'].append(ind)

    return pos


def find_divergence(alignment_path, contigs_path, contigs_info,
                    frequency_path, positions_path, div_sum_path,
                    min_aln_rate, platform, num_proc,
                    sub_thresh, del_thresh, ins_thresh):
    """
    Main function: takes in an alignment and finds the divergent positions
    """
    if not os.path.isfile(alignment_path) or not os.path.isfile(contigs_path):
        ctg_profile = []
        positions = _write_frequency_path(frequency_path, ctg_profile,
                                          sub_thresh, del_thresh, ins_thresh)
        total_header = "".join(["Total_positions_{0}_".format(len(positions["total"])),
                              "with_thresholds_sub_{0}".format(sub_thresh),
                              "_del_{0}_ins_{1}".format(del_thresh, ins_thresh)])
        sub_header = "".join(["Sub_positions_{0}_".format(len(positions["sub"])),
                              "with_threshold_sub_{0}".format(sub_thresh)])
        del_header = "".join(["Del_positions_{0}_".format(len(positions["del"])),
                              "with_threshold_del_{0}".format(del_thresh)])
        ins_header = "".join(["Ins_positions_{0}_".format(len(positions["ins"])),
                              "with_threshold_ins_{0}".format(ins_thresh)])
        _write_positions(positions_path, positions, total_header,
                         sub_header, del_header, ins_header)

        window_len = 1000
        sum_header = "Tentative Divergent Position Summary"
        _write_div_summary(div_sum_path, sum_header, positions,
                          len(ctg_profile), window_len)
        return

    contigs_fasta = fp.read_sequence_dict(contigs_path)
    manager = multiprocessing.Manager()
    aln_reader = SynchronizedSamReader(alignment_path, contigs_fasta, manager,
                                       config.vals["max_read_coverage"])
    chunk_feeder = SynchonizedChunkManager(contigs_fasta, manager)

    results_queue = manager.Queue()
    error_queue = manager.Queue()

    process_in_parallel(_thread_worker, (aln_reader, chunk_feeder, contigs_info, platform,
                                         results_queue, error_queue), num_proc)

    if not error_queue.empty():
        raise error_queue.get()

    total_aln_errors = []
    while not results_queue.empty():
        _, ctg_profile, aln_errors = results_queue.get()

        positions = _write_frequency_path(frequency_path, ctg_profile,
                                          sub_thresh, del_thresh, ins_thresh)
        total_header = "".join(["Total_positions_{0}_".format(len(positions["total"])),
                              "with_thresholds_sub_{0}".format(sub_thresh),
                              "_del_{0}_ins_{1}".format(del_thresh, ins_thresh)])
        sub_header = "".join(["Sub_positions_{0}_".format(len(positions["sub"])),
                              "with_threshold_sub_{0}".format(sub_thresh)])
        del_header = "".join(["Del_positions_{0}_".format(len(positions["del"])),
                              "with_threshold_del_{0}".format(del_thresh)])
        ins_header = "".join(["Ins_positions_{0}_".format(len(positions["ins"])),
                              "with_threshold_ins_{0}".format(ins_thresh)])
        _write_positions(positions_path, positions, total_header,
                         sub_header, del_header, ins_header)

        window_len = 1000
        sum_header = "Tentative Divergent Position Summary"
        _write_div_summary(div_sum_path, sum_header, positions,
                          len(ctg_profile), window_len)

        logger.debug("Total positions: %d", len(positions["total"]))
        total_aln_errors.extend(aln_errors)

    mean_aln_error = sum(total_aln_errors) / (len(total_aln_errors) + 1)
    logger.debug("Alignment error rate: %f", mean_aln_error)


def _write_frequency_path(frequency_path, ctg_profile, sub_thresh,
                          del_thresh, ins_thresh):
    #The set of called positions for each category
    positions = {"total":[], "sub":[], "del":[], "ins":[]}
    with open(frequency_path, 'w') as f:
        f.write("Index\tCov\tMatch\tCount\tSub\tCount\tDel\tCount\tIns\tCount\n")
        for index, elem in enumerate(ctg_profile):
            counts = _count_freqs(elem)
            f.write("{0}\t{c[cov]}\t{c[mat_base]}\t{c[mat_ct]}\t".format(index,c=counts))
            f.write("{c[sub_base]}\t{c[sub_ct]}\t{c[del_base]}\t".format(c=counts))
            f.write("{c[del_ct]}\t{c[ins_base]}\t{c[ins_ct]}\n".format(c=counts))

            #Adds this element to positions if it passes any threshold
            #and updates total_called appropriately
            positions = _call_position(index, counts, positions,
                                       sub_thresh, del_thresh, ins_thresh)
    return positions


def read_frequency_path(frequency_path):
    header = []
    freqs = []
    int_inds = [0,1,3,5,7,9]
    with open(frequency_path, "r") as f:
        for i,line in enumerate(f):
            line = line.strip()
            if i == 0:
                header = line.split("\t")
            else:
                vals = line.split("\t")
                for j in int_inds:
                    vals[j] = int(vals[j])
                freqs.append(vals)
    return header, freqs


def _write_positions(positions_path, positions, total_header, sub_header,
                     del_header, ins_header):
    with open(positions_path, 'w') as f:
        f.write(">{0}\n".format(total_header))
        f.write(",".join([str(x) for x in sorted(positions["total"])]) + "\n")
        f.write(">{0}\n".format(sub_header))
        f.write(",".join([str(x) for x in sorted(positions["sub"])]) + "\n")
        f.write(">{0}\n".format(del_header))
        f.write(",".join([str(x) for x in sorted(positions["del"])]) + "\n")
        f.write(">{0}\n".format(ins_header))
        f.write(",".join([str(x) for x in sorted(positions["ins"])]) + "\n")


def _write_div_summary(div_sum_path, sum_header, positions,
                      seq_len, window_len):
    pos_list = sorted(positions["total"])
    av_div = 0.0
    if seq_len != 0:
        av_div = len(pos_list) / float(seq_len)

    position_gaps = [0 for _ in range(len(pos_list)+1)]
    curr_pos = 0
    for i,p in enumerate(pos_list):
        position_gaps[i] = p-curr_pos
        curr_pos = p
    position_gaps[-1] = seq_len-curr_pos

    mean_position_gap = _mean(position_gaps)
    max_position_gap = max(position_gaps)

    window_len = 1000
    position_counts = [0 for _ in range(((seq_len - 1) // window_len) + 1)]
    window_divs = [0.0 for _ in range(((seq_len - 1) // window_len) + 1)]
    curr_p_i = 0
    for i in range(len(window_divs)):
        start = i*window_len
        end = (i+1)*window_len-1
        if i == len(window_divs)-1:
            end = seq_len-1

        curr_window_len = end-start+1

        if curr_p_i < len(pos_list) and pos_list[curr_p_i] < start:
            raise PositionIOError('Problem with position indices')
        while curr_p_i < len(pos_list) and pos_list[curr_p_i] <= end:
            position_counts[i] += 1
            curr_p_i += 1

        window_divs[i] = 0.0
        if curr_window_len != 0:
            window_divs[i] = position_counts[i] / float(curr_window_len)

    mean_window_div = _mean(window_divs)
    median_window_div = _get_median(window_divs)
    min_window_div = min(window_divs)

    with open(div_sum_path, 'w') as f:
        f.write("{0}\n\n".format(sum_header))

        f.write("{0:33}\t{1}\n".format("Sequence Length:", seq_len))
        f.write("{0:33}\t{1:.4f}\n\n".format("Average Divergence:", av_div))

        f.write("{0:33}\t{1}\n".format("Total Substitution Positions:",
                                            len(positions["sub"])))
        f.write("{0:33}\t{1}\n".format("Total Deletion Positions:",
                                            len(positions["del"])))
        f.write("{0:33}\t{1}\n".format("Total Insertion Positions:",
                                            len(positions["ins"])))
        f.write("{0:33}\t{1}\n".format("Total Positions:",
                                              len(positions["total"])))
        mixed_count = (len(positions["sub"]) + len(positions["del"]) +
                    len(positions["ins"])) - len(positions["total"])
        f.write("{0:33}\t{1}\n\n".format("Mixed Positions:", mixed_count))

        f.write("{0:33}\t{1:.2f}\n".format("Mean Position Gap:",
                                           mean_position_gap))
        f.write("{0:33}\t{1}\n".format("Max Position Gap:", max_position_gap))

        f.write("{0:33}\t{1}\n".format("Window Length:", window_len))
        f.write("{0:33}\t{1:.5f}\n".format("Mean Window Divergence:",
                                           mean_window_div))
        f.write("{0:33}\t{1:.5f}\n".format("Median Window Divergence:",
                                           median_window_div))
        f.write("{0:33}\t{1:.5f}\n".format("Min Window Divergence:",
                                           min_window_div))


class PositionIOError(Exception):
    pass


def read_positions(positions_file):
    """
    Reads positions file into list
    """
    headers = {"total":"", "sub":"", "ins":"", "del":""}
    positions = {"total":[], "sub":[], "ins":[], "del":[]}
    try:
        with open(positions_file, "r") as f:
            for line_id, line in enumerate(f):
                line = line.strip()
                if line_id == 0 and line.startswith(">") and line:
                    headers["total"] = line[1:]
                elif line_id == 1 and line:
                    pos_parts = line.split(",")
                    positions["total"] = [int(x) for x in pos_parts]
                elif line_id == 2 and line.startswith(">") and line:
                    headers["sub"] = line[1:]
                elif line_id == 3 and line:
                    pos_parts = line.split(",")
                    positions["sub"] = [int(x) for x in pos_parts]
                elif line_id == 4 and line.startswith(">") and line:
                    headers["del"] = line[1:]
                elif line_id == 5 and line:
                    pos_parts = line.split(",")
                    positions["del"] = [int(x) for x in pos_parts]
                elif line_id == 6 and line.startswith(">") and line:
                    headers["ins"] = line[1:]
                elif line_id == 7 and line:
                    pos_parts = line.split(",")
                    positions["ins"] = [int(x) for x in pos_parts]
                elif line:
                    raise PositionIOError("Not a valid positions file")
    except IOError as e:
        raise PositionIOError(e)
    return headers, positions


def _get_median(lst):
    if not lst:
        raise ValueError("_get_median() arg is an empty sequence")
    sorted_list = sorted(lst)
    if len(lst) % 2 == 1:
        return sorted_list[len(lst) // 2]
    else:
        mid1 = sorted_list[(len(lst) // 2) - 1]
        mid2 = sorted_list[(len(lst) // 2)]
        return (mid1 + mid2) / 2


def _mean(lst):
    if not lst:
        return 0
    return sum(lst) / len(lst)
