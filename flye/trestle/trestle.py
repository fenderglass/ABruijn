#(c) 2016-2018 by Authors
#This file is a part of Flye program.
#Released under the BSD license (see LICENSE file)

"""
Created on Wed Jan  4 03:50:31 2017

@author: jeffrey_yuan
"""

import os
import logging
from itertools import combinations, product
import copy
import multiprocessing, signal
import traceback

import flye.polishing.alignment as flye_aln
from flye.utils.sam_parser import SynchronizedSamReader, Alignment
import flye.utils.fasta_parser as fp
import flye.config.py_cfg as config
import flye.polishing.polish as pol

import flye.trestle.divergence as div
import flye.trestle.trestle_config as trestle_config

logger = logging.getLogger()


def resolve_repeats(args, trestle_dir, repeats_info, summ_file,
                    resolved_repeats_seqs):
    all_file_names = define_file_names()
    all_labels, initial_file_names = all_file_names[0], all_file_names[2]

    all_resolved_reps_dict = {}
    all_summaries = []
    init_summary(summ_file)

    #1. Process repeats from graph - generates a folder for each repeat
    logger.debug("Finding unbridged repeats")
    process_outputs = process_repeats(args.reads, repeats_info,
                                      trestle_dir, all_labels,
                                      initial_file_names)
    repeat_list, repeat_edges, all_edge_headers = process_outputs
    logger.info("Simple unbridged repeats: %d", len(repeat_list))

    #if not repeat_list:
    #    return

    #Resolve every repeat in a separate thread
    def _thread_worker(func_args, log_file, results_queue, error_queue):
        try:
            #each thred logs to a separate file
            log_formatter = \
                logging.Formatter("[%(asctime)s] %(name)s: %(levelname)s: "
                                  "%(message)s", "%Y-%m-%d %H:%M:%S")
            file_handler = logging.FileHandler(log_file, mode="a")
            file_handler.setFormatter(log_formatter)
            for handler in logger.handlers[:]:
                logger.removeHandler(handler)
            logger.addHandler(file_handler)

            result = resolve_each_repeat(*func_args)
            results_queue.put(result)

        except Exception as e:
            logger.error("Thread exception")
            logger.error(traceback.format_exc())
            error_queue.put(e)

    job_chunks = [repeat_list[i:i + args.threads]
                  for i in range(0, len(repeat_list), args.threads)]

    for job_chunk in job_chunks:
        manager = multiprocessing.Manager()
        results_queue = manager.Queue()
        error_queue = manager.Queue()

        repeat_threads = max(1, args.threads // len(job_chunk))
        orig_sigint = signal.signal(signal.SIGINT, signal.SIG_IGN)
        threads = []
        for rep_id in sorted(job_chunk):
            func_args = (rep_id, repeat_edges, all_edge_headers, args, trestle_dir,
                         repeats_info, all_file_names, repeat_threads)
            log_file = os.path.join(trestle_dir,
                                    "repeat_{0}".format(rep_id), "log.txt")
            threads.append(multiprocessing.Process(target=_thread_worker,
                                        args=(func_args, log_file,
                                              results_queue, error_queue)))
        signal.signal(signal.SIGINT, orig_sigint)

        for t in threads:
            t.start()
        try:
            for t in threads:
                t.join()
                if t.exitcode == -9:
                    logger.error("Looks like the system ran out of memory")
                if t.exitcode != 0:
                    raise Exception("One of the processes exited with code: {0}"
                                    .format(t.exitcode))

        except KeyboardInterrupt:
            for t in threads:
                t.terminate()
            raise

        while not error_queue.empty():
            logger.warning("Non-critical error in trestle thread: " + str(error_queue.get()))
        #if not error_queue.empty():
        #    raise error_queue.get()

        while not results_queue.empty():
            resolved_dict, summary_list = results_queue.get()
            all_resolved_reps_dict.update(resolved_dict)
            all_summaries.extend(summary_list)

    fp.write_fasta_dict(all_resolved_reps_dict, resolved_repeats_seqs)
    num_resolved = 0
    for summ_items in all_summaries:
        if summ_items[6]:
            num_resolved += 1
        update_summary(summ_items, summ_file)
    logger.info("Resolved: %d", num_resolved)


def resolve_each_repeat(rep_id, repeat_edges, all_edge_headers, args,
                        trestle_dir, repeats_info, all_file_names,
                        num_threads):
    SUB_THRESH = trestle_config.vals["sub_thresh"]
    DEL_THRESH = trestle_config.vals["del_thresh"]
    INS_THRESH = trestle_config.vals["ins_thresh"]
    MAX_ITER = trestle_config.vals["max_iter"]
    MIN_ALN_RATE = trestle_config.vals["min_aln_rate"]
    NUM_POL_ITERS = trestle_config.vals["num_pol_iters"]
    ORIENT_CONFIG = trestle_config.vals["orientations_to_run"]
    zero_it = 0

    (all_labels, pol_dir_names, initial_file_names,
     pre_file_names, div_file_names, aln_names,
     middle_file_names, output_file_names) = all_file_names

    repeat_label, side_labels = all_labels
    pol_temp_name, pol_ext_name, pol_cons_name = pol_dir_names
    (template_name, extended_name, repeat_reads_name,
     pre_partitioning_name) = initial_file_names
    pre_edge_reads_name, pre_read_aln_name, partitioning_name = pre_file_names
    div_freq_name, div_pos_name, div_summ_name = div_file_names
    (reads_template_aln_name, cons_temp_aln_name,
     cut_cons_temp_aln_name, reads_cons_aln_name) = aln_names
    (confirmed_pos_name, edge_reads_name,
     cut_cons_name, cons_vs_cons_name) = middle_file_names
    (side_stats_name, int_stats_name, int_confirmed_pos_name,
     resolved_rep_name, res_vs_res_name) = output_file_names

    logger.info("Resolving repeat %d: %s",
                rep_id, repeats_info[rep_id].repeat_path)
    repeat_dir = os.path.join(trestle_dir,
                              repeat_label.format(rep_id))

    run_orientations = []
    if ORIENT_CONFIG == "forward":
        run_orientations = [("forward", rep_id)]
    elif ORIENT_CONFIG == "reverse":
        run_orientations = [("reverse", -rep_id)]
    elif ORIENT_CONFIG == "both":
        run_orientations = [("forward", rep_id), ("reverse", -rep_id)]
    repeat_bridged = False
    resolved_dict = {}
    summary_list = []
    for orientation, rep in run_orientations:
        logger.debug("Orientation: " + orientation)
        orient_dir = os.path.join(repeat_dir, orientation)
        template = os.path.join(orient_dir, template_name)
        extended = os.path.join(orient_dir, extended_name)
        repeat_reads = os.path.join(orient_dir, repeat_reads_name)
        term_bool = {s:False for s in side_labels}

        #2. Polish template and extended templates
        logger.debug("Polishing templates")
        pol_temp_dir = os.path.join(orient_dir, pol_temp_name)
        if not os.path.isdir(pol_temp_dir):
            os.mkdir(pol_temp_dir)
        polished_template, _ = \
            pol.polish(template, [repeat_reads], pol_temp_dir, NUM_POL_ITERS,
                       num_threads, args.platform, args.read_type, output_progress=False)

        if not os.path.getsize(polished_template):
            for side in side_labels:
                term_bool[side] = True

        polished_extended = {}
        pol_ext_dir = os.path.join(orient_dir, pol_ext_name)
        for side in side_labels:
            for edge_id in repeat_edges[rep][side]:
                if not os.path.isdir(pol_ext_dir.format(side, edge_id)):
                    os.mkdir(pol_ext_dir.format(side, edge_id))
                pol_output, _ = \
                    pol.polish(extended.format(side, edge_id), [repeat_reads],
                               pol_ext_dir.format(side, edge_id), NUM_POL_ITERS,
                               num_threads, args.platform, args.read_type,
                               output_progress=False)
                polished_extended[(side, edge_id)] = pol_output
                if not os.path.getsize(pol_output):
                    term_bool[side] = True

        #3. Find divergent positions
        logger.debug("Estimating divergence")
        frequency_path = os.path.join(orient_dir, div_freq_name)
        position_path = os.path.join(orient_dir, div_pos_name)
        summary_path = os.path.join(orient_dir, div_summ_name)

        #logger.info("running Minimap2")
        alignment_file = os.path.join(orient_dir, reads_template_aln_name)
        template_len = 0.0
        if os.path.getsize(polished_template):
            flye_aln.make_alignment(polished_template, [repeat_reads],
                       num_threads, orient_dir, args.platform,
                       alignment_file, reference_mode=True, sam_output=True)
            template_info = flye_aln.get_contigs_info(polished_template)
            template_len = template_info[str(rep)].length

        logger.debug("Finding tentative divergence positions")
        div.find_divergence(alignment_file, polished_template,
                            template_info, frequency_path, position_path,
                            summary_path, MIN_ALN_RATE,
                            args.platform, num_threads,
                            SUB_THRESH, DEL_THRESH, INS_THRESH)
        read_endpoints = find_read_endpoints(alignment_file,
                                             polished_template)
        avg_cov = find_coverage(frequency_path)

        #4. Initialize paths, variables, and stats
        pre_partitioning = os.path.join(orient_dir, pre_partitioning_name)
        pre_edge_reads = os.path.join(orient_dir, pre_edge_reads_name)
        pre_read_align = os.path.join(orient_dir, pre_read_aln_name)
        partitioning = os.path.join(orient_dir, partitioning_name)
        cons_align = os.path.join(orient_dir, cons_temp_aln_name)
        cut_cons_align = os.path.join(orient_dir, cut_cons_temp_aln_name)
        read_align = os.path.join(orient_dir, reads_cons_aln_name)
        confirmed_pos_path = os.path.join(orient_dir, confirmed_pos_name)
        edge_reads = os.path.join(orient_dir, edge_reads_name)
        cut_cons = os.path.join(orient_dir, cut_cons_name)
        polishing_dir = os.path.join(orient_dir, pol_cons_name)
        cons_vs_cons = os.path.join(orient_dir, cons_vs_cons_name)
        side_stats = os.path.join(orient_dir, side_stats_name)
        integrated_stats = os.path.join(orient_dir, int_stats_name)
        int_confirmed_path = os.path.join(orient_dir,
                                          int_confirmed_pos_name)
        resolved_rep_path = os.path.join(orient_dir, resolved_rep_name)
        res_vs_res = os.path.join(orient_dir, res_vs_res_name)

        #5. Re-align reads to extended and initialize partitioning 0
        logger.debug("Checking initial set of edge reads")
        for side in side_labels:
            for edge_id in repeat_edges[rep][side]:
                write_edge_reads(zero_it, side, edge_id,
                                 repeat_reads,
                                 pre_partitioning.format(side),
                                 pre_edge_reads.format(side, edge_id))
                flye_aln.make_alignment(polished_extended[(side, edge_id)],
                                        [pre_edge_reads.format(side, edge_id)],
                                        num_threads, orient_dir, args.platform,
                                        pre_read_align.format(side, edge_id),
                                        reference_mode=True, sam_output=True)
            init_partitioning(repeat_edges[rep][side],
                              side, pre_partitioning.format(side),
                              pre_read_align, polished_extended,
                              partitioning.format(zero_it, side))

        cut_consensus = {}
        side_it = {s:0 for s in side_labels}
        iter_pairs = []
        edge_below_cov = {s:False for s in side_labels}
        dup_part = {s:False for s in side_labels}
        prev_partitionings = {s:set() for s in side_labels}
        #6. Initialize stats
        for side in side_labels:
            edge_below_cov[side] = init_side_stats(
                                rep, side, repeat_edges, args.min_overlap,
                                position_path,
                                partitioning.format(zero_it, side),
                                prev_partitionings[side],
                                template_len,
                                side_stats.format(side))
        init_int_stats(rep, repeat_edges, zero_it, position_path,
                       partitioning, repeat_reads, template_len,
                       avg_cov, integrated_stats)
        #7. Start iterations
        logger.debug("Iterative procedure")
        for it in range(1, MAX_ITER + 1):
            both_break = True
            for side in side_labels:
                if (edge_below_cov[side] or dup_part[side] or
                    term_bool[side]):
                    continue
                else:
                    logger.debug("Iteration %d, '%s'", it, side)
                    both_break = False
                for edge_id in sorted(repeat_edges[rep][side]):
                    #7a. Call consensus on partitioned reads
                    pol_con_dir = polishing_dir.format(
                                it, side, edge_id)
                    curr_reads = edge_reads.format(it, side, edge_id)
                    write_edge_reads(
                                it, side, edge_id,
                                repeat_reads,
                                partitioning.format(it - 1, side),
                                curr_reads)
                    curr_extended = polished_extended[(side, edge_id)]
                    logger.debug("\tPolishing '%s %s' reads", side, edge_id)

                    if not os.path.isdir(pol_con_dir):
                        os.mkdir(pol_con_dir)
                    pol_con_out, _ = \
                        pol.polish(curr_extended, [curr_reads], pol_con_dir,
                                   NUM_POL_ITERS, num_threads, args.platform, args.read_type,
                                   output_progress=False)
                    #7b. Cut consensus where coverage drops
                    cutpoint = locate_consensus_cutpoint(
                                    side, read_endpoints,
                                    curr_reads)
                    if os.path.getsize(pol_con_out):
                        cons_al_file = cons_align.format(it, side, edge_id)
                        flye_aln.make_alignment(polished_template, [pol_con_out],
                                                num_threads, orient_dir,
                                                args.platform, cons_al_file,
                                                reference_mode=True, sam_output=True)
                    else:
                        term_bool[side] = True
                    curr_cut_cons = cut_cons.format(it, side, edge_id)
                    cut_consensus[(it, side, edge_id)] = curr_cut_cons
                    if os.path.isfile(cons_al_file):
                        truncate_consensus(side, cutpoint, cons_al_file,
                                           polished_template,
                                           pol_con_out, curr_cut_cons)
                    else:
                        term_bool[side] = True
                    #7c. Align consensuses to template 
                    #    and reads to consensuses
                    if os.path.isfile(curr_cut_cons):
                        cut_cons_al_file = cut_cons_align.format(it, side, edge_id)
                        flye_aln.make_alignment(polished_template, [curr_cut_cons],
                                                num_threads, orient_dir,
                                                args.platform, cut_cons_al_file,
                                                reference_mode=True, sam_output=True)
                        read_al_file = read_align.format(it, side, edge_id)
                        flye_aln.make_alignment(curr_cut_cons, [repeat_reads],
                                                num_threads, orient_dir,
                                                args.platform, read_al_file,
                                                reference_mode=True, sam_output=True)
                    else:
                        term_bool[side] = True
                #7d. Partition reads using divergent positions
                logger.debug("\tPartitioning '%s' reads", side)
                partition_reads(repeat_edges[rep][side], it, side,
                                   position_path, cut_cons_align,
                                   polished_template, read_align,
                                   cut_consensus, confirmed_pos_path,
                                   partitioning, all_edge_headers[rep])
                #7e. Write stats file for current iteration
                edge_pairs = sorted(combinations(repeat_edges[rep][side],
                                                 2))
                for edge_one, edge_two in edge_pairs:
                    cons_one = cut_consensus[(it, side, edge_one)]
                    cons_two = cut_consensus[(it, side, edge_two)]
                    if (not os.path.isfile(cons_one) or
                        not os.path.isfile(cons_two)):
                        continue
                    cons_cons_file = cons_vs_cons.format(
                                            it, side, edge_one,
                                            it, side, edge_two)
                    flye_aln.make_alignment(cons_two, [cons_one],
                                            num_threads, orient_dir,
                                            args.platform, cons_cons_file,
                                            reference_mode=True, sam_output=True)
                side_stat_outputs = update_side_stats(
                                    repeat_edges[rep][side], it, side,
                                    cut_cons_align, polished_template,
                                    confirmed_pos_path.format(it, side),
                                    partitioning.format(it, side),
                                    prev_partitionings[side],
                                    side_stats.format(side))
                edge_below_cov[side], dup_part[side] = side_stat_outputs
                side_it[side] = it
            iter_pairs.append((side_it[side_labels[0]],
                               side_it[side_labels[1]]))
            update_int_stats(rep, repeat_edges, side_it, cut_cons_align,
                                polished_template,
                                template_len,
                                confirmed_pos_path, int_confirmed_path,
                                partitioning, integrated_stats)
            if both_break:
                break
        #8. Finalize stats files
        logger.debug("Writing stats files")
        for side in side_labels:
            finalize_side_stats(repeat_edges[rep][side], side_it[side],
                                side, cut_cons_align, polished_template,
                                cons_vs_cons, cut_consensus,
                                confirmed_pos_path.format(side_it[side],
                                                          side),
                                partitioning.format(side_it[side], side),
                                edge_below_cov[side],
                                dup_part[side], term_bool[side],
                                side_stats.format(side))
        final_int_outputs = finalize_int_stats(rep, repeat_edges, side_it,
                                               cut_cons_align,
                                               polished_template,
                                               template_len, cons_vs_cons,
                                               cut_consensus,
                                               int_confirmed_path,
                                               partitioning,
                                               integrated_stats,
                                               resolved_rep_path)
        bridged, repeat_seqs, summ_vals = final_int_outputs
        #9. Generate summary and resolved repeat file
        logger.debug("Generating summary and resolved repeat file")
        avg_div = 0.0
        both_resolved_present = False
        if bridged:
            res_inds = list(range(len(repeat_edges[rep]["in"])))
            for res_one, res_two in sorted(combinations(res_inds, 2)):
                res_one_path = resolved_rep_path.format(rep, res_one)
                res_two_path = resolved_rep_path.format(rep, res_two)
                if (os.path.isfile(res_one_path) and
                    os.path.isfile(res_two_path)):
                    both_resolved_present = True
                    repeat_bridged = True
                    flye_aln.make_alignment(res_two_path, [res_one_path],
                                        num_threads, orient_dir,
                                        args.platform,
                                        res_vs_res.format(rep, res_one, res_two),
                                        reference_mode=True, sam_output=True)
            if both_resolved_present:
                avg_div = int_stats_postscript(rep, repeat_edges,
                                               integrated_stats,
                                               resolved_rep_path,
                                               res_vs_res)
        if both_resolved_present:
            resolved_dict.update(repeat_seqs)
        summary_list.append((rep, repeats_info[rep].repeat_path, template_len,
                             avg_cov, summ_vals, avg_div,
                             both_resolved_present))
        remove_unneeded_files(repeat_edges, rep, side_labels, side_it,
                              orient_dir, template, extended, pol_temp_dir,
                              pol_ext_dir, pre_edge_reads,
                              pre_partitioning, pre_read_align,
                              partitioning, cons_align, cut_cons_align,
                              read_align, confirmed_pos_path, edge_reads,
                              cut_cons, polishing_dir, cons_vs_cons,
                              int_confirmed_path, repeat_reads,
                              frequency_path, alignment_file,
                              NUM_POL_ITERS, iter_pairs)
    if repeat_bridged:
        logger.info("Repeat successfully resolved")
    else:
        logger.info("Repeat not resolved")
    return resolved_dict, summary_list

def define_file_names():
    #Defining directory and file names for trestle output
    repeat_label = "repeat_{0}"
    side_labels = ["in", "out"]
    all_labels = repeat_label, side_labels

    pol_temp_name = "Polishing.Template"
    pol_ext_name = "Polishing.Extended.{0}.{1}"
    pol_cons_name = "Polishing.Consensus.{0}.{1}.{2}"
    pol_dir_names = pol_temp_name, pol_ext_name, pol_cons_name

    template_name = "template.fasta"
    extended_name = "extended_templates.{0}.{1}.fasta"
    repeat_reads_name = "repeat_reads.fasta"
    pre_partitioning_name = "pre_partitioning.{0}.txt"
    initial_file_names = (template_name, extended_name, repeat_reads_name,
                          pre_partitioning_name)


    pre_edge_reads_name = "pre_edge_reads.{0}.{1}.txt"
    pre_read_aln_name = "pre_edge_reads.{0}.{1}.vs.extended.minimap.bam"
    partitioning_name = "partitioning.{0}.{1}.txt"
    pre_file_names = pre_edge_reads_name, pre_read_aln_name, partitioning_name

    div_freq_name = "divergence_frequencies.txt"
    div_pos_name = "divergent_positions.txt"
    div_summ_name = "divergence_summary.txt"
    div_file_names = div_freq_name, div_pos_name, div_summ_name

    reads_template_aln_name = "reads.vs.template.minimap.bam"
    cons_temp_aln_name = "uncut_consensus.{0}.{1}.{2}.vs.template.minimap.bam"
    cut_cons_temp_aln_name = "consensus.{0}.{1}.{2}.vs.template.minimap.bam"
    reads_cons_aln_name = "reads.vs.consensus.{0}.{1}.{2}.minimap.bam"
    aln_names = (reads_template_aln_name, cons_temp_aln_name,
                 cut_cons_temp_aln_name, reads_cons_aln_name)

    confirmed_pos_name = "confirmed_positions.{0}.{1}.txt"
    edge_reads_name = "edge_reads.{0}.{1}.{2}.fasta"
    cut_cons_name = "consensus.{0}.{1}.{2}.fasta"
    cons_vs_cons_name = "".join(["consensus.{0}.{1}.{2}.vs.",
                                 "consensus.{3}.{4}.{5}.minimap.bam"])
    middle_file_names = (confirmed_pos_name, edge_reads_name,
                         cut_cons_name, cons_vs_cons_name)

    side_stats_name = "stats_from_{0}.txt"
    int_stats_name = "stats_integrated.txt"
    int_confirmed_pos_name = "integrated_confirmed_positions.{0}.{1}.txt"
    resolved_rep_name = "resolved_repeat_{0}.copy.{1}.fasta"
    res_vs_res_name = "resolved_repeat_{0}.copy.{1}.vs.{2}.minimap.bam"
    output_file_names = (side_stats_name, int_stats_name,
                         int_confirmed_pos_name, resolved_rep_name,
                         res_vs_res_name)

    all_file_names = (all_labels, pol_dir_names, initial_file_names,
                      pre_file_names, div_file_names, aln_names,
                      middle_file_names, output_file_names)
    return all_file_names


#Process Repeats functions
class ProcessingException(Exception):
    pass


def process_repeats(reads, repeats_dict, work_dir, all_labels,
                    initial_file_names):
    """
    Generates repeat dirs and files given reads, repeats_dump and
    graph_edges files. Only returns repeats between min_mult and max_mult
    """

    if not repeats_dict:
        return [], {}, {}

    #creates a separate process to make sure that
    #read dictionary is released after the function exits
    manager = multiprocessing.Manager()
    return_queue = manager.Queue()

    orig_sigint = signal.signal(signal.SIGINT, signal.SIG_IGN)
    thread = multiprocessing.Process(target=_process_repeats_impl,
                                        args=(reads, repeats_dict,
                                              work_dir, all_labels,
                                              initial_file_names,
                                              return_queue))
    signal.signal(signal.SIGINT, orig_sigint)
    thread.start()
    try:
        thread.join()
        if thread.exitcode == -9:
            logger.error("Looks like the system ran out of memory")
        if thread.exitcode != 0:
            raise Exception("One of the processes exited with code: {0}"
                                .format(thread.exitcode))
    except KeyboardInterrupt:
        thread.terminate()
        raise

    return return_queue.get()


def _process_repeats_impl(reads, repeats_dict, work_dir, all_labels,
                          initial_file_names, return_queue):
    """
    This function is called in a separate process
    """
    MIN_MULT = trestle_config.vals["min_mult"]
    MAX_MULT = trestle_config.vals["max_mult"]
    FLANKING_LEN = trestle_config.vals["flanking_len"]
    ORIENT_CONFIG = trestle_config.vals["orientations_to_run"]

    repeat_label, side_labels = all_labels
    (template_name, extended_name, repeat_reads_name,
     pre_partitioning_name) = initial_file_names

    reads_dict = {}
    for read_file in reads:
        reads_dict.update(fp.read_sequence_dict(read_file))
    #orig_graph = fp.read_sequence_dict(graph_edges)
    #graph_dict = {int(h.split('_')[1]):orig_graph[h] for h in orig_graph}

    if not reads_dict:
        raise ProcessingException("No reads found from {0}".format(reads))
    #if not graph_dict:
    #    raise ProcessingException("No edges found from {0}".format(
    #        graph_edges))

    repeat_list = []
    repeat_edges = {}
    all_edge_headers = {}
    for rep in sorted(repeats_dict, reverse=True):
        #Checks multiplicity of repeat and presence of reverse strand
        #One run processes both forward and reverse strand of repeat

        if rep <= 0:
            continue

        valid_repeat = True
        if -rep not in repeats_dict:
            logger.debug("Repeat %s missing reverse strand", rep)
            valid_repeat = False
        elif (repeats_dict[rep].multiplicity < MIN_MULT or
                 repeats_dict[rep].multiplicity > MAX_MULT or
                 repeats_dict[-rep].multiplicity < MIN_MULT or
                 repeats_dict[-rep].multiplicity > MAX_MULT):
            logger.debug("Repeat %s multiplicity not in range: %s",
                         rep, repeats_dict[rep].multiplicity)
            valid_repeat = False
        #if rep not in graph_dict:
        #    logger.debug("Repeat {0} missing from graph file".format(rep))
        #    valid_repeat = False
        if not valid_repeat:
            continue

        #Makes repeat dirs
        repeat_dir = os.path.join(work_dir, repeat_label.format(rep))
        if not os.path.isdir(repeat_dir):
            os.mkdir(repeat_dir)
        repeat_list.append(rep)

        run_orientations = []
        if ORIENT_CONFIG == "forward":
            run_orientations = [("forward", rep)]
        elif ORIENT_CONFIG == "reverse":
            run_orientations = [("reverse", -rep)]
        elif ORIENT_CONFIG == "both":
            run_orientations = [("forward", rep), ("reverse", -rep)]
        for curr_label, curr_rep in run_orientations:
            orient_path = os.path.join(repeat_dir, curr_label)
            if not os.path.isdir(orient_path):
                os.mkdir(orient_path)
            template_path = os.path.join(orient_path, template_name)
            extended_path = os.path.join(orient_path, extended_name)
            repeat_reads_path = os.path.join(orient_path, repeat_reads_name)
            partitioning_path = os.path.join(orient_path,
                                             pre_partitioning_name)

            in_label = side_labels[0]
            out_label = side_labels[1]
            repeat_edges[curr_rep] = {in_label:[], out_label:[]}

            #(mult, all_reads_list, inputs_dict,
            # outputs_dict) = repeats_dict[curr_rep]
            #mult = repeats_dict[curr_rep].multiplicity
            all_reads_list = repeats_dict[curr_rep].all_reads
            inputs_dict = repeats_dict[curr_rep].in_reads
            outputs_dict = repeats_dict[curr_rep].out_reads

            template_dict = {}
            extended_dicts = {}
            repeat_reads_dict = {}
            #Partitioning parts: id_num, Partitioned/Tied/None, 
            #edge_id, top_score, total_score, Header
            partitioning = {in_label:[], out_label:[]}
            read_id = 0

            template_seq = repeats_dict[curr_rep].sequences["template"]
            #if curr_label == "reverse":
            #    template_seq = fp.reverse_complement(graph_dict[rep])
            template_dict[curr_rep] = template_seq

            all_edge_headers[curr_rep] = {}
            out_headers = set()
            #Headers will be in the form -h or +h,
            #edge_dict is in the form >[Input,Output]_edge##_h,
            #rev_comp of read will be written if the header is -h
            for edge_id in inputs_dict:
                repeat_edges[curr_rep][in_label].append(edge_id)
                extended_dicts[(in_label, edge_id)] = {}

                headers = inputs_dict[edge_id]
                for header in headers:
                    if (not header) or (header[0] != '+' and header[0] != '-'):
                        raise ProcessingException(
                            "Input read format not recognized: {0}".format(
                                header))
                    if header[1:] not in reads_dict:
                        raise ProcessingException(
                            "Read header {0} not in any of {1}".format(
                                header[1:], reads))

                    if header[1:] not in all_edge_headers[curr_rep]:
                        status_label = "Partitioned"
                        edge_label = str(edge_id)
                        score = 1
                        total_score = 0
                        partitioning[in_label].append((read_id, status_label,
                                                       edge_label, score,
                                                       total_score,
                                                       header[1:]))
                        all_edge_headers[curr_rep][header[1:]] = read_id
                        read_id += 1

                extend_in_header = "Extended_Template_Input_{0}".format(
                    edge_id)
                #if edge_id > 0:
                #    edge_seq = graph_dict[edge_id]
                #elif edge_id < 0:
                #    edge_seq = fp.reverse_complement(graph_dict[-edge_id])
                edge_seq = repeats_dict[curr_rep].sequences[edge_id]
                extended_seq = edge_seq[-FLANKING_LEN:]
                extended_dicts[(in_label, edge_id)][extend_in_header] = (
                                        extended_seq + template_seq)

            for edge_id in outputs_dict:
                repeat_edges[curr_rep][out_label].append(edge_id)
                extended_dicts[(out_label, edge_id)] = {}

                headers = outputs_dict[edge_id]
                for header in headers:
                    if (not header) or (header[0] != '+' and header[0] != '-'):
                        raise ProcessingException(
                        "Output read format not recognized: {0}".format(
                            header))
                    if header[1:] not in reads_dict:
                        raise ProcessingException(
                        "Read header {0} not in any of {1}".format(
                            header[1:], reads))

                    curr_read_id = read_id
                    if header[1:] not in all_edge_headers[curr_rep]:
                        status_label = "None"
                        edge_label = "NA"
                        score = 0
                        total_score = 0
                        partitioning[in_label].append((read_id, status_label,
                                                       edge_label, score,
                                                       total_score,
                                                       header[1:]))

                        all_edge_headers[curr_rep][header[1:]] = read_id
                        read_id += 1
                    else:
                        curr_read_id = all_edge_headers[curr_rep][header[1:]]

                    if header[1:] not in out_headers:
                        status_label = "Partitioned"
                        edge_label = str(edge_id)
                        score = 1
                        total_score = 0
                        partitioning[out_label].append((curr_read_id,
                                                        status_label,
                                                        edge_label, score,
                                                        total_score,
                                                        header[1:]))
                        out_headers.add(header[1:])

                extend_out_header = "Extended_Template_Output_{0}".format(
                    edge_id)
                #if edge_id > 0:
                #    edge_seq = graph_dict[edge_id]
                #elif edge_id < 0:
                #    edge_seq = fp.reverse_complement(graph_dict[-edge_id])
                edge_seq = repeats_dict[curr_rep].sequences[edge_id]
                extended_seq = edge_seq[:FLANKING_LEN]
                extended_dicts[(out_label, edge_id)][extend_out_header] = (
                                        template_seq + extended_seq)

            #Need to reiterate over in_headers to add in_headers to 
            #out-partitioning while avoiding double-adding ones in both
            for edge_id in inputs_dict:
                headers = inputs_dict[edge_id]
                for header in headers:
                    if header[1:] not in out_headers:
                        curr_read_id = all_edge_headers[curr_rep][header[1:]]
                        status_label = "None"
                        edge_label = "NA"
                        score = 0
                        total_score = 0
                        partitioning[out_label].append((curr_read_id,
                                                        status_label,
                                                        edge_label, score,
                                                        total_score,
                                                        header[1:]))

            for header in all_reads_list:
                if (not header) or (header[0] != '+' and header[0] != '-'):
                    raise ProcessingException(
                        "All reads format not recognized: {0}".format(header))
                if header[1:] not in reads_dict:
                    raise ProcessingException(
                        "Read header {0} not in any of {1}".format(
                            header[1:], reads))

                seq = reads_dict[header[1:]]
                if header[0] == '-':
                    seq = fp.reverse_complement(seq)
                repeat_reads_dict[header[1:]] = seq

                curr_read_id = read_id
                if header[1:] not in all_edge_headers[curr_rep]:
                    all_edge_headers[curr_rep][header[1:]] = read_id
                    read_id += 1

                    status_label = "None"
                    edge_label = "NA"
                    score = 0
                    total_score = 0
                    partitioning[in_label].append((curr_read_id, status_label,
                                                   edge_label, score,
                                                   total_score, header[1:]))

                    status_label = "None"
                    edge_label = "NA"
                    score = 0
                    total_score = 0
                    partitioning[out_label].append((curr_read_id, status_label,
                                                   edge_label, score,
                                                   total_score, header[1:]))

            if template_dict and list(template_dict.values())[0]:
                fp.write_fasta_dict(template_dict, template_path)
            for edge in extended_dicts:
                if extended_dicts[edge] and list(extended_dicts[edge].values())[0]:
                    extended_edge_path = extended_path.format(edge[0],
                                                              edge[1])
                    fp.write_fasta_dict(extended_dicts[edge],
                                        extended_edge_path)
            if repeat_reads_dict and list(repeat_reads_dict.values())[0]:
                fp.write_fasta_dict(repeat_reads_dict, repeat_reads_path)
            for side in side_labels:
                _write_partitioning_file(partitioning[side],
                                         partitioning_path.format(side))

            if not template_dict:
                raise ProcessingException("No template {0} found".format(
                                                curr_rep))
            for edge in extended_dicts:
                if not template_dict:
                    raise ProcessingException(
                        "No extended template {0} {1} {2} found".format(
                            curr_rep, edge[0], edge[1]))
            if not repeat_reads_dict:
                raise ProcessingException("No repeat reads {0} found".format(
                                                curr_rep))
            for side in side_labels:
                if not partitioning[side]:
                    raise ProcessingException(
                        "Empty partitioning file {0}".format(
                            partitioning_path.format(side)))

    return_queue.put((repeat_list, repeat_edges, all_edge_headers))


def _write_partitioning_file(part_list, part_path):
    with open(part_path, "w") as f:
        header_labels = ["Read_ID", "Status", "Edge", "Top Score",
                         "Total Score", "Header"]
        spaced_header = ["{:11}".format(h) for h in header_labels]
        f.write("\t".join(spaced_header))
        f.write("\n")
        for read_label in sorted(part_list):
            spaced_label = ["{:11}".format(h) for h in read_label]
            f.write("\t".join(spaced_label))
            f.write("\n")


def _read_partitioning_file(partitioning_file):
    part_list = []
    with open(partitioning_file, "r") as f:
        for i, line in enumerate(f):
            if i > 0:
                line = line.strip()
                tokens = [t.strip() for t in line.split("\t")]
                for int_ind in [0, 3, 4]:
                    tokens[int_ind] = int(tokens[int_ind])
                part_list.append(tuple(tokens))
    return part_list


def find_coverage(frequency_file):
    coverage = 0.0
    if os.path.isfile(frequency_file):
        header, freqs = div.read_frequency_path(frequency_file)
        cov_ind = header.index("Cov")
        all_covs = [f[cov_ind] for f in freqs]
        coverage = _mean(all_covs)
        #print min(all_covs), _mean(all_covs), max(all_covs)
    return coverage


def write_edge_reads(it, side, edge_id, all_reads, partitioning, out_file):
    all_reads_dict = fp.read_sequence_dict(all_reads)
    part_list = _read_partitioning_file(partitioning)
    edge_header_name = "Read_{0}|Iter_{1}|Side_{2}|Edge_{3}|{4}"
    edge_reads = {}
    for read_id, status, edge, _, _, header in part_list:
        if status == "Partitioned" and edge != "NA" and int(edge) == edge_id:
            edge_seq = all_reads_dict[header]
            edge_header = edge_header_name.format(read_id, it,
                                                  side, edge_id, header)
            edge_reads[edge_header] = edge_seq
    if edge_reads and list(edge_reads.values())[0]:
        fp.write_fasta_dict(edge_reads, out_file)


def init_partitioning(edges, side, pre_partitioning, pre_read_align, extended,
                      partitioning):
    FLANKING_LEN = trestle_config.vals["flanking_len"]
    CONS_ALN_RATE = trestle_config.vals["cons_aln_rate"]
    #dict from read_header to edge
    extend_overlap_reads = {}
    for edge in edges:
        non_overlap_reads = 0
        aligns = _read_alignment(pre_read_align.format(side, edge),
                                 extended[(side, edge)], CONS_ALN_RATE)
        if aligns and aligns[0]:
            for aln in aligns[0]:
                edge_header = aln.qry_id
                read_header = edge_header.split("|")[-1]
                if ((side == "in" and
                        aln.trg_start < FLANKING_LEN) or
                    (side == "out" and
                        aln.trg_end >= aln.trg_len - FLANKING_LEN)):
                    extend_overlap_reads[read_header] = str(edge)
                else:
                    non_overlap_reads += 1
        logger.debug("Side %s, edge %s, non-overlap reads = %d",
                     side, edge, non_overlap_reads)
    partitioned_reads = []
    part_list = _read_partitioning_file(pre_partitioning)
    for read_id, _, edge, _, _, header in part_list:
        if header in extend_overlap_reads:
            partitioned_reads.append((read_id, "Partitioned",
                                      extend_overlap_reads[header],
                                      1, 0, header))
        else:
            partitioned_reads.append((read_id, "None", "NA", 0, 0, header))
    _write_partitioning_file(partitioned_reads, partitioning)


#Cut Consensus Functions
def find_read_endpoints(alignment_file, template):
    CONS_ALN_RATE = trestle_config.vals["cons_aln_rate"]
    read_endpoints = {}
    aligns = _read_alignment(alignment_file, template, CONS_ALN_RATE)
    if aligns and aligns[0]:
        for aln in aligns[0]:
            read_header = aln.qry_id
            start = aln.trg_start
            end = aln.trg_end
            if read_header not in read_endpoints:
                read_endpoints[read_header] = (start, end)
    else:
        logger.debug("No read alignment to template, no read_endpoints")
    return read_endpoints


def locate_consensus_cutpoint(side, read_endpoints, edge_read_file):
    MIN_EDGE_COV = trestle_config.vals["min_edge_cov"]
    all_endpoints = []
    max_endpoint = 0
    edge_reads = fp.read_sequence_dict(edge_read_file)
    for edge_header in edge_reads:
        parts = edge_header.split("|")
        read_header = parts[-1]
        if read_header in read_endpoints:
            endpoint = read_endpoints[read_header]
            if max(endpoint) > max_endpoint:
                max_endpoint = max(endpoint)
            all_endpoints.append(endpoint)
    coverage = [0 for _ in range(max_endpoint + 1)]
    for start, end in all_endpoints:
        for x in range(start, end):
            coverage[x] += 1
    window_len = 100
    cutpoint = -1
    for i in range(len(coverage) - window_len):
        if side == "in":
            window_start = (len(coverage) - window_len) - i
            window_end = len(coverage) - i
            if window_start < 0:
                window_start = 0
            if window_end > len(coverage):
                window_end = len(coverage)
            avg_cov = _mean(coverage[window_start:window_end])
            if avg_cov >= MIN_EDGE_COV:
                cutpoint = window_end
                break
        elif side == "out":
            window_start = i
            window_end = i + window_len
            if window_start < 0:
                window_start = 0
            if window_end > len(coverage):
                window_end = len(coverage)
            avg_cov = _mean(coverage[window_start:window_end])
            if avg_cov >= MIN_EDGE_COV:
                cutpoint = window_start
                break
    return cutpoint


def truncate_consensus(side, cutpoint, cons_al_file, template,
                       polished_consensus, cut_cons_file):
    if cutpoint == -1:
        logger.debug("No cutpoint for consensus file")
        return

    CONS_ALN_RATE = trestle_config.vals["cons_aln_rate"]
    cons_al = _read_alignment(cons_al_file, template, CONS_ALN_RATE)
    consensus_endpoint = -1
    if cons_al and cons_al[0]:
        consensus_endpoint = _find_consensus_endpoint(cutpoint, cons_al, side)
    else:
        logger.debug("No cons alignment to template, no cut consensus")
        return

    if consensus_endpoint != -1:
        cons_seqs = fp.read_sequence_dict(polished_consensus)
        cons_head = list(cons_seqs.keys())[0]
        consensus = list(cons_seqs.values())[0]
        if side == "in":
            start = 0
            end = consensus_endpoint
        elif side == "out":
            start = consensus_endpoint
            end = len(consensus)
        cut_head = "".join([cons_head, "|{0}_{1}".format(start, end)])
        cut_dict = {cut_head:consensus[start:end]}
        fp.write_fasta_dict(cut_dict, cut_cons_file)


def _find_consensus_endpoint(cutpoint, aligns, side):
    consensus_endpoint = -1
    #first try collapsing
    coll_aln = _collapse_cons_aln(aligns)
    if cutpoint >= coll_aln.trg_start and cutpoint < coll_aln.trg_end:
        trg_aln, _ = _index_mapping(coll_aln.trg_seq)
        _, aln_qry = _index_mapping(coll_aln.qry_seq)
        cutpoint_minus_start = cutpoint - coll_aln.trg_start
        aln_ind = trg_aln[cutpoint_minus_start]
        qry_ind = aln_qry[aln_ind]
        consensus_endpoint = qry_ind + coll_aln.qry_start
    else:
        #otherwise try each alignment
        MIN_SUPP_ALN_LEN = trestle_config.vals["min_supp_align_len"]
        #save tuples of cutpoint distance, cutpoint
        aln_endpoints = []
        for i, aln in enumerate(aligns[0]):
            if i == 0 or len(aln.trg_seq) >= MIN_SUPP_ALN_LEN:
                if cutpoint >= aln.trg_start and cutpoint < aln.trg_end:
                    trg_aln, _ = _index_mapping(aln.trg_seq)
                    _, aln_qry = _index_mapping(aln.qry_seq)
                    cutpoint_minus_start = cutpoint - aln.trg_start
                    if cutpoint_minus_start < 0:
                        logger.warning("%s %s %s %s %s", aln.qry_id, aln.trg_id, side,
                                       cutpoint, cutpoint_minus_start)
                        aln_ind = trg_aln[0]
                    elif cutpoint_minus_start >= len(trg_aln):
                        logger.warning("%s %s %s %s %s", aln.qry_id, aln.trg_id, side,
                                       cutpoint, cutpoint_minus_start)
                        aln_ind = trg_aln[-1]
                    else:
                        aln_ind = trg_aln[cutpoint_minus_start]
                    qry_ind = aln_qry[aln_ind]
                    endpoint = qry_ind + coll_aln.qry_start
                    aln_endpoints.append((0, endpoint))
                elif side == "in" and cutpoint >= aln.trg_end:
                    endpoint = aln.qry_end
                    distance = cutpoint - aln.trg_end
                    aln_endpoints.append((distance, endpoint))
                elif side == "out" and cutpoint < aln.trg_start:
                    endpoint = aln.qry_start
                    distance = aln.trg_start - cutpoint
                    aln_endpoints.append((distance, endpoint))
        if aln_endpoints:
            consensus_endpoint = sorted(aln_endpoints)[0][1]
    return consensus_endpoint

#Partition Reads Functions


def partition_reads(edges, it, side, position_path, cons_align_path,
                    template, read_align_path, consensuses,
                    confirmed_pos_path, part_file,
                    headers_to_id):
    CONS_ALN_RATE = trestle_config.vals["cons_aln_rate"]
    BUFFER_COUNT = trestle_config.vals["buffer_count"]

    skip_bool = False
    _, pos = div.read_positions(position_path)

    cons_aligns = {}
    for edge_id in edges:
        if not os.path.isfile(cons_align_path.format(it, side, edge_id)):
            skip_bool = True
        else:
            cons_aligns[edge_id] = _read_alignment(cons_align_path.format(it,
                                                                          side,
                                                                          edge_id),
                                                   template,
                                                   CONS_ALN_RATE)
        if (skip_bool or
                not cons_aligns or
                not cons_aligns[edge_id] or
                not cons_aligns[edge_id][0]):
            logger.debug("No cons alignment found for edge %s", edge_id)
            skip_bool = True
    if skip_bool:
        if it <= 1:
            confirmed_pos = {"total":[], "sub":[], "ins":[], "del":[]}
            rejected_pos = {"total":[], "sub":[], "ins":[], "del":[]}
            consensus_pos = pos
        else:
            previous_pos = _read_confirmed_positions(
                                confirmed_pos_path.format(it - 1, side))
            confirmed_pos, rejected_pos, consensus_pos = previous_pos
    else:
        #collapse multiple consensus alignments to the template
        coll_cons_aligns = {}
        for edge_id in cons_aligns:
            aln = cons_aligns[edge_id]
            coll_cons_aligns[edge_id] = _collapse_cons_aln(aln)
        curr_pos = _evaluate_positions(pos, coll_cons_aligns, side)
        confirmed_pos, rejected_pos, consensus_pos = curr_pos
    _write_confirmed_positions(confirmed_pos, rejected_pos, pos,
                               confirmed_pos_path.format(it, side))
    read_aligns = {}
    for edge_id in edges:
        if (not os.path.isfile(read_align_path.format(it, side, edge_id)) or
            not os.path.isfile(consensuses[(it, side, edge_id)])):
            skip_bool = True
        elif not skip_bool:
            read_aligns[edge_id] = _read_alignment(
                                        read_align_path.format(it, side,
                                                               edge_id),
                                        consensuses[(it, side, edge_id)],
                                        CONS_ALN_RATE)
        if (skip_bool or
                not read_aligns or
                not read_aligns[edge_id] or
                not read_aligns[edge_id][0]):
            logger.debug("No read alignment found for edge %s", edge_id)
            skip_bool = True
    if skip_bool:
        partitioning = _read_partitioning_file(part_file.format(it - 1, side))
    else:
        partitioning = _classify_reads(read_aligns, consensus_pos,
                                       headers_to_id, BUFFER_COUNT)
    _write_partitioning_file(partitioning, part_file.format(it, side))


def _read_alignment(alignment, target_path, min_aln_rate):
    alignments = []
    fasta_dict = fp.read_sequence_dict(target_path)
    aln_reader = SynchronizedSamReader(alignment, fasta_dict, multiprocessing.Manager(),
                                       config.vals["max_read_coverage"])
    for ctg_id in fasta_dict:
        ctg_aln = aln_reader.get_alignments(ctg_id)
        if len(ctg_aln) == 0:
            continue
        alignments.append(ctg_aln)

    return alignments


def _collapse_cons_aln(cons_aligns):
    MAX_SUPP_ALIGN_OVERLAP = trestle_config.vals["max_supp_align_overlap"]
    coll_aln = None
    for aln in cons_aligns[0]:
        if coll_aln is None:
            coll_aln = aln
        elif _overlap(coll_aln, aln) <= MAX_SUPP_ALIGN_OVERLAP:
            coll_aln = _collapse(coll_aln, aln)
    return coll_aln


def _overlap(aln_one, aln_two):
    qry_overlap_lens = []
    if (aln_one.qry_start >= aln_two.qry_start and
            aln_one.qry_start < aln_two.qry_end):
        if aln_one.qry_end >= aln_two.qry_end:
            qry_overlap_lens.append(aln_two.qry_end - aln_one.qry_start)
        else:
            qry_overlap_lens.append(aln_one.qry_end - aln_one.qry_start)
    if (aln_one.qry_end > aln_two.qry_start and
            aln_one.qry_end <= aln_two.qry_end):
        if aln_one.qry_start <= aln_two.qry_start:
            qry_overlap_lens.append(aln_one.qry_end - aln_two.qry_start)
        else:
            qry_overlap_lens.append(aln_one.qry_end - aln_one.qry_start)
    if (aln_two.qry_start >= aln_one.qry_start and
            aln_two.qry_start < aln_one.qry_end):
        if aln_two.qry_end >= aln_one.qry_end:
            qry_overlap_lens.append(aln_one.qry_end - aln_two.qry_start)
        else:
            qry_overlap_lens.append(aln_two.qry_end - aln_two.qry_start)
    if (aln_two.qry_end > aln_one.qry_start and
            aln_two.qry_end <= aln_one.qry_end):
        if aln_two.qry_start <= aln_one.qry_start:
            qry_overlap_lens.append(aln_two.qry_end - aln_one.qry_start)
        else:
            qry_overlap_lens.append(aln_two.qry_end - aln_two.qry_start)
    qry_len = 0
    if qry_overlap_lens:
        qry_len = min(qry_overlap_lens)
    trg_overlap_lens = []
    if (aln_one.trg_start >= aln_two.trg_start and
            aln_one.trg_start < aln_two.trg_end):
        if aln_one.trg_end >= aln_two.trg_end:
            trg_overlap_lens.append(aln_two.trg_end - aln_one.trg_start)
        else:
            trg_overlap_lens.append(aln_one.trg_end - aln_one.trg_start)
    if (aln_one.trg_end > aln_two.trg_start and
            aln_one.trg_end <= aln_two.trg_end):
        if aln_one.trg_start <= aln_two.trg_start:
            trg_overlap_lens.append(aln_one.trg_end - aln_two.trg_start)
        else:
            trg_overlap_lens.append(aln_one.trg_end - aln_one.trg_start)
    if (aln_two.trg_start >= aln_one.trg_start and
            aln_two.trg_start < aln_one.trg_end):
        if aln_two.trg_end >= aln_one.trg_end:
            trg_overlap_lens.append(aln_one.trg_end - aln_two.trg_start)
        else:
            trg_overlap_lens.append(aln_two.trg_end - aln_two.trg_start)
    if (aln_two.trg_end > aln_one.trg_start and
            aln_two.trg_end <= aln_one.trg_end):
        if aln_two.trg_start <= aln_one.trg_start:
            trg_overlap_lens.append(aln_two.trg_end - aln_one.trg_start)
        else:
            trg_overlap_lens.append(aln_two.trg_end - aln_two.trg_start)
    trg_len = 0
    if trg_overlap_lens:
        trg_len = min(trg_overlap_lens)
    return max([qry_len, trg_len])


def _collapse(aln_one, aln_two):
    MAX_SUPP_ALIGN_OVERLAP = trestle_config.vals["max_supp_align_overlap"]
    out_aln = copy.deepcopy(aln_one)
    if (aln_one.qry_sign == "-" or aln_two.qry_sign == "-" or
            _overlap(aln_one, aln_two) > MAX_SUPP_ALIGN_OVERLAP):
        return out_aln
    if (aln_one.qry_start <= aln_two.qry_start and
            aln_one.trg_start <= aln_two.trg_start):
        qry_merge_outs = _merge_alns(aln_one.qry_start, aln_one.qry_end,
                                     aln_one.qry_seq, aln_two.qry_start,
                                     aln_two.qry_end, aln_two.qry_seq)
        one_qry_seq, two_qry_seq, out_qry_end = qry_merge_outs
        trg_merge_outs = _merge_alns(aln_one.trg_start, aln_one.trg_end,
                                     aln_one.trg_seq, aln_two.trg_start,
                                     aln_two.trg_end, aln_two.trg_seq)
        one_trg_seq, two_trg_seq, out_trg_end = trg_merge_outs
        fill_qry = ""
        fill_trg = ""
        qry_lens = len(one_qry_seq) + len(two_qry_seq)
        trg_lens = len(one_trg_seq) + len(two_trg_seq)
        if qry_lens > trg_lens:
            diff = qry_lens - trg_lens
            fill_trg = "-" * diff
        elif trg_lens > qry_lens:
            diff = trg_lens - qry_lens
            fill_qry = "-" * diff
        out_qry_seq = "".join([one_qry_seq, fill_qry, two_qry_seq])
        out_trg_seq = "".join([one_trg_seq, fill_trg, two_trg_seq])
        out_err_rate = ((aln_one.err_rate * len(aln_one.trg_seq) +
                         aln_two.err_rate * len(aln_two.trg_seq)) /
                         (len(aln_one.trg_seq) + len(aln_two.trg_seq)))
        out_aln = Alignment(aln_one.qry_id, aln_one.trg_id, aln_one.qry_start,
                            out_qry_end, aln_one.qry_sign, aln_one.qry_len,
                            aln_one.trg_start, out_trg_end, aln_one.trg_sign,
                            aln_one.trg_len, out_qry_seq, out_trg_seq,
                            out_err_rate, is_secondary=False, is_supplementary=False, map_qv=0)
        return out_aln
    elif (aln_two.qry_start <= aln_one.qry_start and
            aln_two.trg_start <= aln_one.trg_start):
        qry_merge_outs = _merge_alns(aln_two.qry_start, aln_two.qry_end,
                                     aln_two.qry_seq, aln_one.qry_start,
                                     aln_one.qry_end, aln_one.qry_seq)
        two_qry_seq, one_qry_seq, out_qry_end = qry_merge_outs
        trg_merge_outs = _merge_alns(aln_two.trg_start, aln_two.trg_end,
                                     aln_two.trg_seq, aln_one.trg_start,
                                     aln_one.trg_end, aln_one.trg_seq)
        two_trg_seq, one_trg_seq, out_trg_end = trg_merge_outs
        fill_qry = ""
        fill_trg = ""
        qry_lens = len(two_qry_seq) + len(one_qry_seq)
        trg_lens = len(two_trg_seq) + len(one_trg_seq)
        if qry_lens > trg_lens:
            diff = qry_lens - trg_lens
            fill_trg = "-" * diff
        elif trg_lens > qry_lens:
            diff = trg_lens - qry_lens
            fill_qry = "-" * diff
        out_qry_seq = "".join([two_qry_seq, fill_qry, one_qry_seq])
        out_trg_seq = "".join([two_trg_seq, fill_trg, one_trg_seq])
        out_err_rate = ((aln_one.err_rate * len(aln_one.trg_seq) +
                         aln_two.err_rate * len(aln_two.trg_seq)) /
                         (len(aln_one.trg_seq) + len(aln_two.trg_seq)))
        out_aln = Alignment(aln_one.qry_id, aln_one.trg_id, aln_two.qry_start,
                            out_qry_end, aln_one.qry_sign, aln_one.qry_len,
                            aln_two.trg_start, out_trg_end, aln_one.trg_sign,
                            aln_one.trg_len, out_qry_seq, out_trg_seq,
                            out_err_rate, is_secondary=False, is_supplementary=False, map_qv=0)
        return out_aln
    return out_aln


def _merge_alns(first_start, first_end, first_seq,
                second_start, second_end, second_seq):
    first_out_seq = first_seq
    second_out_seq = second_seq
    out_end = second_end
    if first_end <= second_start:
        fill_qry_seq = "N" * (second_start - first_end)
        first_out_seq = "".join([first_seq, fill_qry_seq])
        second_out_seq = second_seq
    else:
        if first_end < second_end:
            overlap = first_end - second_start
            two_cut_ind = _overlap_to_aln_ind(overlap, second_seq)
            first_out_seq = first_seq
            second_out_seq = second_seq[two_cut_ind:]
        else:
            first_out_seq = first_seq
            second_out_seq = ""
            out_end = first_end
    return first_out_seq, second_out_seq, out_end


def _overlap_to_aln_ind(overlap, aln):
    num_bases = 0
    for i, base in enumerate(aln):
        if base != "-":
            num_bases += 1
        if num_bases == overlap:
            return i + 1
    return len(aln)


class EdgeAlignment(object):
    __slots__ = ("edge_id", "qry_seq", "trg_seq", "qry_start", "trg_start",
                 "trg_end", "in_alignment", "curr_aln_ind", "curr_qry_ind",
                 "curr_qry_nuc", "curr_trg_nuc", "curr_ins_nuc")

    def __init__(self, edge_id, qry_seq, trg_seq, qry_start, trg_start,
                 trg_end):
        self.edge_id = edge_id
        self.qry_seq = flye_aln.shift_gaps(trg_seq, qry_seq)
        self.trg_seq = flye_aln.shift_gaps(self.qry_seq, trg_seq)
        self.qry_start = qry_start
        self.trg_start = trg_start
        self.trg_end = trg_end
        self.in_alignment = False
        self.curr_aln_ind = -1
        self.curr_qry_ind = -1
        self.curr_qry_nuc = ""
        self.curr_trg_nuc = ""
        self.curr_ins_nuc = ""

    def reset_nucs(self):
        self.curr_qry_nuc = ""
        self.curr_trg_nuc = ""
        self.curr_ins_nuc = ""


def _evaluate_positions(pos, cons_aligns, side):
    #Includes insertions!
    confirmed_pos = {"total":[], "sub":[], "ins":[], "del":[]}
    rejected_pos = {"total":[], "sub":[], "ins":[], "del":[]}
    consensus_pos = {e:[] for e in cons_aligns}

    alns = {}
    for edge_id in cons_aligns:
        orig_aln = cons_aligns[edge_id]
        alns[edge_id] = EdgeAlignment(edge_id, orig_aln.qry_seq,
                                      orig_aln.trg_seq, orig_aln.qry_start,
                                      orig_aln.trg_start, orig_aln.trg_end)

    min_start_edge = min([alns[e].trg_start for e in alns])
    max_end_edge = max([alns[e].trg_end for e in alns])
    #end indices for conservatively defining confirmed positions
    min_end_edge = min([alns[e].trg_end for e in alns])
    max_start_edge = max([alns[e].trg_start for e in alns])

    for trg_ind in range(min_start_edge, max_end_edge):
        for edge_id in alns:
            aln = alns[edge_id]
            if aln.trg_start == trg_ind:
                aln.curr_aln_ind = 0
                aln.curr_qry_ind = aln.qry_start
                aln.in_alignment = True

            if aln.trg_start > trg_ind or aln.trg_end <= trg_ind:
                aln.in_alignment = False

            if aln.in_alignment:
                while aln.trg_seq[aln.curr_aln_ind] == "-":
                    if aln.qry_seq[aln.curr_aln_ind] != "-":
                        aln.curr_ins_nuc += aln.qry_seq[aln.curr_aln_ind]
                        aln.curr_qry_ind += 1
                    aln.curr_aln_ind += 1
                aln.curr_qry_nuc = aln.qry_seq[aln.curr_aln_ind]
                aln.curr_trg_nuc = aln.trg_seq[aln.curr_aln_ind]

        if trg_ind in pos["total"]:
            if ((side == "in" and trg_ind < min_end_edge) or
                (side == "out" and trg_ind >= max_start_edge)):
                ins_confirmed = False
                del_confirmed = False
                sub_confirmed = False
                qry_nuc = ""
                trg_nuc = ""
                for edge_id in alns:
                    aln = alns[edge_id]
                    if aln.in_alignment:
                        #Directly add positions only to consensuses 
                        # where insertions occur
                        #Add the position prior to curr_qry_ind to 
                        # account for insertion
                        if aln.curr_ins_nuc:
                            ins_confirmed = True
                            consensus_pos[edge_id].append(aln.curr_qry_ind - 1)

                        if qry_nuc and qry_nuc != aln.curr_qry_nuc:
                            if qry_nuc != "N" and aln.curr_qry_nuc != "N":
                                if qry_nuc == "-":
                                    del_confirmed = True
                                else:
                                    sub_confirmed = True
                        else:
                            qry_nuc = aln.curr_qry_nuc
                        if (trg_nuc and trg_nuc != aln.curr_trg_nuc and
                                trg_nuc != "N" and aln.curr_trg_nuc != "N"):
                            logger.debug("Inconsistent trg_nuc, %s %s %s %s",
                                         edge_id, trg_ind, trg_nuc,
                                         aln.curr_trg_nuc)
                        trg_nuc = aln.curr_trg_nuc
                if ins_confirmed or del_confirmed or sub_confirmed:
                    confirmed_pos["total"].append(trg_ind)
                    #Add positions to consensuses for only subs/deletions
                    if del_confirmed or sub_confirmed:
                        for edge_id in alns:
                            aln = alns[edge_id]
                            if aln.in_alignment:
                                consensus_pos[edge_id].append(aln.curr_qry_ind)
                    if trg_ind in pos["ins"]:
                        if ins_confirmed:
                            confirmed_pos["ins"].append(trg_ind)
                        else:
                            rejected_pos["ins"].append(trg_ind)
                    if trg_ind in pos["del"]:
                        if del_confirmed:
                            confirmed_pos["del"].append(trg_ind)
                        else:
                            rejected_pos["del"].append(trg_ind)
                    if trg_ind in pos["sub"]:
                        if sub_confirmed:
                            confirmed_pos["sub"].append(trg_ind)
                        else:
                            rejected_pos["sub"].append(trg_ind)
                else:
                    rejected_pos["total"].append(trg_ind)
                    if trg_ind in pos["ins"]:
                        rejected_pos["ins"].append(trg_ind)
                    if trg_ind in pos["del"]:
                        rejected_pos["del"].append(trg_ind)
                    if trg_ind in pos["sub"]:
                        rejected_pos["sub"].append(trg_ind)

        for edge_id in alns:
            aln = alns[edge_id]
            if aln.in_alignment:
                if aln.qry_seq[aln.curr_aln_ind] != "-":
                    aln.curr_qry_ind += 1
                aln.curr_aln_ind += 1

                aln.reset_nucs()

    return confirmed_pos, rejected_pos, consensus_pos


def _write_confirmed_positions(confirmed, rejected, pos, out_file):
    with open(out_file, 'w') as f:
        f.write(">Confirmed_total_positions_{0}\n"
                .format(len(confirmed["total"])))
        f.write(",".join([str(x) for x in sorted(confirmed["total"])]) + "\n")
        f.write(">Confirmed_sub_positions_{0}\n".format(len(confirmed["sub"])))
        f.write(",".join([str(x) for x in sorted(confirmed["sub"])]) + "\n")
        f.write(">Confirmed_del_positions_{0}\n".format(len(confirmed["del"])))
        f.write(",".join([str(x) for x in sorted(confirmed["del"])]) + "\n")
        f.write(">Confirmed_ins_positions_{0}\n".format(len(confirmed["ins"])))
        f.write(",".join([str(x) for x in sorted(confirmed["ins"])]) + "\n")
        f.write(">Rejected_total_positions_{0}\n".format(len(rejected["total"])))
        f.write(",".join([str(x) for x in sorted(rejected["total"])]) + "\n")
        f.write(">Rejected_sub_positions_{0}\n".format(len(rejected["sub"])))
        f.write(",".join([str(x) for x in sorted(rejected["sub"])]) + "\n")
        f.write(">Rejected_del_positions_{0}\n".format(len(rejected["del"])))
        f.write(",".join([str(x) for x in sorted(rejected["del"])])+ "\n")
        f.write(">Rejected_ins_positions_{0}\n".format(len(rejected["ins"])))
        f.write(",".join([str(x) for x in sorted(rejected["ins"])]) + "\n")
        f.write(">Tentative_total_positions_{0}\n".format(len(pos["total"])))
        f.write(",".join([str(x) for x in sorted(pos["total"])]) + "\n")
        f.write(">Tentative_sub_positions_{0}\n".format(len(pos["sub"])))
        f.write(",".join([str(x) for x in sorted(pos["sub"])]) + "\n")
        f.write(">Tentative_del_positions_{0}\n".format(len(pos["del"])))
        f.write(",".join([str(x) for x in sorted(pos["del"])]) + "\n")
        f.write(">Tentative_ins_positions_{0}\n".format(len(pos["ins"])))
        f.write(",".join([str(x) for x in sorted(pos["ins"])]) + "\n")


def _read_confirmed_positions(confirmed_file):
    confirmed = {"total":[], "sub":[], "ins":[], "del":[]}
    rejected = {"total":[], "sub":[], "ins":[], "del":[]}
    pos = {"total":[], "sub":[], "ins":[], "del":[]}
    with open(confirmed_file, "r") as f:
        for i, line in enumerate(f):
            line = line.strip()
            if i == 1 and line:
                confirmed["total"] = [int(x) for x in line.split(",")]
            elif i == 3 and line:
                confirmed["sub"] = [int(x) for x in line.split(",")]
            elif i == 5 and line:
                confirmed["del"] = [int(x) for x in line.split(",")]
            elif i == 7 and line:
                confirmed["ins"] = [int(x) for x in line.split(",")]
            elif i == 9 and line:
                rejected["total"] = [int(x) for x in line.split(",")]
            elif i == 11 and line:
                rejected["sub"] = [int(x) for x in line.split(",")]
            elif i == 13 and line:
                rejected["del"] = [int(x) for x in line.split(",")]
            elif i == 15 and line:
                rejected["ins"] = [int(x) for x in line.split(",")]
            elif i == 17 and line:
                pos["total"] = [int(x) for x in line.split(",")]
            elif i == 19 and line:
                pos["sub"] = [int(x) for x in line.split(",")]
            elif i == 21 and line:
                pos["del"] = [int(x) for x in line.split(",")]
            elif i == 23 and line:
                pos["ins"] = [int(x) for x in line.split(",")]
    return confirmed, rejected, pos


def _classify_reads(read_aligns, consensus_pos,
                    headers_to_id, buffer_count):
    #Includes insertion positions where an insertion occurs right before the 
    #position for the read.
    #partitioning format same as above:
    #list of (read_id, status, edge_id, top_score, total_score, header)
    partitioning = []

    read_scores = {}
    for edge_id in read_aligns:
        read_counts = {}
        for aln in read_aligns[edge_id][0]:
            read_header = aln.qry_id
            cons_header = aln.trg_id
            #Unmapped segments will not be scored
            if cons_header == "*":
                continue
            if read_header not in read_scores:
                read_scores[read_header] = {}
            read_scores[read_header][edge_id] = 0
            if read_header not in read_counts:
                read_counts[read_header] = 1
            else:
                read_counts[read_header] += 1
            #Any alignments after the first supplementary will not be scored
            if read_counts[read_header] > 2:
                continue
            positions = consensus_pos[edge_id]
            trg_aln, _ = _index_mapping(aln.trg_seq)
            for pos in positions:
                if pos >= aln.trg_start and pos < aln.trg_end:
                    pos_minus_start = pos - aln.trg_start
                    aln_ind = trg_aln[pos_minus_start]
                    if aln.qry_seq[aln_ind] == aln.trg_seq[aln_ind]:
                        read_scores[read_header][edge_id] += 1
    #Iterate through all read_headers so partitioning will be a complete set
    for read_header in headers_to_id:
        read_id = headers_to_id[read_header]
        if read_header in read_scores:
            tie_bool = False
            top_edge = 0
            top_score = 0
            total_score = 0
            for edge_id in read_scores[read_header]:
                edge_score = read_scores[read_header][edge_id]
                #print edge_id, edge_score, top_score
                if edge_score - buffer_count > top_score:
                    top_edge = edge_id
                    top_score = edge_score
                    tie_bool = False
                elif (edge_score - buffer_count <= top_score and
                      edge_score >= top_score):
                    top_score = edge_score
                    tie_bool = True
                elif (edge_score >= top_score - buffer_count and
                      edge_score < top_score):
                    tie_bool = True
                total_score += edge_score

            if total_score == 0:
                status_label = "None"
                edge_label = "NA"
            elif tie_bool:
                status_label = "Tied"
                edge_label = "NA"
            else:
                status_label = "Partitioned"
                edge_label = str(top_edge)
            partitioning.append((read_id, status_label, edge_label,
                                 top_score, total_score, read_header))
        else:
            status_label = "None"
            edge_label = "NA"
            top_score = 0
            total_score = 0
            partitioning.append((read_id, status_label, edge_label,
                                top_score, total_score, read_header))
    return partitioning


def _index_mapping(aln):
    #Given a genomic index, return the alignment index of the alignment
    al_inds = []
    #Given an alignment index, return the genomic index at that position
    gen_inds = []
    for i,b in enumerate(aln):
        gen_inds.append(len(al_inds))
        if b != '-':
            al_inds.append(i)
    return al_inds, gen_inds


def init_side_stats(rep, side, repeat_edges, min_overlap, position_path,
                    partitioning, prev_parts, template_len, stats_file):
    SUB_THRESH = trestle_config.vals["sub_thresh"]
    DEL_THRESH = trestle_config.vals["del_thresh"]
    INS_THRESH = trestle_config.vals["ins_thresh"]
    FLANKING_LEN = trestle_config.vals["flanking_len"]
    BUFFER_COUNT = trestle_config.vals["buffer_count"]
    MAX_ITER = trestle_config.vals["max_iter"]
    MIN_EDGE_COV = trestle_config.vals["min_edge_cov"]
    CONS_ALN_RATE = trestle_config.vals["cons_aln_rate"]

    _, pos = div.read_positions(position_path)
    #Count partitioned reads
    edge_below_cov = False
    part_list = _read_partitioning_file(partitioning)
    edge_reads, _, _ = _get_partitioning_info(part_list,
                                              repeat_edges[rep][side])
    #Check break condition for iteration loop
    for edge in repeat_edges[rep][side]:
        if edge_reads[edge] < MIN_EDGE_COV:
            edge_below_cov = True
    prev_parts.add(tuple(part_list))
    #Prepare header for iteration stats
    #Iter,Rep Lens,Confirmed/Rejected Pos,Partitioned Reads
    header_labels = ["Iter"]
    for edge in sorted(repeat_edges[rep][side]):
        header_labels.extend(["Rep Len {0}".format(edge)])
    header_labels.extend(["Confirmed Pos", "Rejected Pos"])
    for edge in sorted(repeat_edges[rep][side]):
        header_labels.extend(["#Reads {0}".format(edge)])
    header_labels.extend(["#Tied", "#Unassigned"])
    spaced_header = ["{:11}".format(h) for h in header_labels]
    #Write stats output
    with open(stats_file, 'w') as f:
        f.write("{0:25}\t{1}\n".format("Repeat:", rep))
        f.write("{0:25}\t'{1}'\n".format("Side:", side))
        f.write("{0:25}\t".format("Edges:"))
        f.write(", ".join([str(x) for x in sorted(repeat_edges[rep][side])]) + "\n")
        f.write("{0:25}\t{1}\n\n".format("Template Length:", template_len))
        f.write("Initial Option Values\n")
        f.write("{0:25}\t{1}\n".format("min_overlap:", min_overlap))
        f.write("{0:25}\t{1}\n".format("sub_thresh:", SUB_THRESH))
        f.write("{0:25}\t{1}\n".format("del_thresh:", DEL_THRESH))
        f.write("{0:25}\t{1}\n".format("ins_thresh:", INS_THRESH))
        f.write("{0:25}\t{1}\n".format("flanking_len:", FLANKING_LEN))
        f.write("{0:25}\t{1}\n".format("buffer_count:", BUFFER_COUNT))
        f.write("{0:25}\t{1}\n".format("max_iter:", MAX_ITER))
        f.write("{0:25}\t{1}\n".format("min_edge_cov:", MIN_EDGE_COV))
        f.write("{0:25}\t{1}\n".format("cons_aln_rate:", CONS_ALN_RATE))
        f.write("\n")
        f.write("The following numbers are calculated based on moving ")
        f.write("into the repeat from the '{0}' direction\n\n".format(side))
        f.write("{0}\n".format("Divergent Positions:"))
        f.write("{0:25}\t{1}\n".format("Total", len(pos["total"])))
        f.write("{0:25}\t{1}\n".format("Substitutions", len(pos["sub"])))
        f.write("{0:25}\t{1}\n".format("Deletions", len(pos["del"])))
        f.write("{0:25}\t{1}\n".format("Insertions", len(pos["ins"])))
        f.write("\n")
        f.write("{0:25}\t{1}\n".format("Total Starting Reads:",
                                       sum(edge_reads.values())))
        for edge in sorted(repeat_edges[rep][side]):
            f.write("{0}{1}{2:18}\t{3}\n".format("Edge ", edge,
                                                 " starting reads:",
                                                 edge_reads[edge]))
        f.write("\n\n")
        f.write("\t".join(spaced_header))
        f.write("\n")

    return edge_below_cov


def update_side_stats(edges, it, side, cons_align_path, template,
                      confirmed_pos_path, partitioning, prev_parts,
                      stats_file):
    CONS_ALN_RATE = trestle_config.vals["cons_aln_rate"]
    MIN_EDGE_COV = trestle_config.vals["min_edge_cov"]
    #Write stats for each iteration
    #Iter,Rep Lens,Confirmed/Rejected Pos,Partitioned Reads
    stats_out = [str(it)]
    for edge_id in sorted(edges):
        rep_len = 0
        if os.path.isfile(cons_align_path.format(it, side, edge_id)):
            cons_align = _read_alignment(cons_align_path.format(it, side,
                                                                edge_id),
                                         template,
                                         CONS_ALN_RATE)
            if cons_align and cons_align[0]:
                if side == "in":
                    rep_len = (cons_align[0][0].qry_len -
                                cons_align[0][0].qry_start)
                elif side == "out":
                    rep_len = cons_align[0][0].qry_end
        stats_out.extend([str(rep_len)])
    confirmed_total = 0
    rejected_total = 0
    if it > 0:
        confirmed, rejected, _ = _read_confirmed_positions(confirmed_pos_path)
        confirmed_total = len(confirmed["total"])
        rejected_total = len(rejected["total"])
    stats_out.extend([str(confirmed_total),
                      str(rejected_total)])
    edge_below_cov = False
    dup_part = False
    part_list = _read_partitioning_file(partitioning)
    edge_reads, tied_reads, unassigned_reads = _get_partitioning_info(part_list, edges)
    for edge_id in sorted(edges):
        stats_out.extend([str(edge_reads[edge_id])])
    stats_out.extend([str(tied_reads), str(unassigned_reads)])
    #Check break conditions for iteration loop
    for edge in edges:
        if edge_reads[edge] < MIN_EDGE_COV:
            edge_below_cov = True
    if tuple(part_list) in prev_parts:
        dup_part = True
    else:
        prev_parts.add(tuple(part_list))
    spaced_header = ["{:11}".format(x) for x in stats_out]
    with open(stats_file, "a") as f:
        f.write("\t".join(spaced_header))
        f.write("\n")

    return edge_below_cov, dup_part


def finalize_side_stats(edges, it, side, cons_align_path, template,
                 cons_vs_cons_path, consensuses, confirmed_pos_path,
                 partitioning, edge_below_cov, dup_part, term_bool,
                 stats_file):
    CONS_ALN_RATE = trestle_config.vals["cons_aln_rate"]
    MAX_ITER = trestle_config.vals["max_iter"]

    with open(stats_file, "a") as f:
        f.write("\n\n")
        f.write("{0:26}\t{1}\n\n".format("Final Iter:", it))
        f.write("Iteration terminated because:\n")
        if it == MAX_ITER:
            f.write("Max iter reached\n")
        if edge_below_cov:
            f.write("Edge coverage fell below min_edge_cov\n")
        if dup_part:
            f.write("Partitioning was identical to a previous iteration\n")
        if term_bool:
            f.write("Encountered empty consensus sequence or alignment\n")
        f.write("\n")
        #Write out alignment indices for edges vs template
        limit_ind = None
        limit_label = ""
        if side == "in":
            limit_label = "Min Template End"
        elif side == "out":
            limit_label = "Max Template Start"
        for edge_id in sorted(edges):
            qry_start = 0
            qry_end = 0
            qry_len = 0
            trg_start = 0
            trg_end = 0
            trg_len = 0
            curr_cons_path = cons_align_path.format(it, side, edge_id)
            if os.path.isfile(curr_cons_path):
                cons_align = _read_alignment(curr_cons_path,
                                             template,
                                             CONS_ALN_RATE)
                if cons_align and cons_align[0]:
                    #collapse multiple consensus alignments
                    coll_cons = _collapse_cons_aln(cons_align)
                    qry_start = coll_cons.qry_start
                    qry_end = coll_cons.qry_end
                    qry_len = coll_cons.qry_len
                    trg_start = coll_cons.trg_start
                    trg_end = coll_cons.trg_end
                    trg_len = coll_cons.trg_len
                    if limit_ind is None or (
                            (side == "in" and trg_end < limit_ind) or
                            (side == "out" and trg_start >= limit_ind)):
                        if side == "in":
                            limit_ind = trg_end
                        elif side == "out":
                            limit_ind = trg_start
            f.write("Edge {0}|Template Alignment\n".format(edge_id))
            f.write("{0}{1}{2:20}\t{3:5}-{4:5} of {5:5}\n".format(
                    "Edge ", edge_id, ":",
                    qry_start, qry_end, qry_len))
            f.write("{0:26}\t{1:5}-{2:5} of {3:5}\n".format("Template:",
                    trg_start, trg_end, trg_len))
        f.write("\n")
        f.write("{0:26}\t{1}\n".format(limit_label, limit_ind))
        f.write("(End of positions considered)\n\n")
        #Write out alignment indices for edges vs edges
        edge_pairs = sorted(combinations(edges, 2))
        for edge_one, edge_two in edge_pairs:
            qry_start = 0
            qry_end = 0
            qry_len = 0
            trg_start = 0
            trg_end = 0
            trg_len = 0
            qry_seq = ""
            trg_seq = ""
            if (os.path.isfile(cons_vs_cons_path.format(it, side, edge_one,
                                                       it, side, edge_two)) and
                os.path.isfile(consensuses[(it, side, edge_two)])):
                cons_vs_cons = _read_alignment(cons_vs_cons_path.format(
                                                    it, side, edge_one,
                                                    it, side, edge_two),
                                               consensuses[(it, side,
                                                            edge_two)],
                                               CONS_ALN_RATE)
                if cons_vs_cons and cons_vs_cons[0]:
                    qry_start = cons_vs_cons[0][0].qry_start
                    qry_end = cons_vs_cons[0][0].qry_end
                    qry_len = cons_vs_cons[0][0].qry_len
                    trg_start = cons_vs_cons[0][0].trg_start
                    trg_end = cons_vs_cons[0][0].trg_end
                    trg_len = cons_vs_cons[0][0].trg_len
                    qry_seq = cons_vs_cons[0][0].qry_seq
                    trg_seq = cons_vs_cons[0][0].trg_seq
            f.write("Edge {0}|Edge {1} Alignment\n".format(edge_one, edge_two))
            f.write("{0}{1}{2:20}\t{3:5}-{4:5} of {5:5}\n".format(
                    "Edge ", edge_one, ":",
                    qry_start, qry_end, qry_len))
            f.write("{0}{1}{2:20}\t{3:5}-{4:5} of {5:5}\n".format(
                    "Edge ", edge_two, ":",
                    trg_start, trg_end, trg_len))
            div_rate = _calculate_divergence(qry_seq, trg_seq)
            f.write("{0:26}\t{1:.4f}\n".format("Divergence Rate:", div_rate))
        f.write("\n")
        #Write overall position stats
        types = ["total", "sub", "del", "ins"]
        confirmed = {t:[] for t in types}
        rejected = {t:[] for t in types}
        pos = {t:[] for t in types}
        if it > 0:
            confirmed_pos_output = _read_confirmed_positions(confirmed_pos_path)
            confirmed, rejected, pos = confirmed_pos_output
        if side == "in":
            largest_pos = -1
            if confirmed["total"]:
                largest_pos = max(confirmed["total"])
            f.write("{0:26}\t{1}\n".format("Largest Confirmed Position:",
                                           largest_pos))
        elif side == "out":
            smallest_pos = -1
            if confirmed["total"]:
                smallest_pos = min(confirmed["total"])
            f.write("{0:26}\t{1}\n".format("Smallest Confirmed Position:",
                                           smallest_pos))
        remainings = {}
        for typ in types:
            remainings[typ] = len(pos[typ]) - (len(confirmed[typ]) +
                                               len(rejected[typ]))
        type_strs = ["Total", "Sub", "Del", "Ins"]
        for typ, typ_str in zip(types, type_strs):
            confirmed_frac = 0.0
            rejected_frac = 0.0
            remaining_frac = 0.0
            if len(pos[typ]) != 0:
                confirmed_frac = len(confirmed[typ]) / float(len(pos[typ]))
                rejected_frac = len(rejected[typ]) / float(len(pos[typ]))
                remaining_frac = remainings[typ] / float(len(pos[typ]))
            f.write("{0:26}\t{1}/{2} = {3:.3f}\n".format(
                            "Confirmed {0} Positions:".format(typ_str),
                            len(confirmed[typ]),
                            len(pos[typ]),
                            confirmed_frac))
            f.write("{0:26}\t{1}/{2} = {3:.3f}\n".format(
                            "Rejected {0} Positions:".format(typ_str),
                            len(rejected[typ]),
                            len(pos[typ]),
                            rejected_frac))
            f.write("{0:26}\t{1}/{2} = {3:.3f}\n".format(
                            "Remaining {0} Positions:".format(typ_str),
                            remainings[typ],
                            len(pos[typ]),
                            remaining_frac))
            f.write("\n")
        f.write("\n")
        #Write overall partitioning stats
        part_list = _read_partitioning_file(partitioning)
        edge_reads = {edge:0 for edge in edges}
        tied_reads = 0
        unassigned_reads = 0
        total_reads = len(part_list)
        for _, status, edge, _, _, _  in part_list:
            if status == "Partitioned" and edge != "NA":
                edge_reads[int(edge)] += 1
            elif status == "Tied":
                tied_reads += 1
            elif status == "None":
                unassigned_reads += 1
            else:
                exception_str = "Unknown status {0} in partitioning file {1}"
                raise Exception(exception_str.format(status, partitioning))
        for edge_id in sorted(edges):
            f.write("{0}{1}{2:13}\t{3}/{4} = {5:.4f}\n".format(
                                    "Total Edge ", edge_id, " Reads:",
                                    edge_reads[edge_id], total_reads,
                                    edge_reads[edge_id] / float(total_reads)))
        f.write("{0:26}\t{1}/{2} = {3:.4f}\n".format("Total Tied Reads:",
                                              tied_reads, total_reads,
                                              tied_reads / float(total_reads)))
        f.write("{0:26}\t{1}/{2} = {3:.4f}\n".format("Total Unassigned Reads:",
                                      unassigned_reads, total_reads,
                                      unassigned_reads / float(total_reads)))
        f.write("\n")


def init_int_stats(rep, repeat_edges, zero_it, position_path, partitioning,
                   all_reads_file, template_len, cov, int_stats_file):
    #Count edge reads
    side_reads = {}
    total_reads = 0
    all_side_reads = 0
    internal_reads = 0
    for side in sorted(repeat_edges[rep]):
        part_list = _read_partitioning_file(partitioning.format(zero_it, side))
        total_reads = len(part_list)
        partitioning_outputs = _get_partitioning_info(part_list,
                                                      repeat_edges[rep][side])
        side_reads[side], _, _ = partitioning_outputs
        all_side_reads += sum(side_reads[side].values())
    internal_reads = total_reads - all_side_reads
    all_reads_n50 = _n50(all_reads_file)
    #Prepare header for iterative integrated stats
    #in/out Iter,Mean in/out/gap Len,Confirmed/Rejected Pos,Bridging Reads
    header_labels = []
    for side in sorted(repeat_edges[rep]):
        header_labels.extend(["{0} Iter".format(side)])
    header_labels.extend(["in Len", "Gap Len", "out Len"])
    header_labels.extend(["Confirmed", "Rejected"])
    side_edges = []
    for side in sorted(repeat_edges[rep]):
        side_edges.append([])
        for edge in sorted(repeat_edges[rep][side]):
            side_edges[-1].append("{0}{1}".format(side,edge))
    for edge_pair in sorted(product(*side_edges)):
        header_labels.extend(["{0}".format("|".join(edge_pair))])
    spaced_header = ["{:8}".format(x) for x in header_labels]
    #Write to file
    with open(int_stats_file, 'w') as f:
        f.write("{0:16}\t{1}\n".format("Repeat:", rep))
        f.write("{0:16}\t{1}\n".format("Template Length:", template_len))
        f.write("{0:16}\t{1:.2f}\n".format("Avg Coverage:", cov))
        f.write("{0:16}\t{1}\n".format("# All Reads:", total_reads))
        f.write("{0:16}\t{1}\n\n".format("All Reads N50:", all_reads_n50))
        edge_headers = ["Side", "    Edge", "# Reads"]
        spaced_edge_header = ["{:5}".format(h) for h in edge_headers]
        f.write("\t".join(spaced_edge_header))
        f.write("\n")
        for side in sorted(repeat_edges[rep]):
            for edge_id in sorted(repeat_edges[rep][side]):
                edge_values = [side, edge_id, side_reads[side][edge_id]]
                spaced_values = ["{:6}".format(h) for h in edge_values]
                f.write("\t".join(spaced_values))
                f.write("\n")
        f.write("{0:12}\t  {1}\n".format("Internal", internal_reads))
        f.write("\n\n")
        f.write("\t".join(spaced_header))
        f.write("\n")


def update_int_stats(rep, repeat_edges, side_it, cons_align_path, template,
                     template_len, confirmed_pos_path, int_confirmed_path,
                     partitioning, int_stats_file):
    CONS_ALN_RATE = trestle_config.vals["cons_aln_rate"]

    stats_out = []
    #Add side iters
    for side in sorted(repeat_edges[rep]):
        stats_out.extend([str(side_it[side])])
    #Find median in, out, and gap lengths
    medians = {s:0 for s in repeat_edges[rep]}
    for side in sorted(repeat_edges[rep]):
        trg_limits = []
        for edge_id in sorted(repeat_edges[rep][side]):
            curr_cons_path = cons_align_path.format(side_it[side],
                                                    side, edge_id)
            if os.path.isfile(curr_cons_path):
                cons_align = _read_alignment(curr_cons_path,
                                             template,
                                             CONS_ALN_RATE)
                if cons_align and cons_align[0]:
                    if side == "in":
                        trg_limits.append(cons_align[0][0].trg_end)
                    elif side == "out":
                        trg_limits.append(template_len -
                                          cons_align[0][0].trg_start)
        if trg_limits:
            medians[side] = _get_median(trg_limits)
    gap_len = template_len - (medians["in"] + medians["out"])
    stats_out.extend([str(medians["in"]), str(gap_len), str(medians["out"])])
    #Add confirmed and rejected reads
    in_confirmed_path = confirmed_pos_path.format(side_it["in"], "in")
    out_confirmed_path = confirmed_pos_path.format(side_it["out"], "out")
    types = ["total", "sub", "del", "ins"]
    int_confirmed = {t:[] for t in types}
    int_rejected = {t:[] for t in types}
    pos = {t:[] for t in types}
    if side_it["in"] > 0 and side_it["out"] > 0:
        all_in_pos = _read_confirmed_positions(in_confirmed_path)
        all_out_pos = _read_confirmed_positions(out_confirmed_path)
        confirmed_pos_outputs = _integrate_confirmed_pos(all_in_pos,
                                                         all_out_pos)
        int_confirmed, int_rejected, pos = confirmed_pos_outputs
    elif side_it["in"] > 0:
        all_in_pos = _read_confirmed_positions(in_confirmed_path)
        int_confirmed, int_rejected, pos = all_in_pos
    elif side_it["out"] > 0:
        all_out_pos = _read_confirmed_positions(out_confirmed_path)
        int_confirmed, int_rejected, pos = all_out_pos
    _write_confirmed_positions(int_confirmed, int_rejected, pos,
                               int_confirmed_path.format(side_it["in"],
                                                         side_it["out"]))
    stats_out.extend([str(len(int_confirmed["total"])),
                      str(len(int_rejected["total"]))])
    #Get bridging reads for each pair of in/out edges
    side_headers_dict = {}
    all_headers = set()
    for side in sorted(repeat_edges[rep]):
        side_headers_dict[side] = {}
        part_list = _read_partitioning_file(partitioning.format(side_it[side],
                                                                side))
        for _, status, edge, _, _, header in part_list:
            all_headers.add(header)
            if status == "Partitioned" and edge != "NA":
                side_headers_dict[side][header] = (side, int(edge))
    bridging_reads = {}
    side_edges = []
    for side in sorted(repeat_edges[rep]):
        side_edges.append([])
        for edge in sorted(repeat_edges[rep][side]):
            side_edges[-1].append((side, edge))
    for edge_pair in sorted(product(*side_edges)):
        bridging_reads[edge_pair] = 0
    for header in all_headers:
        if (header in side_headers_dict["in"] and
            header in side_headers_dict["out"]):
            in_edge = side_headers_dict["in"][header]
            out_edge = side_headers_dict["out"][header]
            bridging_reads[(in_edge, out_edge)] += 1
    for edge_pair in sorted(bridging_reads):
        #stats_out.extend(["{0}".format(edge_pair)])
        stats_out.extend([str(bridging_reads[edge_pair])])
    spaced_header = ["{:8}".format(x) for x in stats_out]
    #Write to file
    with open(int_stats_file, "a") as f:
        f.write("\t".join(spaced_header))
        f.write("\n")


def finalize_int_stats(rep, repeat_edges, side_it, cons_align_path, template,
                       template_len, cons_vs_cons_path, consensuses,
                       int_confirmed_path, partitioning, int_stats_file,
                       resolved_seq_file):
    CONS_ALN_RATE = trestle_config.vals["cons_aln_rate"]
    MIN_BRIDGE_COUNT = trestle_config.vals["min_bridge_count"]
    MIN_BRIDGE_FACTOR = trestle_config.vals["min_bridge_factor"]

    #Resolved repeat seqs to be returned, NOT written
    resolved_repeats = {}
    summ_vals = []
    with open(int_stats_file, "a") as f:
        f.write("\n\n")
        for side in sorted(repeat_edges[rep]):
            f.write("{0}'{1}'{2:8}\t{3}\n"
                    .format("Final ", side, " Iter:", side_it[side]))
        f.write("\n\n")
        #Overall confirmed and rejected positions
        types = ["total", "sub", "del", "ins"]
        int_confirmed = {t:[] for t in types}
        int_rejected = {t:[] for t in types}
        pos = {t:[] for t in types}
        if side_it["in"] > 0 or side_it["out"] > 0:
            int_confirmed, int_rejected, pos = _read_confirmed_positions(
                int_confirmed_path.format(side_it["in"], side_it["out"]))
        remainings = {}
        for typ in types:
            remainings[typ] = len(pos[typ]) - (len(int_confirmed[typ]) +
                                               len(int_rejected[typ]))
        type_strs = ["Total", "Sub", "Del", "Ins"]
        for typ, typ_str in zip(types, type_strs):
            confirmed_frac = 0.0
            rejected_frac = 0.0
            remaining_frac = 0.0
            if len(pos[typ]) != 0:
                confirmed_frac = len(int_confirmed[typ]) / float(len(pos[typ]))
                rejected_frac = len(int_rejected[typ]) / float(len(pos[typ]))
                remaining_frac = remainings[typ] / float(len(pos[typ]))
            f.write("{0:26}\t{1}/{2} = {3:.3f}\n".format(
                            "Confirmed {0} Positions:".format(typ_str),
                            len(int_confirmed[typ]),
                            len(pos[typ]),
                            confirmed_frac))
            f.write("{0:26}\t{1}/{2} = {3:.3f}\n".format(
                            "Rejected {0} Positions:".format(typ_str),
                            len(int_rejected[typ]),
                            len(pos[typ]),
                            rejected_frac))
            f.write("{0:26}\t{1}/{2} = {3:.3f}\n".format(
                            "Remaining {0} Positions:".format(typ_str),
                            remainings[typ],
                            len(pos[typ]),
                            remaining_frac))
            f.write("\n")
        f.write("\n")
        #Basic stats for confirmed positions
        av_div = 0.0
        if template_len != 0:
            av_div = len(int_confirmed["total"]) / float(template_len)
        position_gaps = [0 for _ in range(len(int_confirmed["total"]) + 1)]
        curr_pos = 0
        for i, p in enumerate(int_confirmed["total"]):
            position_gaps[i] = p - curr_pos
            curr_pos = p
        position_gaps[-1] = template_len - curr_pos
        mean_position_gap = _mean(position_gaps)
        max_position_gap = max(position_gaps)
        f.write("{0:26}\t{1}\n".format("Template Length:", template_len))
        f.write("{0:26}\t{1}\n".format("# Confirmed Positions:",
                                       len(int_confirmed["total"])))
        f.write("{0:26}\t{1:.4f}\n".format("Confirmed Pos Avg Divergence:",
                                           av_div))
        f.write("{0:26}\t{1:.2f}\n".format("Mean Confirmed Pos Gap:",
                                           mean_position_gap))
        f.write("{0:26}\t{1}\n".format("Max Confirmed Pos Gap:",
                                       max_position_gap))
        f.write("\n\n")
        summ_vals.extend([len(int_confirmed["total"]), max_position_gap])
        #Write bridging reads
        side_headers_dict = {}
        all_headers = set()
        for side in sorted(repeat_edges[rep]):
            side_headers_dict[side] = {}
            part_list = _read_partitioning_file(partitioning.format(
                                                    side_it[side], side))
            for _, status, edge, _, _, header in part_list:
                all_headers.add(header)
                if status == "Partitioned" and edge != "NA":
                    side_headers_dict[side][header] = (side, int(edge))
        bridging_reads = {}
        side_edges = []
        for side in sorted(repeat_edges[rep]):
            side_edges.append([])
            for edge in sorted(repeat_edges[rep][side]):
                side_edges[-1].append((side, edge))
        for edge_pair in sorted(product(*side_edges)):
            bridging_reads[edge_pair] = 0
        for header in all_headers:
            if (header in side_headers_dict["in"] and
                header in side_headers_dict["out"]):
                in_edge = side_headers_dict["in"][header]
                out_edge = side_headers_dict["out"][header]
                bridging_reads[(in_edge, out_edge)] += 1
        for edge_pair in sorted(bridging_reads):
            pair_label = "|".join(["{0}{1}".format(x[0], x[1]) for x in edge_pair])
            f.write("{0}{1:21}\t{2}\n".format(pair_label, " Bridging Reads:",
                                              bridging_reads[edge_pair]))
        f.write("\n\n")
        #Write combos which are sets of bridging reads
        all_combos = _get_combos(side_edges[0], side_edges[1])
        combo_support = [0 for _ in all_combos]
        for i, combo in enumerate(all_combos):
            for edge_pair in combo:
                if edge_pair in bridging_reads:
                    combo_support[i] += bridging_reads[edge_pair]
        for i, combo in enumerate(all_combos):
            f.write("{0} {1}\n".format("Combo", i))
            coms = ["|".join(["".join([str(z) for z in x]) for x in y]) for y in combo]
            combo_edges = " + ".join(coms)
            f.write("{0:12}\t{1}\n".format("Resolution:", combo_edges))
            f.write("{0:12}\t{1}\n\n".format("Support:", combo_support[i]))
        #Bridging conditions 
        bridged = False
        bridged_edges = None
        combo_inds = list(zip(combo_support, list(range(len(combo_support)))))
        sorted_combos = sorted(combo_inds, reverse=True)
        if (len(sorted_combos) > 1 and
            sorted_combos[0][0] >= MIN_BRIDGE_COUNT and
            sorted_combos[0][0] >= sorted_combos[1][0] * MIN_BRIDGE_FACTOR):
            bridged = True
            bridged_edges = all_combos[sorted_combos[0][1]]
        best_combo = sorted_combos[0][1]
        best_support = sorted_combos[0][0]
        best_against = 0
        second_combo = -1
        second_support = 0
        if len(sorted_combos) > 1:
            for support, _ in sorted_combos[1:]:
                best_against += support
            second_combo = sorted_combos[1][1]
            second_support = sorted_combos[1][0]
        if bridged:
            f.write("BRIDGED\n")
            f.write("Bridging Combo: {0}\n".format(best_combo))
            br_ct_str = "{0} (min_bridge_count)".format(MIN_BRIDGE_COUNT)
            br_diff_str = "{0} * {1} (Combo {2} * min_bridge_factor)".format(
                second_support, MIN_BRIDGE_FACTOR, second_combo)
            f.write("Support = {0}\t> {1}\n{2:12}\t> {3}\n".format(
                best_support, br_ct_str, "", br_diff_str))
            f.write("Resolution:\n")
            for edge_pair in sorted(bridged_edges):
                f.write("{0[0]} {0[1]:2}  {1:3} {2[0]} {2[1]}\n"
                        .format(edge_pair[0], "->", edge_pair[1]))
            f.write("\n\n")
        else:
            f.write("UNBRIDGED\n")
            f.write("Best combo {0}\n".format(best_combo))
            f.write("{0:20}\t{1}\n".format("min_bridge_count",
                                           MIN_BRIDGE_COUNT))
            f.write("{0:20}\t{1}\n\n\n".format("min_bridge_factor",
                                               MIN_BRIDGE_FACTOR))
        summ_vals.extend([bridged, best_support, best_against])
        #If not bridged, find in/gap/out lengths and divergence rates
        if not bridged:
            #Write median in, out, and gap lengths
            side_lens = {s:0 for s in repeat_edges[rep]}
            for side in sorted(repeat_edges[rep]):
                trg_limits = []
                for edge_id in sorted(repeat_edges[rep][side]):
                    curr_cons_path = cons_align_path.format(side_it[side],
                                                             side, edge_id)
                    if os.path.isfile(curr_cons_path):
                        cons_align = _read_alignment(curr_cons_path,
                                                     template,
                                                     CONS_ALN_RATE)
                        if cons_align and cons_align[0]:
                            if side == "in":
                                trg_limits.append(cons_align[0][0].trg_end)
                            elif side == "out":
                                trg_limits.append(template_len -
                                                  cons_align[0][0].trg_start)
                if trg_limits:
                    side_lens[side] = _get_median(trg_limits)
            gap_len = template_len - (side_lens["in"] + side_lens["out"])
            f.write("{0:30}\t{1}\n".format("Median in Sequence Length:",
                                           side_lens["in"]))
            f.write("{0:30}\t{1}\n".format("Median out Sequence Length:",
                                           side_lens["out"]))
            f.write("{0:30}\t{1}\n\n".format("Median Gap/Overlap Length:",
                                             gap_len))

            #Write mean in and out divergence rates
            div_rates = {s:[] for s in repeat_edges[rep]}
            for side in sorted(repeat_edges[rep]):
                side_pairs = sorted(combinations(repeat_edges[rep][side], 2))
                for edge_one, edge_two in side_pairs:
                    cons_cons_file = cons_vs_cons_path.format(
                                        side_it[side], side, edge_one,
                                        side_it[side], side, edge_two)
                    if (os.path.isfile(cons_cons_file) and
                        os.path.isfile(consensuses[(side_it[side],
                                                    side, edge_two)])):
                        cons_vs_cons = _read_alignment(cons_cons_file,
                                                       consensuses[(side_it[side],
                                                       side, edge_two)],
                                                       CONS_ALN_RATE)
                        if cons_vs_cons and cons_vs_cons[0]:
                            edge_rate = _calculate_divergence(
                                            cons_vs_cons[0][0].qry_seq,
                                            cons_vs_cons[0][0].trg_seq)
                            div_rates[side].append(edge_rate)
            mean_in_div = 0.0
            if div_rates["in"]:
                mean_in_div = _mean(div_rates["in"])
            mean_out_div = 0.0
            if div_rates["out"]:
                mean_out_div = _mean(div_rates["out"])
            weighted_mean_div = 0.0
            if side_lens["in"] + side_lens["out"] != 0:
                weighted_mean_div = ((mean_in_div*side_lens["in"] +
                                     mean_out_div*side_lens["out"]) /
                                     float(side_lens["in"] + side_lens["out"]))
            f.write("{0:30}\t{1}\n".format("Mean in Divergence Rate:",
                                            mean_in_div))
            f.write("{0:30}\t{1}\n".format("Mean out Divergence Rate:",
                                            mean_out_div))
            f.write("{0:30}\t{1}\n\n".format("Weighted Mean Divergence Rate:",
                                          weighted_mean_div))
            res_str = "No resolution so no resolved file for repeat {0}\n\n"
            f.write(res_str.format(rep))
            #for i, edge in enumerate(sorted(repeat_edges[rep]["in"])):
                #header = "Repeat_{0}_unbridged_copy_{1}".format(rep, i)
                #resolved_repeats[header] = ""
                #seq_dict = {header:""}
                #fp.write_fasta_dict(seq_dict, resolved_seq_file.format(i))
            summ_vals.extend(["*", "*"])
        #If bridged, find overlap and construct repeat copy sequences
        else:
            #Find end of repeat as min/max of in/out cons_vs_cons alignments
            edge_limits = {}
            for side in sorted(repeat_edges[rep]):
                side_pairs = sorted(combinations(repeat_edges[rep][side], 2))
                for edge_one, edge_two in side_pairs:
                    cons_cons_file = cons_vs_cons_path.format(
                                        side_it[side], side, edge_one,
                                        side_it[side], side, edge_two)
                    if (os.path.isfile(cons_cons_file) and
                        os.path.isfile(consensuses[(side_it[side],
                                                    side, edge_two)])):
                        cons_vs_cons = _read_alignment(cons_cons_file,
                                                       consensuses[(side_it[side],
                                                       side, edge_two)],
                                                       CONS_ALN_RATE)
                        if cons_vs_cons and cons_vs_cons[0]:
                            #collapse multiple consensus alignments
                            coll_cons = _collapse_cons_aln(cons_vs_cons)
                            one_start = coll_cons.qry_start
                            one_end = coll_cons.qry_end
                            two_start = coll_cons.trg_start
                            two_end = coll_cons.trg_end
                            if side == "in":
                                if (side, edge_one) not in edge_limits:
                                    edge_limits[(side, edge_one)] = one_start
                                elif one_start < edge_limits[(side, edge_one)]:
                                    edge_limits[(side, edge_one)] = one_start
                                if (side, edge_two) not in edge_limits:
                                    edge_limits[(side, edge_two)] = two_start
                                elif two_start < edge_limits[(side, edge_two)]:
                                    edge_limits[(side, edge_two)] = two_start
                            elif side == "out":
                                if (side, edge_one) not in edge_limits:
                                    edge_limits[(side, edge_one)] = one_end
                                elif one_end > edge_limits[(side, edge_one)]:
                                    edge_limits[(side, edge_one)] = one_end
                                if (side, edge_two) not in edge_limits:
                                    edge_limits[(side, edge_two)] = two_end
                                elif two_end > edge_limits[(side, edge_two)]:
                                    edge_limits[(side, edge_two)] = two_end
            #For each edge_pair, find starting and ending indices of
            #in, out, and template sequences to construct sequences
            summ_resolution = []
            resolved_sequences = []
            for i, edge_pair in enumerate(sorted(bridged_edges)):
                f.write("Repeat Copy {0}\n".format(i))
                f.write("{0[0]} {0[1]:2}  {1:3} {2[0]} {2[1]}\n".format(
                                             edge_pair[0],
                                             "->",
                                             edge_pair[1]))
                in_start = None
                out_end = None
                out_align = None
                in_align = None
                for side, edge_id in edge_pair:
                    if side == "in" and (side, edge_id) in edge_limits:
                        in_start = edge_limits[(side, edge_id)]
                    elif side == "out" and (side, edge_id) in edge_limits:
                        out_end = edge_limits[(side, edge_id)]
                    if os.path.isfile(cons_align_path.format(side_it[side],
                                                             side,
                                                             edge_id)):
                        cons_align = _read_alignment(
                                cons_align_path.format(side_it[side],
                                                       side,
                                                       edge_id),
                                template,
                                CONS_ALN_RATE)
                        if cons_align and cons_align[0]:
                            #collapse multiple consensus alignments
                            coll_cons_align = _collapse_cons_aln(cons_align)
                            if side == "in":
                                in_align = coll_cons_align
                            elif side == "out":
                                out_align = coll_cons_align
                if not in_align:
                    in_start = 0
                    in_end = 0
                    temp_start = 0
                    #if in_start is None:
                    #    in_start = 0
                else:
                    in_start = in_align.qry_start
                    in_end = in_align.qry_end
                    temp_start = in_align.trg_end
                    #if in_start is None:
                    #    in_start = in_align.qry_start
                    #f.write("CHECK: in qry {0} - {1} of {2}\n".format(in_align.qry_start, 
                    #                in_align.qry_end, in_align.qry_len))
                    #f.write("CHECK: in trg {0} - {1} of {2}\n".format(in_align.trg_start, 
                    #                in_align.trg_end, in_align.trg_len))
                if not out_align:
                    temp_end = 0
                    out_start = 0
                    out_end = 0
                    #if out_end is None:
                    #    out_end = 0
                    out_qry_seq = ""
                    out_trg_seq = ""
                    out_trg_end = 0
                    out_qry_end = 0
                else:
                    temp_end = out_align.trg_start
                    out_start = out_align.qry_start
                    out_end = out_align.qry_end
                    #if out_end is None:
                    #    out_end = out_align.qry_end
                    out_qry_seq = out_align.qry_seq
                    out_trg_seq = out_align.trg_seq
                    out_trg_end = out_align.trg_end
                    out_qry_end = out_align.qry_end
                    #f.write("CHECK: out qry {0} - {1} of {2}\n".format(out_align.qry_start, 
                    #                out_align.qry_end, out_align.qry_len))
                    #f.write("CHECK: out trg {0} - {1} of {2}\n".format(out_align.trg_start, 
                    #                out_align.trg_end, out_align.trg_len))
                f.write("Alignment Indices:\n")
                f.write("{0:10}\t{1:5} - {2:5}\n".format("in",
                                                         in_start, in_end))
                #f.write("{0:10}\t{1:5} - {2:5}\n".format("Template", 
                                                          #temp_start, 
                                                          #temp_end))
                f.write("{0:10}\t{1:5} - {2:5}\n".format("out",
                                                         out_start, out_end))
                #Report gap/overlap length
                gap_len = temp_end - temp_start
                if gap_len >= 0:
                    f.write("{0}\t{1}\n".format("Gap between edges:", gap_len))
                else:
                    f.write("{0}\t{1}\n\n".format("Overlap between edges:",
                                                  -gap_len))
                    #in sequence used to represent overlapping segment
                    #print check of overlapping segment
                    new_temp_end = temp_start
                    new_out_start = None
                    _, out_aln_qry = _index_mapping(out_qry_seq)
                    out_trg_aln, _ = _index_mapping(out_trg_seq)

                    in_edge = edge_pair[0][1]
                    out_edge = edge_pair[1][1]
                    if temp_start >= out_trg_end:
                        #f.write("CHECK, unhelpful case, temp_start {0}\n".format(temp_start))
                        new_out_start = out_qry_end
                    else:
                        #f.write("CHECK: temp_start {0}, len(out_trg_aln) {1}\n".format(temp_start, len(out_trg_aln)))
                        temp_trg_start = temp_start - temp_end
                        if temp_trg_start < len(out_trg_aln):
                            out_aln_ind = out_trg_aln[temp_trg_start]
                            #f.write("CHECK: out_aln_ind {0}, len(out_aln_qry) {1}\n".format(out_aln_ind, len(out_aln_qry)))
                            if out_aln_ind < len(out_aln_qry):
                                new_out_start = (out_start +
                                                 out_aln_qry[out_aln_ind])
                                #f.write("CHECK: new_out_start {0}\n".format(new_out_start))

                    #_check_overlap(
                    #        consensuses[(side_it["in"], "in", in_edge)], 
                    #        template,
                    #        consensuses[(side_it["out"], "out", out_edge)], 
                    #        -gap_len, in_start, in_end, temp_start, temp_end, 
                    #        out_start, out_end,
                    #        new_out_start, in_align.qry_seq, in_align.trg_seq, 
                    #        out_align.qry_seq, out_align.trg_seq, out_trg_aln, 
                    #        out_aln_trg, out_qry_aln, out_aln_qry, 
                    #        out_align.trg_end, out_align.qry_end, 
                    #        in_align, out_align)

                    temp_end = new_temp_end
                    if new_out_start:
                        out_start = new_out_start
                    f.write("Adjusted Alignment Indices:\n")
                    f.write("{0:10}\t{1:5} - {2:5}\n".format("in",
                                                        in_start, in_end))
                    if temp_start != new_temp_end:
                        f.write("{0:10}\t{1:5} - {2:5}\n".format("Template",
                                                        temp_start,
                                                        new_temp_end))
                    f.write("{0:10}\t{1:5} - {2:5}\n\n\n".format("out",
                                                        new_out_start,
                                                        out_end))

                in_edge = edge_pair[0][1]
                out_edge = edge_pair[1][1]
                #header = "_".join(["Repeat_{0}".format(rep), 
                #                        "bridged_copy_{0}".format(i), 
                #                        "in_{0}_{1}_{2}".format(in_edge, 
                #                                                in_start, 
                #                                                in_end), 
                #                        "template_{0}_{1}".format(temp_start, 
                #                                                  temp_end), 
                #                        "out_{0}_{1}_{2}".format(out_edge, 
                #                                                 out_start, 
                #                                                 out_end)])
                header = "repeat_{0}_path_{1}_{2}".format(rep, in_edge, out_edge)
                copy_seq = ""
                if side_it["in"] > 0 and side_it["out"] > 0:
                    copy_seq = _construct_repeat_copy(
                            consensuses[(side_it["in"], "in", in_edge)],
                            template,
                            consensuses[(side_it["out"], "out", out_edge)],
                            in_start, in_end,
                            temp_start, temp_end,
                            out_start, out_end)
                resolved_repeats[header] = copy_seq
                if copy_seq:
                    seq_dict = {header:copy_seq}
                    fp.write_fasta_dict(seq_dict,
                                        resolved_seq_file.format(rep, i))
                #in_str = "".join(["in", str(in_edge)])
                #out_str = "".join(["out", str(out_edge)])
                #summ_resolution.append("|".join([in_str, out_str]))
                summ_resolution.append("{0},{1}".format(in_edge,out_edge))
                resolved_sequences.append(header)
            #summ_vals.extend(["+".join(summ_resolution)])
            summ_vals.append(":".join(summ_resolution))
            summ_vals.append(":".join(resolved_sequences))
    return bridged, resolved_repeats, summ_vals


def int_stats_postscript(rep, repeat_edges, integrated_stats,
                         resolved_rep_path, res_vs_res):
    CONS_ALN_RATE = trestle_config.vals["cons_aln_rate"]

    divs = []
    with open(integrated_stats, "a") as f:
        res_inds = list(range(len(repeat_edges[rep]["in"])))
        f.write("Resolved Repeat Sequence Alignments\n")
        for res_one, res_two in sorted(combinations(res_inds, 2)):
            qry_start = 0
            qry_end = 0
            qry_len = 0
            trg_start = 0
            trg_end = 0
            trg_len = 0
            qry_seq = ""
            trg_seq = ""
            if os.path.isfile(res_vs_res.format(rep, res_one, res_two) and
                resolved_rep_path.format(rep, res_two)):
                res_align = _read_alignment(res_vs_res.format(rep, res_one,
                                                              res_two),
                                            resolved_rep_path.format(rep,
                                                                     res_two),
                                            CONS_ALN_RATE)
                if res_align and res_align[0]:
                    qry_start = res_align[0][0].qry_start
                    qry_end = res_align[0][0].qry_end
                    qry_len = res_align[0][0].qry_len
                    trg_start = res_align[0][0].trg_start
                    trg_end = res_align[0][0].trg_end
                    trg_len = res_align[0][0].trg_len
                    qry_seq = res_align[0][0].qry_seq
                    trg_seq = res_align[0][0].trg_seq
            f.write("Copy {0}|Copy {1}\n".format(res_one, res_two))
            f.write("{0}{1}{2:16}\t{3:5}-{4:5} of {5:5}\n".format(
                    "Copy ", res_one, ":",
                    qry_start, qry_end, qry_len))
            f.write("{0}{1}{2:16}\t{3:5}-{4:5} of {5:5}\n".format(
                    "Copy ", res_two, ":",
                    trg_start, trg_end, trg_len))
            div_rate = _calculate_divergence(qry_seq, trg_seq)
            divs.append(div_rate)
            f.write("{0:26}\t{1:.4f}\n".format("Divergence Rate:", div_rate))
        f.write("\n")
    return _mean(divs)


def _get_partitioning_info(part_list, edges):
    edge_reads = {edge:0 for edge in edges}
    tied_reads = 0
    unassigned_reads = 0
    for _, status, edge, _, _, _ in part_list:
        if status == "Partitioned" and edge != "NA":
            edge_reads[int(edge)] += 1
        elif status == "Tied":
            tied_reads += 1
        elif status == "None":
            unassigned_reads += 1
        else:
            exception_str = "Unknown status {0} in partitioning file"
            raise Exception(exception_str.format(status))
    return edge_reads, tied_reads, unassigned_reads


def _calculate_divergence(qry_seq, trg_seq):
    if not qry_seq or not trg_seq:
        return 0.0

    curr_del = 0
    curr_ins = 0

    match_count = 0
    mis_count = 0
    del_count = 0
    ins_count = 0

    for q, t in zip(qry_seq, trg_seq):
        if q == t:
            match_count += 1
            if curr_del != 0:
                del_count += 1
                curr_del = 0
            if curr_ins != 0:
                ins_count += 1
                curr_ins = 0
        elif q == "-" and t != "-":
            curr_del += 1
            if curr_ins != 0:
                ins_count += 1
                curr_ins = 0
        elif q != "-" and t == "-":
            curr_ins += 1
            if curr_del != 0:
                del_count += 1
                curr_del = 0
        elif q != t:
            mis_count += 1
            if curr_del != 0:
                del_count += 1
                curr_del = 0
            if curr_ins != 0:
                ins_count += 1
                curr_ins = 0
        else:
            raise Exception("No alignment conditions fit, {0} {1}".format(q, t))
    if curr_del != 0:
        del_count += 1
        curr_del = 0
    if curr_ins != 0:
        ins_count += 1
        curr_ins = 0

    indel_sim_rate = 0.0
    total = match_count + mis_count + del_count + ins_count
    if total != 0:
        indel_sim_rate = match_count / float(total)
    return 1 - indel_sim_rate


def _n50(reads_file):
    reads_dict = fp.read_sequence_dict(reads_file)
    read_lengths = sorted([len(x) for x in reads_dict.values()], reverse=True)
    summed_len = 0
    n50 = 0
    for l in read_lengths:
        summed_len += l
        if summed_len >= sum(read_lengths) // 2:
            n50 = l
            break
    return n50


def _get_median(lst):
    if not lst:
        raise ValueError("_get_median() arg is an empty sequence")
    sorted_list = sorted(lst)
    if len(lst) % 2 == 1:
        return sorted_list[len(lst) // 2]
    else:
        mid1 = sorted_list[(len(lst) // 2) - 1]
        mid2 = sorted_list[(len(lst) // 2)]
        return mid1 + mid2 // 2


def _integrate_confirmed_pos(all_in_pos, all_out_pos):
    in_conf, in_rej, in_pos = all_in_pos
    out_conf, out_rej, _ = all_out_pos

    integrated_confirmed = {"total":[], "sub":[], "ins":[], "del":[]}
    integrated_rejected = {"total":[], "sub":[], "ins":[], "del":[]}

    for pos in sorted(in_pos["total"]):
        for pos_type in in_conf:
            if pos in in_conf[pos_type] or pos in out_conf[pos_type]:
                integrated_confirmed[pos_type].append(pos)
            elif pos in in_rej[pos_type] or pos in out_rej[pos_type]:
                integrated_rejected[pos_type].append(pos)
    return integrated_confirmed, integrated_rejected, in_pos


def _get_combos(in_list, out_list):
    all_combos = []
    for combo in _combo_helper(in_list, out_list):
        all_combos.append(combo)
    return all_combos


def _combo_helper(in_list, out_list):
    if not in_list or not out_list:
        yield []
        return
    else:
        in1 = in_list[0]
        for j in range(len(out_list)):
            combo = (in1, out_list[j])
            for rest in _combo_helper(in_list[1:],
                                      out_list[:j] + out_list[j + 1:]):
                yield [combo] + rest


def _get_aln_end(aln_start, aln_seq):
    return aln_start+len(aln_seq.replace("-",""))


"""
def _check_overlap(in_file, temp_file, out_file, overlap, in_start, in_end,
                   temp_start, temp_end, out_start, out_end, new_out_start,
                   in_qry, in_trg, out_qry, out_trg, out_trg_aln, out_aln_trg,
                   out_qry_aln, out_aln_qry, out_trg_end, out_qry_end,
                   in_align, out_align):
    in_dict = fp.read_sequence_dict(in_file)
    in_seq = in_dict.values()[0]
    temp_dict = fp.read_sequence_dict(temp_file)
    temp_seq = temp_dict.values()[0]
    out_dict = fp.read_sequence_dict(out_file)
    out_seq = out_dict.values()[0]
    for i in range(len(out_qry)/50-1, len(out_qry)/50+1):
        aln_ind_st = i*50
        aln_ind_end = (i+1)*50
        if aln_ind_end > len(out_qry):
            aln_ind_end = len(out_qry)
        print 'ALN inds', aln_ind_st, aln_ind_end
        qry_ind_st = out_aln_qry[aln_ind_st]
        if aln_ind_end < len(out_aln_qry):
            qry_ind_end = out_aln_qry[aln_ind_end]
        else:
            qry_ind_end = out_aln_qry[-1]
        print 'QRY inds', qry_ind_st, qry_ind_end
        trg_ind_st = out_aln_trg[aln_ind_st]
        if aln_ind_end < len(out_aln_trg):
            trg_ind_end = out_aln_trg[aln_ind_end]
        else:
            trg_ind_end = out_aln_trg[-1]
        print 'TRG inds', trg_ind_st, trg_ind_end

        print "TRG ALN", out_trg_aln[trg_ind_st:trg_ind_end]
        print "ALN TRG", out_aln_trg[aln_ind_st:aln_ind_end]
        print "QRY ALN", out_qry_aln[qry_ind_st:qry_ind_end]
        print "ALN QRY", out_aln_qry[aln_ind_st:aln_ind_end]
        print "QRY SEQ", out_qry[aln_ind_st:aln_ind_end]
        print "TRG SEQ", out_trg[aln_ind_st:aln_ind_end]
        print
    print 'In end, in template end',in_end,temp_start
    print 'AR In qry end',in_qry[-10:]
    print 'AR In trg end',in_trg[-10:]
    print 'Out old start, old end, new start, out template start', out_start,
    print out_end, new_out_start, temp_end
    print "Out_trg_end", out_trg_end
    print "Out_qry_end", out_qry_end
    print "In align qry inds", in_align.qry_start, in_align.qry_end,
    print in_align.qry_len
    print "In align trg inds", in_align.trg_start, in_align.trg_end,
    print in_align.trg_len
    print "Out align qry inds", out_align.qry_start, out_align.qry_end,
    print out_align.qry_len
    print "Out align trg inds", out_align.trg_start, out_align.trg_end,
    print out_align.trg_len
    print
    print "Overlap:\t{0}".format(overlap)
    print "In seq(-30 to end):\t{0}".format(in_seq[in_end-30:in_end])
    temp_end_seq = temp_seq[temp_start-30:temp_start]
    print "Template seq(-30 to end):\t{0}".format(temp_end_seq)
    #print "Out seq:\t{0}".format(out_seq[out_start:out_end])
    #print "AR In seq:\t{0}".format(in_seq[in_start-10:in_end+10])
    #print "AR Template seq:\t{0}".format(temp_seq[temp_end:temp_start+10])
    #print "AR Out seq:\t{0}".format(out_seq[out_start:out_end+10])
    pre_new_out = out_seq[new_out_start-30:new_out_start]
    post_new_out = out_seq[new_out_start:new_out_start+30]
    print "New out seq(-30 to new start):\t{0}".format(pre_new_out)
    print "New out seq(new_start to +30):\t{0}".format(post_new_out)
    print
"""


def _construct_repeat_copy(in_file, temp_file, out_file, in_start, in_end,
                           temp_start, temp_end, out_start, out_end):
    if (not os.path.isfile(in_file) or
        not os.path.isfile(temp_file) or
        not os.path.isfile(out_file)):
        return ""
    in_dict = fp.read_sequence_dict(in_file)
    temp_dict = fp.read_sequence_dict(temp_file)
    out_dict = fp.read_sequence_dict(out_file)
    if not in_dict or not temp_dict or not out_dict:
        return ""

    in_seq = list(in_dict.values())[0]
    temp_seq = list(temp_dict.values())[0]
    out_seq = list(out_dict.values())[0]
    seq = ''.join([in_seq[in_start:in_end],
                   temp_seq[temp_start:temp_end],
                   out_seq[out_start:out_end]])
    return seq


def init_summary(summary_file):
    with open(summary_file, "w") as f:
        summ_header_labels = ["Repeat_ID", "Path", "Template", "Cov",
                              "#Conf_Pos", "Max_Pos_Gap", "Bridged?",
                              "Support", "Against", "Avg_Div", "Resolution",
                              "Sequences"]
        #spaced_header = map("{:13}".format, summ_header_labels)
        f.write(" ".join(["{:<13}".format(str(x)) for x in summ_header_labels]) + "\n")


def update_summary(summ_items, summary_file):
    (rep_id, graph_path, template_len, avg_cov, summ_vals,
     avg_div, both_resolved_present) = summ_items

    (confirmed_pos, max_pos_gap, bridged,
     support, against, resolution, sequences) = tuple(summ_vals)

    avg_cov = "{:.4f}".format(avg_cov)
    avg_div = "{:.4f}".format(avg_div)
    graph_path = ",".join([str(p) for p in graph_path])
    bridged = bridged and both_resolved_present

    summ_out = [rep_id, graph_path, template_len, avg_cov, confirmed_pos,
                max_pos_gap, bridged, support, against, avg_div, resolution,
                sequences]
    with open(summary_file, "a") as f:
        f.write(" ".join(["{:<13}".format(str(x)) for x in summ_out]) + "\n")

def remove_unneeded_files(repeat_edges, rep, side_labels, side_it, orient_dir,
                          template, extended, pol_temp_dir, pol_ext_dir,
                          pre_edge_reads, pre_partitioning, pre_read_align,
                          partitioning, cons_align, cut_cons_align,
                          read_align, confirmed_pos_path, edge_reads,
                          cut_cons, polishing_dir, cons_vs_cons,
                          int_confirmed_path, repeat_reads, frequency_path,
                          alignment_file, num_pol_iters, iter_pairs):
    add_dir_name = "additional_output"
    add_dir = os.path.join(orient_dir, add_dir_name)
    if not os.path.isdir(add_dir):
        os.mkdir(add_dir)
    pol_name = "polished_{0}.fasta".format(num_pol_iters)
    pol_template = "polished_template.fasta"
    pol_ext = "polished_extended.{0}.{1}.fasta"
    pol_temp_file = os.path.join(pol_temp_dir, pol_name)
    if os.path.exists(pol_temp_file):
        os.rename(pol_temp_file, os.path.join(add_dir, pol_template))
    for side in side_labels:
        for edge_id in repeat_edges[rep][side]:
            pol_ext_file = os.path.join(pol_ext_dir.format(side, edge_id),
                                        pol_name)
            if os.path.exists(pol_ext_file):
                os.rename(pol_ext_file,
                          os.path.join(add_dir,
                                       pol_ext.format(side, edge_id)))

    files_to_remove = [template]
    dirs_to_remove = [pol_temp_dir]
    files_to_move = [repeat_reads, frequency_path, alignment_file]
    if os.path.exists(pol_temp_dir):
        for fil in os.listdir(pol_temp_dir):
            files_to_remove.append(os.path.join(pol_temp_dir, fil))

    for side in side_labels:
        for edge_id in repeat_edges[rep][side]:
            files_to_remove.append(extended.format(side, edge_id))
            curr_pol_ext_dir = pol_ext_dir.format(side, edge_id)
            dirs_to_remove.append(curr_pol_ext_dir)
            if os.path.exists(curr_pol_ext_dir):
                for fil in os.listdir(curr_pol_ext_dir):
                    files_to_remove.append(os.path.join(curr_pol_ext_dir, fil))
            files_to_remove.append(pre_edge_reads.format(side, edge_id))
            files_to_remove.append(pre_read_align.format(side, edge_id))
            for it in range(1, side_it[side] + 1):
                files_to_remove.append(cons_align.format(it, side, edge_id))
                files_to_remove.append(read_align.format(it, side, edge_id))
                files_to_remove.append(edge_reads.format(it, side, edge_id))
                pol_cons = polishing_dir.format(it, side, edge_id)
                dirs_to_remove.append(pol_cons)
                if os.path.exists(pol_cons):
                    for fil in os.listdir(pol_cons):
                        files_to_remove.append(os.path.join(pol_cons, fil))
            for it in range(1, side_it[side]):
                files_to_remove.append(cut_cons_align.format(it, side, edge_id))
                files_to_remove.append(cut_cons.format(it, side, edge_id))
            it = side_it[side]
            files_to_move.append(cut_cons_align.format(it, side, edge_id))
            files_to_move.append(cut_cons.format(it, side, edge_id))

        edge_pairs = sorted(combinations(repeat_edges[rep][side], 2))
        for edge_one, edge_two in edge_pairs:
            for it in range(1, side_it[side]):
                cons_cons_file = cons_vs_cons.format(it, side, edge_one,
                                                     it, side, edge_two)
                files_to_remove.append(cons_cons_file)
            it = side_it[side]
            cons_cons_file = cons_vs_cons.format(it, side, edge_one,
                                                 it, side, edge_two)
            files_to_move.append(cons_cons_file)
        files_to_remove.append(pre_partitioning.format(side))
        for it in range(1, side_it[side]):
            files_to_remove.append(partitioning.format(it, side))
            files_to_remove.append(confirmed_pos_path.format(it, side))
        for it in [0, side_it[side]]:
            files_to_move.append(partitioning.format(it, side))
        it = side_it[side]
        files_to_move.append(confirmed_pos_path.format(it, side))

    last_conf_pos = int_confirmed_path.format(side_it[side_labels[0]],
                                             side_it[side_labels[1]])
    for it1, it2 in iter_pairs:
        curr_conf_pos = int_confirmed_path.format(it1, it2)
        if curr_conf_pos != last_conf_pos:
            files_to_remove.append(curr_conf_pos)
        else:
            files_to_move.append(curr_conf_pos)

    for f in files_to_remove:
        if os.path.exists(f):
            os.remove(f)
    for d in dirs_to_remove:
        if os.path.exists(d):
            os.rmdir(d)
    for f in files_to_move:
        if os.path.exists(f):
            split_path = os.path.split(f)
            new_file = os.path.join(split_path[0], add_dir_name, split_path[1])
            os.rename(f, new_file)


def _mean(lst):
    if not lst:
        return 0
    return sum(lst) / len(lst)
