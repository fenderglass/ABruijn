#(c) 2016-2018 by Authors
#This file is a part of Flye program.
#Released under the BSD license (see LICENSE file)

import os
import logging

import flye.utils.fasta_parser as fp
import flye.short_plasmids.unmapped_reads as unmapped
import flye.short_plasmids.circular_sequences as circular
from flye.polishing.alignment import make_alignment
import flye.polishing.polish as pol
from flye.repeat_graph.repeat_graph import EdgeSequence, RgEdge


logger = logging.getLogger()


def assemble_short_plasmids(args, work_dir, contigs_path):
    logger.debug("Extracting unmapped reads")
    reads2contigs_mapping = os.path.join(work_dir, "reads2contigs.paf")
    make_alignment(contigs_path, args.reads, args.threads,
                   work_dir, args.platform, reads2contigs_mapping,
                   reference_mode=True, sam_output=False)

    unmapped_reads_path = os.path.join(work_dir, "unmapped_reads.fasta")
    unmapped.extract_unmapped_reads(args, reads2contigs_mapping,
                                    unmapped_reads_path,
                                    mapping_rate_threshold=0.5)

    logger.debug("Finding self-mappings for unmapped reads")
    unmapped_reads_mapping = os.path.join(work_dir, "unmapped_ava.paf")
    make_alignment(unmapped_reads_path, [unmapped_reads_path], args.threads,
                   work_dir, args.platform, unmapped_reads_mapping,
                   reference_mode=False, sam_output=False)

    logger.debug("Extracting circular reads")
    circular_reads = circular.extract_circular_reads(unmapped_reads_mapping)
    logger.debug("Extracted %d circular reads", len(circular_reads))

    logger.debug("Extracing circular pairs")
    circular_pairs = circular.extract_circular_pairs(unmapped_reads_mapping)
    logger.debug("Extracted %d circular pairs", len(circular_pairs))

    #extracting only the necesssary subset of reads (the entire file could be pretty big)
    interesting_reads = {}
    for read in circular_reads:
        interesting_reads[read] = None
    for pair in circular_pairs:
        interesting_reads[pair[0].query] = None
        interesting_reads[pair[0].target] = None
    for hdr, seq in fp.stream_sequence(unmapped_reads_path):
        if hdr in interesting_reads:
            interesting_reads[hdr] = seq

    trimmed_circular_reads = \
        circular.trim_circular_reads(circular_reads, interesting_reads)
    trimmed_circular_pairs = \
        circular.trim_circular_pairs(circular_pairs, interesting_reads)
    trimmed_sequences_path = os.path.join(work_dir, "trimmed_sequences.fasta")
    fp.write_fasta_dict(dict(list(trimmed_circular_reads.items()) +
                             list(trimmed_circular_pairs.items())),
                        trimmed_sequences_path)

    logger.debug("Clustering circular sequences")
    trimmed_sequences_mapping = os.path.join(work_dir, "trimmed.paf")
    make_alignment(trimmed_sequences_path, [trimmed_sequences_path], args.threads,
                   work_dir, args.platform, trimmed_sequences_mapping,
                   reference_mode=False, sam_output=False)

    plasmids = \
        circular.extract_unique_plasmids(trimmed_sequences_mapping,
                                         trimmed_sequences_path)

    plasmids_raw = os.path.join(work_dir, "plasmids_raw.fasta")
    fp.write_fasta_dict(plasmids, plasmids_raw)
    _, polished_stats = \
        pol.polish(plasmids_raw, [unmapped_reads_path], work_dir, 1,
                   args.threads, args.platform, args.read_type, output_progress=False)

    #extract coverage
    plasmids_with_coverage = {}
    if os.path.isfile(polished_stats):
        with open(polished_stats, "r") as f:
            for line in f:
                if line.startswith("#"): continue
                tokens = line.strip().split()
                seq_id, coverage = tokens[0], int(tokens[2])
                if coverage > 0:
                    plasmids_with_coverage[seq_id] = plasmids[seq_id], coverage

    logger.info("Added %d extra contigs", len(plasmids_with_coverage))

    # remove all unnecesarry files
    os.remove(reads2contigs_mapping)
    os.remove(unmapped_reads_path)
    os.remove(unmapped_reads_mapping)
    os.remove(trimmed_sequences_path)
    os.remove(trimmed_sequences_mapping)

    return plasmids_with_coverage


def update_graph(repeat_graph, plasmids_dict):
    for num, (plasmid, coverage) in enumerate(plasmids_dict.values()):
        new_seq_name = "plasmid_{0}".format(num)
        repeat_graph.edges_fasta[new_seq_name] = plasmid
        new_edge_seq = EdgeSequence("+" + new_seq_name, len(plasmid))
        compl_edge_seq = EdgeSequence("-" + new_seq_name, len(plasmid))

        new_edge_id = max(repeat_graph.edges.keys()) + 1
        node_fwd = repeat_graph.add_node()
        edge_fwd = RgEdge(node_fwd, node_fwd, new_edge_id)
        edge_fwd.edge_sequences.append(new_edge_seq)
        edge_fwd.mean_coverage = coverage
        repeat_graph.add_edge(edge_fwd)

        node_rev = repeat_graph.add_node()
        edge_rev = RgEdge(node_rev, node_rev, -new_edge_id)
        edge_rev.edge_sequences.append(compl_edge_seq)
        edge_rev.mean_coverage = coverage
        repeat_graph.add_edge(edge_rev)

        logger.debug("Added edge %d, length %d, coverage %d",
                     new_edge_id, len(plasmid), coverage)
