#(c) 2016-2018 by Authors
#This file is a part of Flye program.
#Released under the BSD license (see LICENSE file)

"""
Modifies repeat graph using the Tresle output
"""

import logging
from itertools import chain
from collections import defaultdict

import flye.utils.fasta_parser as fp
from flye.repeat_graph.graph_alignment import iter_alignments
from flye.six import iteritems

logger = logging.getLogger()


class Connection(object):
    __slots__ = ("id", "path", "sequence")
    def __init__(self, id=None, path=None, sequence=""):
        self.id = id
        self.path = path
        self.sequence = sequence


class RepeatInfo(object):
    __slots__ = ("id", "repeat_path", "all_reads", "in_reads", "out_reads",
                 "sequences", "multiplicity")

    def __init__(self, id, repeat_path, all_reads, in_reads, out_reads,
                 sequences, multiplicity):
        self.id = id
        self.repeat_path = repeat_path
        self.all_reads = all_reads
        self.in_reads = in_reads
        self.out_reads = out_reads
        self.sequences = sequences
        self.multiplicity = multiplicity


def get_simple_repeats(repeat_graph, alignments_file, edge_seqs):
    next_path_id = 1
    path_ids = {}
    repeats_dict = {}
    MULT = 2

    paths_to_resolve = []
    interesting_edges = set()
    for path in repeat_graph.get_unbranching_paths():
        if not path[0].repetitive or path[0].self_complement:
            continue

        is_simple = True
        inputs = set()
        for in_edge in path[0].node_left.in_edges:
            inputs.add(in_edge.edge_id)
            if in_edge.repetitive:
                is_simple = False

        outputs = set()
        for out_edge in path[-1].node_right.out_edges:
            outputs.add(out_edge.edge_id)
            if out_edge.repetitive:
                is_simple = False

        if not is_simple or len(inputs) != MULT or len(outputs) != MULT:
            continue

        paths_to_resolve.append((path, inputs, outputs))
        interesting_edges.update(set([e.edge_id for e in path]))

    interesting_alignments = []
    for read_aln in iter_alignments(alignments_file):
        repeat_read = False
        for edge_aln in read_aln:
            if edge_aln.edge_id in interesting_edges:
                repeat_read = True
        if repeat_read:
            interesting_alignments.append(read_aln)

    for path, inputs, outputs in paths_to_resolve:
        if path[0].edge_id not in path_ids:
            path_ids[path[0].edge_id] = next_path_id
            path_ids[-path[-1].edge_id] = -next_path_id
            next_path_id += 1
        path_id = path_ids[path[0].edge_id]

        repeat_edge_ids = set([e.edge_id for e in path])
        inner_reads = []
        input_reads = defaultdict(list)
        output_reads = defaultdict(list)
        for read_aln in interesting_alignments:
            repeat_read = False
            for edge_aln in read_aln:
                if edge_aln.edge_id in repeat_edge_ids:
                    repeat_read = True
            if not repeat_read:
                continue

            inner_reads.append(read_aln[0].overlap.cur_id)
            for prev_edge, next_edge in zip(read_aln[:-1], read_aln[1:]):
                if (prev_edge.edge_id in inputs and
                        next_edge.edge_id == path[0].edge_id):
                    input_reads[prev_edge.edge_id].append(prev_edge.overlap.cur_id)

                if (prev_edge.edge_id == path[-1].edge_id and
                        next_edge.edge_id in outputs):
                    output_reads[next_edge.edge_id].append(next_edge.overlap.cur_id)

        if (not len(inner_reads) or len(input_reads) != MULT or
                len(output_reads) != MULT):
            continue

        #add edges sequences:
        sequences = {}
        for edge in chain(input_reads, output_reads):
            seq_id = repeat_graph.edges[edge].edge_sequences[0].edge_seq_name
            seq = edge_seqs[seq_id[1:]]
            if seq_id[0] == "-":
                seq = fp.reverse_complement(seq)
            sequences[edge] = seq

        template_seq = ""
        for edge in path:
            seq_id = edge.edge_sequences[0].edge_seq_name
            seq = edge_seqs[seq_id[1:]]
            if seq_id[0] == "-":
                seq = fp.reverse_complement(seq)
            template_seq += seq
        sequences["template"] = template_seq

        #print path_id
        #for h, s in sequences.items():
        #    print h, s[:100]

        repeats_dict[path_id] = RepeatInfo(path_id, [e.edge_id for e in path],
                                           inner_reads, input_reads, output_reads,
                                           sequences, MULT)

    return repeats_dict


def dump_repeats(repeats_info, filename):
    with open(filename, "w") as f:
        for repeat_id, info in iteritems(repeats_info):
            f.write("#Repeat {0}\t{1}\n\n".format(repeat_id, info.multiplicity))

            f.write("#All reads\t{0}\n".format(len(info.all_reads)))
            for read in info.all_reads:
                f.write(read + "\n")
            f.write("\n")

            for in_edge in info.in_reads:
                f.write("#Input {0}\t{1}\n".format(in_edge, len(info.in_reads[in_edge])))
                for read in info.in_reads[in_edge]:
                    f.write(read + "\n")
                f.write("\n")

            for out_edge in info.out_reads:
                f.write("#Output {0}\t{1}\n".format(out_edge, len(info.out_reads[out_edge])))
                for read in info.out_reads[out_edge]:
                    f.write(read + "\n")
                f.write("\n")


def apply_changes(repeat_graph, trestle_results,
                  resolved_repeats_fasta):
    #repeat_graph.output_dot("before.dot")
    connections = _get_connections(trestle_results)
    edges_to_remove = set()
    for conn in connections:

        #FIXME: sometimes trestle outputs empty resolved seque. need a proper fix later
        if not len(resolved_repeats_fasta[conn.sequence]):
            continue

        repeat_graph.separate_path(conn.path, conn.id,
                                   resolved_repeats_fasta[conn.sequence])
        edges_to_remove.update(conn.path[1:-1])

    for edge_id in edges_to_remove:
        edge = repeat_graph.edges[edge_id]
        if not edge.self_complement:
            repeat_graph.remove_edge(repeat_graph.complement_edge(edge))
        repeat_graph.remove_edge(edge)
    #repeat_graph.output_dot("after.dot")


def _get_connections(trestle_results):
    connections = []
    resolved_repeats = set()
    with open(trestle_results, "r") as f:
        for line in f:
            if line.startswith("Repeat"): continue

            tokens = line.strip().split()
            repeat_id, bridged = int(tokens[0]), tokens[6]
            if bridged == "True" and abs(repeat_id) not in resolved_repeats:
                resolved_repeats.add(abs(repeat_id))

                repeat_path = [int(x) for x in tokens[1].split(",")]
                res_1, res_2 = tokens[10].split(":")
                in_1, out_1 = res_1.split(",")
                in_2, out_2 = res_2.split(",")

                seq_1, seq_2 = tokens[11].split(":")
                connection_1 = [int(in_1)] + repeat_path + [int(out_1)]
                connection_2 = [int(in_2)] + repeat_path + [int(out_2)]

                logger.debug("Repeat %s: %s, %s", repeat_id,
                             connection_1, connection_2)

                new_seq_id = "trestle_resolved_" + str(repeat_id) + "_copy_"
                connections.extend([Connection(new_seq_id + "1", connection_1, seq_1),
                                    Connection(new_seq_id + "2", connection_2, seq_2)])

    return connections
