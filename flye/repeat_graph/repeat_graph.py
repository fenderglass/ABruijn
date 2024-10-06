#(c) 2016 by Authors
#This file is a part of ABruijn program.
#Released under the BSD license (see LICENSE file)

"""
This module provides repeat graph parsing/serializing functions,
(as output by repeat graph construction module)
as well as some basic operations
"""

class RgEdge(object):
    __slots__ = ("node_left", "node_right", "edge_id", "repetitive",
                 "self_complement", "resolved", "mean_coverage",
                 "alt_group", "edge_sequences")

    def __init__(self, node_left=None, node_right=None,
                 edge_id=None):
        self.node_left = node_left
        self.node_right = node_right
        self.edge_id = edge_id

        self.repetitive = False
        self.self_complement = False
        self.resolved = False
        self.mean_coverage = 0
        self.edge_sequences = []
        self.alt_group = -1

    def length(self):
        if not self.edge_sequences:
            return 0

        return sum([s.edge_seq_len
                    for s in self.edge_sequences]) // len(self.edge_sequences)

    def __repr__(self):
        return "(id={0}, len={1}, cov={2} rep={3})" \
                .format(self.edge_id, self.length(),
                        self.mean_coverage, self.repetitive)


class EdgeSequence(object):
    __slots__ = ("edge_seq_name", "edge_seq_len", "orig_seq_id", "orig_seq_len",
                 "orig_seq_start", "orig_seq_end")

    def __init__(self, edge_seq_name=None, edge_seq_len=0, orig_seq_id="*",
                 orig_seq_len=0, orig_seq_start=0, orig_seq_end=0):
        self.edge_seq_name = edge_seq_name
        self.edge_seq_len = edge_seq_len

        self.orig_seq_id = orig_seq_id
        self.orig_seq_len = orig_seq_len
        self.orig_seq_start = orig_seq_start
        self.orig_seq_end = orig_seq_end


class RgNode(object):
    __slots__ = ("in_edges", "out_edges")

    def __init__(self):
        self.in_edges = []
        self.out_edges = []

    def is_bifurcation(self):
        return len(self.in_edges) != 1 or len(self.out_edges) != 1


class RepeatGraph(object):
    __slots__ = ("nodes", "edges", "edges_fasta")

    def __init__(self, edges_fasta):
        self.nodes = []
        self.edges = {} #key = edge id
        self.edges_fasta = edges_fasta

    def add_node(self):
        self.nodes.append(RgNode())
        return self.nodes[-1]

    def add_edge(self, edge):
        self.edges[edge.edge_id] = edge
        edge.node_left.out_edges.append(edge)
        edge.node_right.in_edges.append(edge)

    def remove_edge(self, edge):
        _remove_from_list(edge.node_left.out_edges, edge)
        _remove_from_list(edge.node_right.in_edges, edge)
        del self.edges[edge.edge_id]

    def complement_edge(self, edge):
        if edge.self_complement:
            return edge
        return self.edges[-edge.edge_id]

    def get_unbranching_paths(self):
        unbranching_paths = []
        visited_edges = set()

        for edge in self.edges.values():
            if edge in visited_edges:
                continue

            traversed = [edge]
            if not edge.self_complement:
                cur_node = edge.node_left
                while (not cur_node.is_bifurcation() and
                       len(cur_node.in_edges) > 0 and
                       cur_node.in_edges[0] not in visited_edges and
                       not cur_node.in_edges[0].self_complement):
                    traversed.append(cur_node.in_edges[0])
                    visited_edges.add(cur_node.in_edges[0])
                    cur_node = cur_node.in_edges[0].node_left

                traversed = traversed[::-1]
                cur_node = edge.node_right

                while (not cur_node.is_bifurcation() and
                       len(cur_node.out_edges) > 0 and
                       cur_node.out_edges[0] not in visited_edges and
                       not cur_node.out_edges[0].self_complement):
                    traversed.append(cur_node.out_edges[0])
                    visited_edges.add(cur_node.out_edges[0])
                    cur_node = cur_node.out_edges[0].node_right

            unbranching_paths.append(traversed)

        return unbranching_paths

    def load_from_file(self, filename):
        id_to_node = {}
        cur_edge = None
        with open(filename, "r") as f:
            for line in f:
                tokens = line.strip().split()
                if tokens[0] == "Edge":
                    (edge_id, left_node, right_node, repetitive,
                     self_complement, resolved, mean_coverage, alt_group) = tokens[1:]
                    if left_node not in id_to_node:
                        id_to_node[left_node] = self.add_node()
                    if right_node not in id_to_node:
                        id_to_node[right_node] = self.add_node()

                    cur_edge = RgEdge(id_to_node[left_node],
                                      id_to_node[right_node],
                                      _to_signed_id(int(edge_id)))
                    cur_edge.repetitive = bool(int(repetitive))
                    cur_edge.self_complement = bool(int(self_complement))
                    cur_edge.resolved = bool(int(resolved))
                    cur_edge.mean_coverage = int(mean_coverage)
                    cur_edge.alt_group = int(alt_group)
                    self.add_edge(cur_edge)

                elif tokens[0] == "Sequence":
                    (edge_seq_name, edge_seq_len, orig_seq_id,
                    orig_seq_len, orig_seq_start, orig_seq_end) = tokens[1:]
                    edge_seq = EdgeSequence(edge_seq_name, int(edge_seq_len),
                                            orig_seq_id, orig_seq_len,
                                            orig_seq_start, orig_seq_end)
                    cur_edge.edge_sequences.append(edge_seq)

                else:
                    raise Exception("Error parsing " + filename)

    def dump_to_file(self, filename):
        next_node_id = 0
        node_ids = {}
        for node in self.nodes:
            node_ids[node] = next_node_id
            next_node_id += 1

        with open(filename, "w") as f:
            for edge in self.edges.values():
                f.write("Edge\t{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\n"
                    .format(_to_unsigned_id(edge.edge_id), node_ids[edge.node_left],
                            node_ids[edge.node_right], int(edge.repetitive),
                            int(edge.self_complement), int(edge.resolved),
                            int(edge.mean_coverage), int(edge.alt_group)))

                for seq in edge.edge_sequences:
                    f.write("\tSequence\t{0} {1} {2} {3} {4} {5}\n"
                        .format(seq.edge_seq_name, seq.edge_seq_len,
                                seq.orig_seq_id, seq.orig_seq_len,
                                seq.orig_seq_start, seq.orig_seq_end))

    def output_dot(self, filename):
        next_node_id = 0
        node_ids = {}
        for node in self.nodes:
            node_ids[node] = next_node_id
            next_node_id += 1

        with open(filename, "w") as f:
            f.write("digraph {\nnodesep = 0.5;\n"
                    "node [shape = circle, label = \"\", height = 0.3];")
            for edge in self.edges.values():
                f.write("{0} -> {1} [label = \"{2}\", color = \"{3}\"]\n"
                    .format(node_ids[edge.node_left], node_ids[edge.node_right],
                            edge.edge_id, "red" if edge.repetitive else "black"))
            f.write("}")

    def separate_path(self, graph_path, new_seq_id, new_seq_seq):
        """
        Separates the path (and its complement) on the graph.
        First and last edges in the path are disconnected from the graph,
        and then connected by a new edge. For example,
        a path (A -> X -> Y -> ... -> Z -> B) will be transformed
        into a new unbranching path A -> N -> B,
        where N represents a new edge with the given sequence.
        The intermediate path edges remain in the graph (their mean coverage 
        is modified accordingly) and they acquire the attribute 'resolved'.
        Resolved edges could later be cleaned up by using XXX function.
        """

        def separate_one(edges_path, new_edge_id, new_edge_seq):
            left_node = self.add_node()
            _remove_from_list(edges_path[0].node_right.in_edges, edges_path[0])
            edges_path[0].node_right = left_node
            left_node.in_edges.append(edges_path[0])

            path_coverage = (edges_path[0].mean_coverage +
                             edges_path[-1].mean_coverage) // 2
            for mid_edge in edges_path[1:-1]:
                mid_edge.resolved = True
                mid_edge.mean_coverage -= path_coverage

            right_node = left_node
            if len(edges_path) > 2:
                right_node = self.add_node()
                new_edge = RgEdge(left_node, right_node, new_edge_id)
                self.add_edge(new_edge)
                new_edge.mean_coverage = path_coverage
                new_edge.edge_sequences.append(new_edge_seq)

            _remove_from_list(edges_path[-1].node_left.out_edges, edges_path[-1])
            edges_path[-1].node_left = right_node
            right_node.out_edges.append(edges_path[-1])

        if len(graph_path) < 2:
            raise Exception("Path is too short")

        fwd_edges = []
        rev_edges = []
        for e in graph_path:
            if e not in self.edges:
                raise Exception("Nonexistent edge")
            fwd_edges.append(self.edges[e])
            rev_edges.append(self.complement_edge(self.edges[e]))
        rev_edges = rev_edges[::-1]

        new_edge_seq = EdgeSequence("+" + new_seq_id, len(new_seq_seq))
        compl_edge_seq = EdgeSequence("-" + new_seq_id, len(new_seq_seq))
        self.edges_fasta[new_seq_id] = new_seq_seq

        new_edge_id = max(self.edges.keys()) + 1
        separate_one(fwd_edges, new_edge_id, new_edge_seq)
        separate_one(rev_edges, -new_edge_id, compl_edge_seq)


def _remove_from_list(lst, elem):
    lst = [x for x in lst if x != elem]


def _to_signed_id(unsigned_id):
    return -(unsigned_id + 1) // 2 if unsigned_id % 2 else unsigned_id // 2 + 1


def _to_unsigned_id(signed_id):
    unsigned_id = abs(signed_id) * 2 - 2
    return unsigned_id + int(signed_id < 0)
