#(c) 2016 by Authors
#This file is a part of ABruijn program.
#Released under the BSD license (see LICENSE file)

"""
This module provides parsing/dumping of reads alignment
to the repreat graph (as used internally in Flye)
"""


class OverlapRange(object):
    __slots__ = ("cur_id", "cur_len", "cur_start", "cur_end",
                 "ext_id", "ext_len", "ext_start", "ext_end",
                 "left_shift", "right_shift", "score", "divergence")

    def __init__(self, cur_id, cur_len, cur_start, cur_end,
                 ext_id, ext_len, ext_start, ext_end,
                 left_shift, right_shift, score, divergence):
        self.cur_id = cur_id
        self.cur_len = cur_len
        self.cur_start = cur_start
        self.cur_end = cur_end
        self.ext_id = ext_id
        self.ext_len = ext_len
        self.ext_start = ext_start
        self.ext_end = ext_end
        self.left_shift = left_shift
        self.right_shift = right_shift
        self.score = score
        self.divergence = divergence


class GraphAlignment(object):
    __slots__ = ("edge_id", "overlap")

    def __init__(self, edge_id, overlap):
        self.edge_id = edge_id
        self.overlap = overlap


def iter_alignments(filename):
    """
    Returns alignment generator
    """
    #alignments = []
    current_chain = []
    with open(filename, "r") as f:
        for line in f:
            if not line: continue

            tokens = line.strip().split()
            if tokens[0] == "Chain":
                if current_chain:
                    yield current_chain
                    #alignments.append(current_chain)
                    current_chain = []

            elif tokens[0] == "Aln":
                (edge_id, cur_id, cur_start, cur_end, cur_len,
                ext_id, ext_start, ext_end, ext_len, left_shift,
                right_shift, score, divergence) = tokens[1:]

                ovlp = OverlapRange(cur_id, int(cur_len), int(cur_start), int(cur_end),
                                    ext_id, int(ext_len), int(ext_start), int(ext_end),
                                    int(left_shift), int(right_shift), int(score),
                                    float(divergence))
                current_chain.append(GraphAlignment(_to_signed_id(int(edge_id)), ovlp))

            else:
                raise Exception("Error parsing " + filename)

        if current_chain:
            yield current_chain


#TODO:
#def write_alignments(alignments, filename):
#    pass


def _to_signed_id(unsigned_id):
    return -(unsigned_id + 1) // 2 if unsigned_id % 2 else unsigned_id // 2 + 1


def _to_unsigned_id(signed_id):
    unsigned_id = abs(signed_id) * 2 - 2
    return unsigned_id + int(signed_id < 0)
