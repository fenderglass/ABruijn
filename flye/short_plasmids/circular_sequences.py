#(c) 2016-2018 by Authors
#This file is a part of Flye program.
#Released under the BSD license (see LICENSE file)


import flye.short_plasmids.utils as utils
import flye.short_plasmids.unmapped_reads as unmapped
import flye.utils.fasta_parser as fp
from flye.utils.sam_parser import read_paf, read_paf_grouped
import logging
from flye.six import iteritems

logger = logging.getLogger()

def is_circular_read(hit, max_overhang=150):
    if hit.query != hit.target:
        return False

    if not hit.query_start < hit.query_end < hit.target_start < hit.target_end:
        return False

    if not hit.query_left_overhang() < max_overhang:
        return False

    if not hit.target_right_overhang() < max_overhang:
        return False

    return True


def extract_circular_reads(unmapped_reads_mapping, max_overhang=150):
    circular_reads = {}

    for current_hit in read_paf(unmapped_reads_mapping):
        if is_circular_read(current_hit, max_overhang):
            hit = circular_reads.get(current_hit.query)
            if hit is None or current_hit.query_mapping_length() > \
               hit.query_mapping_length():
                circular_reads[current_hit.query] = current_hit
                #logger.debug("\t" + current_hit.query)

    return circular_reads


def trim_circular_reads(circular_reads, unmapped_reads):
    trimmed_circular_reads = dict()

    for i, (read, hit) in enumerate(iteritems(circular_reads)):
        sequence = unmapped_reads[read][:hit.target_start].upper()
        trimmed_circular_reads["circular_read_" + str(i)] = sequence

    return trimmed_circular_reads


def trim_circular_pairs(circular_pairs, unmapped_reads):
    trimmed_circular_pairs = dict()

    for i, pair in enumerate(circular_pairs):
        lhs = unmapped_reads[pair[0].query]
        rhs = unmapped_reads[pair[0].target]
        trimmed_seq = lhs[pair[1].query_end:pair[0].query_end]
        trimmed_seq += rhs[pair[0].target_end:]
        trimmed_circular_pairs["circular_pair_" + str(i)] = trimmed_seq.upper()

    return trimmed_circular_pairs


def mapping_segments_without_intersection(circular_pair):
    if not circular_pair[1].query_start < circular_pair[1].query_end < \
           circular_pair[0].query_start < circular_pair[0].query_end:
        return False

    if not circular_pair[0].target_start < circular_pair[0].target_end < \
           circular_pair[1].target_start < circular_pair[1].target_end:
        return False

    return True


def extract_circular_pairs(unmapped_reads_mapping, max_overhang=300):
    circular_pairs = []
    used_reads = set()

    #each hit group stores alginmemnts for each (query, target) pair
    for hit_group in read_paf_grouped(unmapped_reads_mapping):
        if hit_group[0].query == hit_group[0].target:
            continue

        if hit_group[0].query in used_reads or hit_group[0].target in used_reads:
            continue

        circular_pair = [None, None]
        has_overlap = False
        is_circular = False

        for hit in hit_group:
            if (not has_overlap and hit.query_right_overhang() < max_overhang and
                    hit.target_left_overhang() < max_overhang):
                has_overlap = True
                circular_pair[0] = hit
                continue

            if (not is_circular and hit.query_left_overhang() < max_overhang and
                   hit.target_right_overhang() < max_overhang):
                is_circular = True
                circular_pair[1] = hit

        if (has_overlap and is_circular and
                mapping_segments_without_intersection(circular_pair)):
            circular_pairs.append(circular_pair)
            used_reads.add(circular_pair[0].target)
            used_reads.add(circular_pair[0].query)

            #logger.debug("\t" + circular_pair[0].target + "_" + circular_pair[0].query)

    return circular_pairs


def extract_unique_plasmids(trimmed_reads_mapping, trimmed_reads_path,
                            mapping_rate_threshold=0.8,
                            max_length_difference=500,
                            min_sequence_length=1000):
    trimmed_reads = set()
    for hit in read_paf(trimmed_reads_mapping):
        trimmed_reads.add(hit.query)
        trimmed_reads.add(hit.target)

    trimmed_reads = list(trimmed_reads)
    n_trimmed_reads = len(trimmed_reads)
    read2int = dict()
    int2read = dict()

    for i in range(n_trimmed_reads):
        read2int[trimmed_reads[i]] = i
        int2read[i] = trimmed_reads[i]

    similarity_graph = [[] for _ in range(n_trimmed_reads)]

    #each hit group stores alginmemnts for each (query, target) pair
    for hit_group in read_paf_grouped(trimmed_reads_mapping):
        if hit_group[0].query == hit_group[0].target:
            continue

        query_mapping_segments = []
        target_mapping_segments = []
        for hit in hit_group:
            query_mapping_segments.append(unmapped.MappingSegment(hit.query_start,
                                                                  hit.query_end))
            target_mapping_segments.append(unmapped.MappingSegment(hit.target_start,
                                                                   hit.target_end))

        query_length = hit_group[0].query_length
        target_length = hit_group[0].target_length
        query_mapping_rate = unmapped.calc_mapping_rate(query_length,
                                                        query_mapping_segments)
        target_mapping_rate = unmapped.calc_mapping_rate(target_length,
                                                         target_mapping_segments)

        if (query_mapping_rate > mapping_rate_threshold or
                target_mapping_rate > mapping_rate_threshold):
            #abs(query_length - target_length) < max_length_difference:
            vertex1 = read2int[hit_group[0].query]
            vertex2 = read2int[hit_group[0].target]
            similarity_graph[vertex1].append(vertex2)
            similarity_graph[vertex2].append(vertex1)

    connected_components, n_components = \
        utils.find_connected_components(similarity_graph)

    groups = [[] for _ in range(n_components)]
    for i in range(len(connected_components)):
        groups[connected_components[i]].append(int2read[i])

    #for g in groups:
    #    logger.debug("Group {0}".format(len(g)))
    #    for s in g:
    #        logger.debug("\t{0}".format(seq_lengths[s]))


    groups = [group for group in groups if len(group) > 1]
    trimmed_reads_dict = fp.read_sequence_dict(trimmed_reads_path)
    unique_plasmids = dict()

    for group in groups:
        sequence = trimmed_reads_dict[group[0]]
        if len(sequence) >= min_sequence_length:
            unique_plasmids[group[0]] = sequence

    return unique_plasmids
