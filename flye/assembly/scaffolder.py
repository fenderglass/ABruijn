#(c) 2017 by Authors
#This file is a part of Flye program.
#Released under the BSD license (see LICENSE file)

"""
Final output generator
"""

import logging

import flye.utils.fasta_parser as fp
import flye.config.py_cfg as cfg
from flye.six import iteritems

logger = logging.getLogger()


def generate_scaffolds(contigs_file, links_file, out_scaffolds):

    contigs_fasta = fp.read_sequence_dict(contigs_file)
    scaffolds_fasta = {}
    used_contigs = set()

    connections = {}
    if links_file is not None:
        with open(links_file, "r") as f:
            for line in f:
                line = line.strip()
                if not line: continue
                ctg_1, sign_1, ctg_2, sign_2 = line.split("\t")
                if ctg_1 in contigs_fasta and ctg_2 in contigs_fasta:
                    connections[sign_1 + ctg_1] = sign_2 + ctg_2
                    connections[rc(sign_2) + ctg_2] = rc(sign_1) + ctg_1

    scaffolds_fasta = {}
    scaffolds_seq = {}
    for ctg in contigs_fasta:
        if ctg in used_contigs: continue

        used_contigs.add(ctg)
        scf = ["-" + ctg]
        #extending right
        while (scf[-1] in connections and
               unsigned(connections[scf[-1]]) not in used_contigs):
            scf.append(connections[scf[-1]])
            used_contigs.add(unsigned(scf[-1]))

        for i, ctg in enumerate(scf):
            scf[i] = rc(ctg[0]) + unsigned(ctg)
        scf = scf[::-1]

        #extending left
        while (scf[-1] in connections and
               unsigned(connections[scf[-1]]) not in used_contigs):
            scf.append(connections[scf[-1]])
            used_contigs.add(unsigned(scf[-1]))

        #generating sequence interleaved by Ns
        if len(scf) == 1:
            scaffolds_fasta[unsigned(ctg)] = contigs_fasta[unsigned(ctg)]
            scaffolds_seq[unsigned(ctg)] = scf
        else:
            scf_name = "scaffold_" + unsigned(scf[0]).strip("contig_")
            scaffolds_seq[scf_name] = scf
            scf_seq = []
            for scf_ctg in scf:
                if scf_ctg[0] == "+":
                    scf_seq.append(contigs_fasta[unsigned(scf_ctg)])
                else:
                    scf_seq.append(fp.reverse_complement(
                                    contigs_fasta[unsigned(scf_ctg)]))
            gap = "N" * cfg.vals["scaffold_gap"]
            scaffolds_fasta[scf_name] = gap.join(scf_seq)

    fp.write_fasta_dict(scaffolds_fasta, out_scaffolds)
    return scaffolds_seq


class SeqStats(object):
    __slots__ = ("name", "length", "coverage", "circular",
                 "repeat", "mult", "telomere", "alt_group", "graph_path")

    def __init__(self, name="", length="", coverage="", circular="N",
                 repeat="N", mult="1", telomere="none",
                 alt_group="*", graph_path=""):
        self.name = name
        self.length = length
        self.coverage = coverage
        self.circular = circular
        self.repeat = repeat
        self.mult = mult
        self.telomere = telomere
        self.alt_group = alt_group
        self.graph_path = graph_path

    def print_out(self, handle):
        handle.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\n"
                     .format(self.name, self.length, self.coverage,
                             self.circular, self.repeat, self.mult,
                             self.alt_group, self.graph_path))


def generate_stats(repeat_file, polished_file, scaffolds, out_stats):
    """
    Compiles information from multiple stages
    """
    contigs_stats = {}
    header_line = "#seq_name\tlength\tcov.\tcirc.\trepeat\tmult.\talt_group\tgraph_path"
    for line in open(repeat_file, "r"):
        if line.startswith("#"): continue
        tokens = line.strip().split("\t")
        contigs_stats[tokens[0]] = SeqStats(*tokens)

    if polished_file is not None:
        for line in open(polished_file, "r"):
            if line.startswith("#"): continue
            tokens = line.strip().split("\t")

            #update multiplicity proportionally
            cov_rate = float(tokens[2]) / (float(contigs_stats[tokens[0]].coverage) + 1)
            contigs_stats[tokens[0]].mult = \
                    max(1, int(float(contigs_stats[tokens[0]].mult) * cov_rate))

            #update length and coverage
            contigs_stats[tokens[0]].length = tokens[1]
            contigs_stats[tokens[0]].coverage = tokens[2]

    scaffolds_stats = {}
    for scf, scf_seq in iteritems(scaffolds):
        scaffolds_stats[scf] = SeqStats(scf)
        scf_length = sum([int(contigs_stats[unsigned(c)].length) for c in scf_seq])
        scf_length += (len(scf_seq) - 1) * cfg.vals["scaffold_gap"]
        scaffolds_stats[scf].length = str(scf_length)

        scf_cov = _mean([int(contigs_stats[unsigned(c)].coverage) for c in scf_seq])
        scaffolds_stats[scf].coverage = str(scf_cov)

        scaffolds_stats[scf].repeat = contigs_stats[unsigned(scf_seq[0])].repeat
        scaffolds_stats[scf].circular = contigs_stats[unsigned(scf_seq[0])].circular
        scaffolds_stats[scf].alt_group = contigs_stats[unsigned(scf_seq[0])].alt_group

        scf_mult = min([int(contigs_stats[unsigned(c)].mult) for c in scf_seq])
        scaffolds_stats[scf].mult = str(scf_mult)

        #telomere information
        telomere_left = contigs_stats[unsigned(scf_seq[0])].telomere
        telomere_right = contigs_stats[unsigned(scf_seq[-1])].telomere
        if scf_seq[0][0] == "+":
            scf_left = telomere_left in ["left", "both"]
        else:
            scf_left = telomere_left in ["right", "both"]
        if scf_seq[-1][0] == "+":
            scf_right = telomere_right in ["right", "both"]
        else:
            scf_right = telomere_right in ["left", "both"]
        #if scf_left and scf_right: scaffolds_stats[scf].telomere = "both"
        #elif scf_left and not scf_right: scaffolds_stats[scf].telomere = "left"
        #elif not scf_left and scf_right: scaffolds_stats[scf].telomere = "right"
        #else: scaffolds_stats[scf].telomere = "none"

        #graph path
        path = []
        for ctg in scf_seq:
            ctg_path = contigs_stats[unsigned(ctg)].graph_path
            if ctg[0] == "-":
                ctg_path = ",".join([str(-int(x))
                                     for x in ctg_path.split(",")][::-1])
            path.append(ctg_path)
        prefix = "*," if scf_left else ""
        suffix = ",*" if scf_right else ""
        scaffolds_stats[scf].graph_path = prefix + ",??,".join(path) + suffix

    with open(out_stats, "w") as f:
        f.write(header_line + "\n")
        #for scf in sorted(scaffolds_stats,
        #                  key=lambda x: int(x.rsplit("_", 1)[-1])):
        for scf in sorted(scaffolds_stats,
                          key=lambda x: int(scaffolds_stats[x].length), reverse=True):
            scaffolds_stats[scf].print_out(f)

    total_length = sum([int(x.length) for x in scaffolds_stats.values()])
    if total_length == 0: return

    num_scaffolds = len(scaffolds_stats)
    num_contigs = sum([len(x) for x in scaffolds.values()])

    scaffold_lengths = [int(s.length) for s in scaffolds_stats.values()]
    contig_lengths = []
    for scf in scaffolds.values():
        for ctg in scf:
            contig_lengths.append(int(contigs_stats[unsigned(ctg)].length))
    largest_scf = max(scaffold_lengths)

    #ctg_n50 = _calc_n50(contig_lengths, total_length)
    scf_n50 = _calc_n50(scaffold_lengths, total_length)

    mean_read_cov = 0
    for scf in scaffolds_stats.values():
        mean_read_cov += int(scf.length) * int(scf.coverage)
    mean_read_cov //= total_length

    logger.info("Assembly statistics:\n\n"
                "\tTotal length:\t%d\n"
                "\tFragments:\t%d\n"
                #"\tContigs N50:\t{2}\n"
                "\tFragments N50:\t%d\n"
                "\tLargest frg:\t%d\n"
                "\tScaffolds:\t%d\n"
                "\tMean coverage:\t%d\n",
                total_length, num_scaffolds, scf_n50, largest_scf,
                num_contigs - num_scaffolds, mean_read_cov)


def short_statistics(fasta_file):
    lengths = list(fp.read_sequence_lengths(fasta_file).values())
    if not lengths:
        return 0, 0
    total_size = sum(lengths)
    return total_size, _calc_n50(lengths, total_size)


def rc(sign):
    return "+" if sign == "-" else "-"


def unsigned(ctg):
    return ctg[1:]


def _mean(lst):
    if not lst: return 0
    return sum(lst) // len(lst)


def _calc_n50(scaffolds_lengths, assembly_len):
    n50 = 0
    sum_len = 0
    for l in sorted(scaffolds_lengths, reverse=True):
        sum_len += l
        if sum_len > assembly_len // 2:
            n50 = l
            break
    return n50
