#(c) 2019 by Authors
#This file is a part of Flye program.
#Released under the BSD license (see LICENSE file)

"""
Provides multithreaded parser for SAM files
"""

import os
import re
import sys
from collections import namedtuple, defaultdict
import subprocess
import logging
import multiprocessing
import ctypes
import time
import random
from copy import copy

#In Python2, everything is bytes (=str)
#In Python3, we are doing IO in bytes, but everywhere else strngs = unicode
if sys.version_info < (3, 0):
    from string import maketrans
    _STR = lambda x: x
    _BYTES = lambda x: x
else:
    maketrans = bytes.maketrans
    _STR = bytes.decode
    _BYTES = str.encode

from flye.six import iteritems

import flye.utils.fasta_parser as fp
from flye.utils.utils import get_median

logger = logging.getLogger()

SAMTOOLS_BIN = "flye-samtools"
Alignment = namedtuple("Alignment", ["qry_id", "trg_id", "qry_start", "qry_end",
                                     "qry_sign", "qry_len", "trg_start",
                                     "trg_end", "trg_sign", "trg_len",
                                     "qry_seq", "trg_seq", "err_rate",
                                     "is_secondary", "is_supplementary", "map_qv"])


ContigRegion = namedtuple("ContigRegion", ["ctg_id", "start", "end"])


class AlignmentException(Exception):
    pass


class PafHit(object):
    """
    Stores paf alignment
    """
    __slots__ = ("query", "query_length", "query_start", "query_end",
                 "target", "target_length", "target_start", "target_end")
    def __init__(self, raw_hit):
        hit = raw_hit.split()

        self.query = hit[0]
        self.query_length = int(hit[1])
        self.query_start = int(hit[2])
        self.query_end = int(hit[3])

        self.target = hit[5]
        self.target_length = int(hit[6])
        self.target_start = int(hit[7])
        self.target_end = int(hit[8])

    def query_mapping_length(self):
        return self.query_end - self.query_start + 1

    def target_mapping_length(self):
        return self.target_end - self.target_start + 1

    def query_left_overhang(self):
        return self.query_start

    def query_right_overhang(self):
        return self.query_length - self.query_end + 1

    def target_left_overhang(self):
        return self.target_start

    def target_right_overhang(self):
        return self.target_length - self.target_end + 1


def read_paf(filename):
    """
    Streams out paf alignments
    """
    with open(filename, "rb") as f:
        for raw_hit in f:
            yield PafHit(_STR(raw_hit))


def read_paf_grouped(filename):
    """
    Outputs chunks of alignments for each (query, target)pair.
    Assumes that PAF alignment is already sorted by query.
    """
    prev_hit = None
    target_hits = defaultdict(list)
    for hit in read_paf(filename):
        if prev_hit is not None and hit.query != prev_hit.query:
            for trg in sorted(target_hits):
                yield target_hits[trg]
            target_hits = defaultdict(list)

        target_hits[hit.target].append(hit)
        prev_hit = hit

    if len(target_hits):
        for trg in sorted(target_hits):
            yield target_hits[trg]


class SynchonizedChunkManager(object):
    """
    Helper class to organize multiprocessing.
    Stores the list of reference segments / chunks and can
    return them in multiple threds
    """
    def __init__(self, reference_fasta, multiproc_manager, chunk_size=None):
        #prepare list of chunks to read
        self.fetch_list = []
        self.chunk_size = chunk_size

        #will be shared between processes
        #self.shared_manager = multiprocessing.Manager()
        self.shared_num_jobs = multiprocessing.Value(ctypes.c_int, 0)
        if multiproc_manager:
            self.shared_lock = multiproc_manager.Lock()
            self.shared_eof = multiprocessing.Value(ctypes.c_bool, False)
        else:
            self.shared_lock = None
            self.shared_eof = ctypes.c_bool(False)


        for ctg_id in reference_fasta:
            ctg_len = len(reference_fasta[ctg_id])
            chunk_size = self.chunk_size if self.chunk_size is not None else ctg_len
            for i in range(0, max(ctg_len // chunk_size, 1)):
                reg_start = i * chunk_size
                reg_end = (i + 1) * chunk_size
                if ctg_len - reg_end < chunk_size:
                    reg_end = ctg_len
                self.fetch_list.append(ContigRegion(ctg_id, reg_start, reg_end))
                #logger.debug("Region: {0} {1} {2}".format(ctg_id, reg_start, reg_end))

        if len(self.fetch_list) == 0:
            self.shared_eof.value = True

    def is_done(self):
        return self.shared_eof.value

    def get_chunk(self):
        job_id = None
        while True:
            if self.shared_lock:
                self.shared_lock.acquire()

            if self.shared_eof.value:
                if self.shared_lock:
                    self.shared_lock.release()
                return None
            
            job_id = self.shared_num_jobs.value
            self.shared_num_jobs.value = self.shared_num_jobs.value + 1
            if self.shared_num_jobs.value == len(self.fetch_list):
                self.shared_eof.value = True

            if self.shared_lock:
                self.shared_lock.release()
            break

            time.sleep(0.01)

        parsed_contig = _BYTES(self.fetch_list[job_id].ctg_id)
        region = self.fetch_list[job_id]
        return region


class SynchronizedSamReader(object):
    """
    Parses SAM file in multiple threads.
    """
    def __init__(self, sam_alignment, reference_fasta, multiproc_manager,
                 max_coverage=None, use_secondary=False):
        #check that alignment exists
        if not os.path.exists(sam_alignment):
            raise AlignmentException("Can't open {0}".format(sam_alignment))
        if not os.path.exists(sam_alignment + ".bai"):
            raise AlignmentException("Bam not indexed: {0}".format(sam_alignment))

        #will not be changed during exceution, each process has its own copy
        self.aln_path = sam_alignment
        self.max_coverage = max_coverage
        self.use_secondary = use_secondary
        self.cigar_parser = re.compile(b"[0-9]+[MIDNSHP=X]")

        #self.shared_manager = multiprocessing.Manager()
        self.ref_fasta = dict() if multiproc_manager == None else multiproc_manager.dict()
        for (h, s) in iteritems(reference_fasta):
            self.ref_fasta[_BYTES(h)] = _BYTES(s)

    def get_region_sequence(self, region_id, region_start=None, region_end=None):
        parsed_contig = _BYTES(region_id)
        contig_str = self.ref_fasta[parsed_contig]
        if region_start is None:
            region_start = 0
        if region_end is None:
            region_end = len(contig_str)

        return _STR(contig_str[region_start:region_end])

    def _parse_cigar(self, cigar_str, read_str, ctg_str, ctg_pos):
        #ctg_str = self.ref_fasta[ctg_name]
        trg_seq = []
        qry_seq = []
        trg_start = ctg_pos - 1
        trg_pos = ctg_pos - 1
        qry_start = 0
        qry_pos = 0

        left_hard = True
        left_soft = True
        hard_clipped_left = 0
        hard_clipped_right = 0
        soft_clipped_left = 0
        soft_clipped_right = 0
        for token in self.cigar_parser.findall(cigar_str):
            size, op = int(token[:-1]), token[-1:]
            if op == b"H":
                if left_hard:
                    qry_start += size
                    hard_clipped_left += size
                else:
                    hard_clipped_right += size
            elif op == b"S":
                qry_pos += size
                if left_soft:
                    soft_clipped_left += size
                else:
                    soft_clipped_right += size
            elif op in b"MX=":
                qry_seq.append(read_str[qry_pos : qry_pos + size].upper())
                trg_seq.append(ctg_str[trg_pos : trg_pos + size].upper())
                qry_pos += size
                trg_pos += size
            elif op == b"I":
                qry_seq.append(read_str[qry_pos : qry_pos + size].upper())
                trg_seq.append(b"-" * size)
                qry_pos += size
            elif op == b"D":
                qry_seq.append(b"-" * size)
                trg_seq.append(ctg_str[trg_pos : trg_pos + size].upper())
                trg_pos += size
            else:
                raise AlignmentException("Unsupported CIGAR operation: " + str(op))
            left_hard = False
            if op != b"H":
                left_soft = False

        trg_seq = b"".join(trg_seq)
        qry_seq = b"".join(qry_seq)
        matches = 0
        for i in range(len(trg_seq)):
            if trg_seq[i] == qry_seq[i]:
                matches += 1
        err_rate = 1 - matches / len(trg_seq)

        trg_end = trg_pos
        qry_end = qry_pos + hard_clipped_left
        qry_len = qry_end + hard_clipped_right
        qry_start += soft_clipped_left
        qry_end -= soft_clipped_right

        return (trg_start, trg_end, len(ctg_str), trg_seq,
                qry_start, qry_end, qry_len, qry_seq, err_rate)

    def trim_and_transpose(_self, alignmens, region_start, region_end):
        """
        Transforms alignments so that the are strictly within the interval,
        and shifts the coordinates relative to this interval
        """
        MIN_ALN = 100

        trimmed_aln = []
        for aln in alignmens:
            if aln.trg_start >= region_start and aln.trg_end <= region_end:
                trimmed_aln.append(copy(aln))
                continue

            #trimming from left
            new_qry_start = aln.qry_start
            new_trg_start = aln.trg_start
            left_offset = None
            for left_offset in range(len(aln.trg_seq)):
                if new_trg_start >= region_start:
                    break
                if aln.trg_seq[left_offset] != "-":
                    new_trg_start += 1
                if aln.qry_seq[left_offset] != "-":
                    new_qry_start += 1

            #trimming from right
            new_qry_end = aln.qry_end
            new_trg_end = aln.trg_end
            right_offset = None
            for right_offset in range(len(aln.trg_seq)):
                if new_trg_end <= region_end:
                    break
                if aln.trg_seq[-1 - right_offset] != "-":
                    new_trg_end -= 1
                if aln.qry_seq[-1 - right_offset] != "-":
                    new_qry_end -= 1

            if new_trg_end - new_qry_end > MIN_ALN:
                new_qry_seq = aln.qry_seq[left_offset : len(aln.qry_seq) - right_offset]
                new_trg_seq = aln.trg_seq[left_offset : len(aln.trg_seq) - right_offset]
                trimmed_aln.append(aln._replace(qry_start=new_qry_start, qry_end=new_qry_end,
                                                trg_start=new_trg_start, trg_end=new_trg_end,
                                                qry_seq=new_qry_seq, trg_seq=new_trg_seq))

            #print("Aln trg", aln.trg_start, aln.trg_end, "qry", aln.qry_start, aln.qry_end)
            #print("Left offset", left_offset, "right offset", right_offset)
            #print("New aln", new_trg_start, new_trg_end, new_qry_start, new_qry_end)
            #print("")

        for i, aln in enumerate(trimmed_aln):
            trimmed_aln[i] = aln._replace(trg_start=aln.trg_start - region_start,
                                          trg_end=aln.trg_end - region_start,
                                          trg_len=region_end - region_start)

        #print(len(alignmens), len(trimmed_aln))

        return trimmed_aln

    def get_median_depth(self, region_id, region_start=None, region_end=None):
        parsed_contig = _BYTES(region_id)
        contig_str = self.ref_fasta[parsed_contig]
        if region_start is None:
            region_start = 0
        if region_end is None:
            region_end = len(contig_str)

        samtools_out = subprocess.Popen("{0} depth {1} -r '{2}:{3}-{4}' -a -m 0 -Q 10 -l 100"
                                        .format(SAMTOOLS_BIN, self.aln_path,
                                                _STR(parsed_contig), region_start, region_end),
                                        shell=True, stdout=subprocess.PIPE).stdout
        all_cov_pos = []
        for line in samtools_out:
            _ctg, _pos, coverage = line.split()
            all_cov_pos.append(float(coverage))

        return get_median(all_cov_pos) if all_cov_pos else 0

    def get_alignments(self, region_id, region_start=None, region_end=None):
        parsed_contig = _BYTES(region_id)
        contig_str = self.ref_fasta[parsed_contig]
        if region_start is None:
            region_start = 0
        if region_end is None:
            region_end = len(contig_str)
        #logger.debug("Reading region: {0} {1} {2}".format(region_id, region_start, region_end))

        aln_file = subprocess.Popen("{0} view {1} '{2}:{3}-{4}'"
                                        .format(SAMTOOLS_BIN, self.aln_path,
                                                _STR(parsed_contig), region_start, region_end),
                                    shell=True, stdout=subprocess.PIPE).stdout

        chunk_buffer = []
        for line in aln_file:
            chunk_buffer.append(line)

        #shuffle alignments so that they uniformly distributed. Needed for
        #max_coverage subsampling. Using the same seed for determinism
        random.Random(42).shuffle(chunk_buffer)

        ###

        sequence_length = 0
        alignments = []
        for line in chunk_buffer:
            tokens = line.strip().split()
            if len(tokens) < 11:
                #raise AlignmentException("Error reading SAM file")
                continue

            flags = int(tokens[1])
            is_unmapped = flags & 0x4
            is_secondary = flags & 0x100
            is_supplementary = flags & 0x800
            is_reversed = flags & 0x16

            if is_unmapped: continue
            if is_secondary and not self.use_secondary: continue

            read_id = tokens[0]
            cigar_str = tokens[5]
            read_str = tokens[9]
            map_qv = int(tokens[4])
            ctg_pos = int(tokens[3])

            if read_str == b"*":
                continue
                #raise Exception("Error parsing SAM: record without read sequence")

            (trg_start, trg_end, trg_len, trg_seq,
            qry_start, qry_end, qry_len, qry_seq, err_rate) = \
                    self._parse_cigar(cigar_str, read_str, contig_str, ctg_pos)

            #OVERHANG = cfg.vals["read_aln_overhang"]
            #if (float(qry_end - qry_start) / qry_len > self.min_aln_rate or
            #        trg_start < OVERHANG or trg_len - trg_end < OVERHANG):
            aln = Alignment(_STR(read_id), _STR(parsed_contig),
                            qry_start, qry_end, "-" if is_reversed else "+", qry_len,
                            trg_start, trg_end, "+", trg_len,
                            _STR(qry_seq), _STR(trg_seq), err_rate,
                            is_secondary, is_supplementary, map_qv)
            alignments.append(aln)

            sequence_length += qry_end - qry_start
            if sequence_length // len(contig_str) > self.max_coverage:
                break

        #finally, sort alignments by read and by score
        alignments.sort(key=lambda a: (a.qry_id, -(a.qry_end - a.qry_start)))

        return alignments

    def get_all_alignments(self):
        aln_file = subprocess.Popen("{0} view {1}".format(SAMTOOLS_BIN, self.aln_path),
                                    shell=True, stdout=subprocess.PIPE).stdout

        alignments = []
        for line in aln_file:
            tokens = line.strip().split()
            if len(tokens) < 11:
                #raise AlignmentException("Error reading SAM file")
                continue

            flags = int(tokens[1])
            is_unmapped = flags & 0x4
            is_secondary = flags & 0x100
            is_supplementary = flags & 0x800
            is_reversed = flags & 0x16

            if is_unmapped: continue
            if is_secondary and not self.use_secondary: continue

            read_id = tokens[0]
            parsed_contig = tokens[2]
            cigar_str = tokens[5]
            read_str = tokens[9]
            map_qv = int(tokens[4])
            ctg_pos = int(tokens[3])

            if read_str == b"*":
                continue
                #raise Exception("Error parsing SAM: record without read sequence")

            contig_str = self.ref_fasta[parsed_contig]

            (trg_start, trg_end, trg_len, trg_seq,
            qry_start, qry_end, qry_len, qry_seq, err_rate) = \
                    self._parse_cigar(cigar_str, read_str, contig_str, ctg_pos)

            aln = Alignment(_STR(read_id), _STR(parsed_contig),
                            qry_start, qry_end, "-" if is_reversed else "+", qry_len,
                            trg_start, trg_end, "+", trg_len,
                            _STR(qry_seq), _STR(trg_seq), err_rate,
                            is_secondary, is_supplementary, map_qv)
            alignments.append(aln)

        return alignments


#def _is_sam_header(line):
#    return line[:3] in [b"@PG", b"@HD", b"@SQ", b"@RG", b"@CO"]
