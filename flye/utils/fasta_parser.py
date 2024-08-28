#(c) 2013-2016 by Authors
#This file is a part of ABruijn program.
#Released under the BSD license (see LICENSE file)

"""
This module provides some basic FASTA I/O
"""

import logging
import gzip
import io
import sys

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


logger = logging.getLogger()


class FastaError(Exception):
    pass


#Imported functions: take and return unicode strings

def read_sequence_dict(filename):
    """
    Reads Fasta/q file (could be gzip'ed) into a dictionary
    """
    seq_dict = {}
    for hdr, seq in stream_sequence(filename):
        seq_dict[hdr] = seq
    return seq_dict


def read_sequence_lengths(filename):
    seq_dict = {}
    for hdr, seq in stream_sequence(filename):
        seq_dict[hdr] = len(seq)
    return seq_dict


def stream_sequence(filename):
    try:
        gzipped, fastq = _is_fastq(filename)

        if not gzipped:
            handle = open(filename, "rb")
        else:
            #handle = os.popen("gunzip -c {0}".format(filename))
            gz = gzip.open(filename, "rb")
            handle = io.BufferedReader(gz)

        if fastq:
            for hdr, seq, _ in _read_fastq(handle):
                if not _validate_seq(seq):
                    raise FastaError("Invalid char while reading {0}"
                                     .format(filename))
                yield _STR(hdr), _STR(_to_acgt_bytes(seq))
        else:
            for hdr, seq in _read_fasta(handle):
                if not _validate_seq(seq):
                    raise FastaError("Invalid char while reading {0}"
                                     .format(filename))
                yield _STR(hdr), _STR(_to_acgt_bytes(seq))

    except IOError as e:
        raise FastaError(e)


def write_fasta_dict(fasta_dict, filename):
    """
    Writes dictionary with fasta to file
    """
    with open(filename, "w") as f:
        for header in sorted(fasta_dict):
            f.write(">{0}\n".format(header))

            for i in range(0, len(fasta_dict[header]), 60):
                f.write(fasta_dict[header][i:i + 60] + "\n")


def reverse_complement(unicode_str):
    return _STR(reverse_complement_bytes(_BYTES(unicode_str)))


def reverse_complement_bytes(bytes_str):
    return bytes_str.translate(reverse_complement_bytes.COMPL)[::-1]
reverse_complement_bytes.COMPL = maketrans(b"ATGCURYKMSWBVDHNXatgcurykmswbvdhnx",
                                           b"TACGAYRMKSWVBHDNXtacgayrmkswvbhdnx")


def to_acgt(unicode_str):
    return _STR(_to_acgt_bytes(_BYTES(unicode_str)))


#Internal functions: use bytes for faster operations

def _is_fastq(filename):
    suffix = filename.rsplit(".")[-1]
    without_gz = filename
    gzipped = False

    if suffix == "gz":
        gzipped = True
        without_gz = filename.rstrip(".gz")

    suffix = without_gz.rsplit(".")[-1]
    if suffix in ["fasta", "fa"]:
        return gzipped, False

    if suffix in ["fastq", "fq"]:
        return gzipped, True

    raise FastaError("Unknown file extension: " + filename)


def _read_fasta(file_handle):
    """
    bytes input / output
    """
    header = None
    seq = []

    for line in file_handle:
        line = line.strip()
        if not line:
            continue

        if line.startswith(b">"):
            if header:
                yield header, b"".join(seq)
                seq = []
            header = line[1:].split()[0]
        else:
            seq.append(line)

    if header and len(seq):
        yield header, b"".join(seq)


def _read_fastq(file_handle):
    """
    bytes input / output
    """
    seq = None
    qual = None
    header = None
    state_counter = 0

    for no, line in enumerate(file_handle):
        line = line.strip()
        if not line:
            continue

        if state_counter == 0:
            if line[0 : 1] != b"@":
                raise FastaError("Fastq format error: {0} at line {1}"
                                    .format(file_handle.name, no))
            header = line[1:].split()[0]

        if state_counter == 1:
            seq = line

        if state_counter == 2:
            if line[0 : 1] != b"+":
                raise FastaError("Fastq format error: {0} at line {1}"
                                    .format(file_handle.name, no))

        if state_counter == 3:
            qual = line
            yield header, seq, qual

        state_counter = (state_counter + 1) % 4


def _validate_seq(sequence):
    """
    sequence : bytes
    """
    #if len(sequence.strip(_validate_seq.VALID_CHARS)) > 0:
    #    return False
    if len(sequence.translate(None, _validate_seq.VALID_CHARS)):
        return False
    return True
_validate_seq.VALID_CHARS = b"ACGTURYKMSWBDHVNXatgcurykmswbvdhnx"


def _to_acgt_bytes(dna_str):
    """
    assumes tha all characters are valid.
    dna_str : bytes
    """
    if len(dna_str.translate(None, _to_acgt_bytes.ACGT_CHARS)) == 0:
        return dna_str
    #if len(dna_str.strip(_to_acgt_bytes.ACGT_CHARS)) == 0:
    #    return dna_str
    else:
        if not _to_acgt_bytes.ACGT_WARN:
            _to_acgt_bytes.ACGT_WARN = True
            logger.warning("Input contain non-ACGT characters - "
                           "they will be converted to arbitrary ACGTs")
        return dna_str.translate(_to_acgt_bytes.TO_ACGT)
_to_acgt_bytes.ACGT_WARN = False
_to_acgt_bytes.ACGT_CHARS = b"ACGTacgt"
_to_acgt_bytes.TO_ACGT = maketrans(b"URYKMSWBVDHNXurykmswbvdhnx",
                                   b"ACGTACGTACGTAacgtacgtacgta")
