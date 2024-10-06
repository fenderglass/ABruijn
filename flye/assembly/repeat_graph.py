#(c) 2016 by Authors
#This file is a part of ABruijn program.
#Released under the BSD license (see LICENSE file)

"""
Runs repeat/contigger binary
"""

import subprocess
import logging
import os

from flye.utils.utils import which

REPEAT_BIN = "flye-modules"
CONTIGGER_BIN = "flye-modules"
logger = logging.getLogger()


class RepeatException(Exception):
    pass


def check_binaries():
    if not which(REPEAT_BIN) or not which(CONTIGGER_BIN):
        raise RepeatException("Repeat/contigger binaries were not found. "
                              "Did you run 'make'?")
    try:
        devnull = open(os.devnull, "w")
        subprocess.check_call([REPEAT_BIN, "repeat", "-h"], stderr=devnull)
    except subprocess.CalledProcessError as e:
        raise RepeatException(str(e))
    except OSError as e:
        raise RepeatException(str(e))


def analyse_repeats(args, run_params, input_assembly, out_folder,
                    log_file, config_file):
    logger.debug("-----Begin repeat analyser log------")

    cmdline = [REPEAT_BIN, "repeat", "--disjointigs", input_assembly,
               "--reads", ",".join(args.reads), "--out-dir", out_folder,
               "--config", config_file, "--log", log_file,
               "--threads", str(args.threads)]
    if args.debug:
        cmdline.append("--debug")
    if args.meta:
        cmdline.append("--meta")
    if args.keep_haplotypes:
        cmdline.append("--keep-haplotypes")
    #if args.kmer_size:
    #    cmdline.extend(["--kmer", str(args.kmer_size)])
    cmdline.extend(["--min-ovlp", str(run_params["min_overlap"])])

    if args.extra_params:
        cmdline.extend(["--extra-params", args.extra_params])

    try:
        logger.debug("Running: " + " ".join(cmdline))
        subprocess.check_call(cmdline)
    except subprocess.CalledProcessError as e:
        if e.returncode == -9:
            logger.error("Looks like the system ran out of memory")
        raise RepeatException(str(e))
    except OSError as e:
        raise RepeatException(str(e))


def generate_contigs(args, run_params, graph_edges, out_folder,
                    log_file, config_file, repeat_graph, reads_alignment):
    logger.debug("-----Begin contigger analyser log------")

    cmdline = [CONTIGGER_BIN, "contigger", "--graph-edges", graph_edges,
               "--reads", ",".join(args.reads), "--out-dir", out_folder,
               "--config", config_file, "--repeat-graph", repeat_graph,
               "--graph-aln", reads_alignment, "--log", log_file,
               "--threads", str(args.threads)]
    if args.debug:
        cmdline.append("--debug")
    if args.keep_haplotypes:
        cmdline.append("--no-scaffold")
    #if args.kmer_size:
    #    cmdline.extend(["--kmer", str(args.kmer_size)])
    cmdline.extend(["--min-ovlp", str(run_params["min_overlap"])])

    if args.extra_params:
        cmdline.extend(["--extra-params", args.extra_params])

    try:
        logger.debug("Running: " + " ".join(cmdline))
        subprocess.check_call(cmdline)
    except subprocess.CalledProcessError as e:
        if e.returncode == -9:
            logger.error("Looks like the system ran out of memory")
        raise RepeatException(str(e))
    except OSError as e:
        raise RepeatException(str(e))
