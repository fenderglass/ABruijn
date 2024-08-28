#!/usr/bin/env python

#(c) 2019 by Authors
#This file is a part of the Flye package.
#Released under the BSD license (see LICENSE file)

"""
Runs simple toy test
"""

import os
import sys
import subprocess
import shutil
from distutils.spawn import find_executable


def test_toy():
    if not find_executable("flye"):
        sys.exit("flye is not installed!")

    print("Running toy test:\n")
    script_dir = os.path.dirname(os.path.realpath(__file__))
    reads_file = os.path.join(script_dir, "data", "ecoli_500kb_reads_hifi.fastq.gz")
    out_dir = "flye_toy_test"
    subprocess.check_call(["flye", "--pacbio-corr", reads_file, "-g", "500k",
                           "-o", out_dir, "-t", "8", "-m", "1000"])
    shutil.rmtree(out_dir)
    print("\nTEST SUCCESSFUL")


def main():
    test_toy()
    return 0


if __name__ == "__main__":
    sys.exit(main())
