#(c) 2016 by Authors
#This file is a part of ABruijn program.
#Released under the BSD license (see LICENSE file)

"""
Configuration file for the Python part of the pipeline
"""

import os

vals = {
        "pkg_root" : os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
        "pipeline_version" : 3,

        #additional configuration files for binary modules
        "bin_cfg" : {
            "raw" : "config/bin_cfg/asm_raw_reads.cfg",
            "corrected" : "config/bin_cfg/asm_corrected_reads.cfg",
            "hifi" : "config/bin_cfg/asm_hifi.cfg",
            "nano_hq" : "config/bin_cfg/asm_nano_hq.cfg",
            "subasm" : "config/bin_cfg/asm_subasm.cfg"
        },

        "min_overlap_range" : {
            "raw" : [1000, 10000],
            "corrected" : [1000, 10000],
            "hifi" : [1000, 10000],
            "nano_hq" : [1000, 10000],
            "subasm" : [1000, 1000]
        },
        "max_meta_overlap" : 10000,

        #polishing
        "simple_kmer_length" : 4,
        "solid_kmer_length" : 10,
        "max_bubble_length" : 500,
        "max_bubble_branches" : 50,
        "max_read_coverage" : 1000,
        "min_polish_aln_len" : 500,

        #final coverage filtering
        "relative_minimum_coverage" : 5,
        "hard_minimum_coverage" : 3,

        "err_modes" : {
            "pacbio" : {
                "subs_matrix" : "config/bin_cfg/pacbio_chm13_substitutions.mat",
                "hopo_matrix" : "config/bin_cfg/pacbio_chm13_homopolymers.mat",
                "solid_missmatch" : 0.3,
                "solid_indel" : 0.3,
                "max_aln_error" : 0.25,
                "hopo_enabled" : False
            },
            "nano" : {
                "subs_matrix" : "config/bin_cfg/nano_r94_substitutions.mat",
                "hopo_matrix" : "config/bin_cfg/nano_r94_g36_homopolymers.mat",
                "solid_missmatch" : 0.3,
                "solid_indel" : 0.3,
                "max_aln_error" : 0.25,
                "hopo_enabled" : False
            },
        },

        #scaffolding
        "scaffold_gap" : 100
    }
