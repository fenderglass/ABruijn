#(c) 2016 by Authors
#This file is a part of ABruijn program.
#Released under the BSD license (see LICENSE file)

"""
Main logic of the package
"""

import sys
import os
import logging
import argparse
import json
import shutil
import subprocess

import flye.polishing.alignment as aln
import flye.polishing.polish as pol
import flye.polishing.consensus as cons
import flye.assembly.assemble as asm
import flye.assembly.repeat_graph as repeat
import flye.assembly.scaffolder as scf
from flye.__version__ import __version__
from flye.__build__ import __build__
import flye.config.py_cfg as cfg
from flye.config.configurator import setup_params, ConfigException
from flye.utils.bytes2human import human2bytes, bytes2human
from flye.utils.sam_parser import AlignmentException
import flye.utils.fasta_parser as fp
#import flye.short_plasmids.plasmids as plas
#import flye.trestle.trestle as tres
#import flye.trestle.graph_resolver as tres_graph
from flye.repeat_graph.repeat_graph import RepeatGraph

if sys.version_info[0] == 2:
    raise Exception("Python 2 is unsupported")

logger = logging.getLogger()

class ResumeException(Exception):
    pass

class Job(object):
    """
    Describes an abstract list of jobs with persistent
    status that can be resumed
    """
    run_params = {"stage_name" : ""}

    def __init__(self):
        self.name = None
        self.args = None
        self.work_dir = None
        self.out_files = {}
        self.log_file = None

    def run(self):
        logger.info(">>>STAGE: %s", self.name)

    def save(self, save_file):
        Job.run_params["stage_name"] = self.name

        with open(save_file, "w") as f:
            json.dump(Job.run_params, f)

    def load(self, save_file):
        with open(save_file, "r") as f:
            data = json.load(f)
            if (not "pipeline_version" in data or
                    data["pipeline_version"] != cfg.vals["pipeline_version"]):
                raise ResumeException("Inconsistent pipeline version")

            Job.run_params = data

    def completed(self, save_file):
        with open(save_file, "r") as f:
            dummy_data = json.load(f)

            for file_path in self.out_files.values():
                if not os.path.exists(file_path):
                    return False

            return True


class JobConfigure(Job):
    def __init__(self, args, work_dir):
        super(JobConfigure, self).__init__()
        self.args = args
        self.work_dir = work_dir
        self.name = "configure"

    def run(self):
        super(JobConfigure, self).run()
        params = setup_params(self.args)
        Job.run_params = params


class JobAssembly(Job):
    def __init__(self, args, work_dir, log_file):
        super(JobAssembly, self).__init__()
        #self.out_assembly = out_assembly
        self.args = args
        self.work_dir = work_dir
        self.log_file = log_file

        self.name = "assembly"
        self.assembly_dir = os.path.join(self.work_dir, "00-assembly")
        self.assembly_filename = os.path.join(self.assembly_dir,
                                              "draft_assembly.fasta")
        self.out_files["assembly"] = self.assembly_filename

    def run(self):
        super(JobAssembly, self).run()
        if not os.path.isdir(self.assembly_dir):
            os.mkdir(self.assembly_dir)
        asm.assemble(self.args, Job.run_params, self.assembly_filename,
                     self.log_file, self.args.asm_config, )
        if os.path.getsize(self.assembly_filename) == 0:
            raise asm.AssembleException("No disjointigs were assembled - "
                                        "please check if the read type and genome "
                                        "size parameters are correct")
        asm_len, asm_n50 = scf.short_statistics(self.assembly_filename)
        logger.debug("Disjointigs length: %d, N50: %d", asm_len, asm_n50)


class JobRepeat(Job):
    def __init__(self, args, work_dir, log_file, disjointigs):
        super(JobRepeat, self).__init__()

        self.args = args
        self.disjointigs = disjointigs
        self.log_file = log_file
        self.name = "repeat"

        self.work_dir = os.path.join(work_dir, "20-repeat")
        self.out_files["repeat_graph"] = os.path.join(self.work_dir,
                                                      "repeat_graph_dump")
        self.out_files["repeat_graph_edges"] = os.path.join(self.work_dir,
                                                        "repeat_graph_edges.fasta")
        self.out_files["reads_alignment"] = os.path.join(self.work_dir,
                                                         "read_alignment_dump")
        #self.out_files["repeats_dump"] = os.path.join(self.work_dir,
        #                                              "repeats_dump")

    def run(self):
        super(JobRepeat, self).run()
        if not os.path.isdir(self.work_dir):
            os.mkdir(self.work_dir)
        logger.info("Building and resolving repeat graph")
        repeat.analyse_repeats(self.args, Job.run_params, self.disjointigs,
                               self.work_dir, self.log_file,
                               self.args.asm_config)


class JobContigger(Job):
    def __init__(self, args, work_dir, log_file, repeat_graph_edges,
                 repeat_graph, reads_alignment):
        super(JobContigger, self).__init__()

        self.args = args
        self.repeat_graph_edges = repeat_graph_edges
        self.repeat_graph = repeat_graph
        self.reads_alignment = reads_alignment
        self.log_file = log_file
        self.name = "contigger"

        self.work_dir = os.path.join(work_dir, "30-contigger")
        self.out_files["contigs"] = os.path.join(self.work_dir,
                                                 "contigs.fasta")
        self.out_files["assembly_graph"] = os.path.join(self.work_dir,
                                                        "graph_final.gv")
        self.out_files["edges_sequences"] = os.path.join(self.work_dir,
                                                         "graph_final.fasta")
        self.out_files["gfa_graph"] = os.path.join(self.work_dir,
                                                   "graph_final.gfa")
        self.out_files["stats"] = os.path.join(self.work_dir, "contigs_stats.txt")
        self.out_files["scaffold_links"] = os.path.join(self.work_dir,
                                                        "scaffolds_links.txt")

    def run(self):
        super(JobContigger, self).run()
        if not os.path.isdir(self.work_dir):
            os.mkdir(self.work_dir)
        logger.info("Generating contigs")
        repeat.generate_contigs(self.args, Job.run_params, self.repeat_graph_edges,
                                self.work_dir, self.log_file, self.args.asm_config,
                                self.repeat_graph, self.reads_alignment)

        if os.path.getsize(self.out_files["contigs"]) == 0:
            raise asm.AssembleException("No contigs were assembled - "
                                        "pipeline stopped")



def _list_files(startpath, maxlevel=1):
    for root, _, files in os.walk(startpath):
        level = root.replace(startpath, "").count(os.sep)
        if level > maxlevel:
            continue
        indent = " " * 4 * (level)
        logger.debug(indent + os.path.basename(root) + "/")
        subindent = " " * 4 * (level + 1)
        for f in files:
            fsize = bytes2human(os.path.getsize(os.path.join(root, f)))
            logger.debug("%s%-12s%s", subindent, fsize, f)


class JobFinalize(Job):
    def __init__(self, args, work_dir, log_file,
                 contigs_file, graph_file, repeat_stats,
                 polished_stats, polished_gfa, scaffold_links):
        super(JobFinalize, self).__init__()

        self.args = args
        self.log_file = log_file
        self.name = "finalize"
        self.contigs_file = contigs_file
        self.graph_file = graph_file
        self.repeat_stats = repeat_stats
        self.polished_stats = polished_stats
        self.scaffold_links = scaffold_links
        self.polished_gfa = polished_gfa
        self.work_dir = work_dir

        #self.out_files["contigs"] = os.path.join(work_dir, "contigs.fasta")
        #self.out_files["scaffolds"] = os.path.join(work_dir, "scaffolds.fasta")
        self.out_files["assembly"] = os.path.join(work_dir, "assembly.fasta")
        self.out_files["stats"] = os.path.join(work_dir, "assembly_info.txt")
        self.out_files["graph"] = os.path.join(work_dir, "assembly_graph.gv")
        self.out_files["gfa"] = os.path.join(work_dir, "assembly_graph.gfa")

    def run(self):
        super(JobFinalize, self).run()
        #shutil.copy(self.contigs_file, self.out_files["contigs"])
        shutil.copy(self.graph_file, self.out_files["graph"])
        shutil.copy(self.polished_gfa, self.out_files["gfa"])

        scf_links = self.scaffold_links if self.args.scaffold else None
        scaffolds = scf.generate_scaffolds(self.contigs_file, scf_links,
                                           self.out_files["assembly"])

        #create the scaffolds.fasta symlink for backward compatability
        #try:
        #    if os.path.lexists(self.out_files["scaffolds"]):
        #        os.remove(self.out_files["scaffolds"])
        #    relative_link = os.path.relpath(self.out_files["assembly"],
        #                                    self.work_dir)
        #    os.symlink(relative_link, self.out_files["scaffolds"])
        #except OSError as e:
        #    logger.debug(e)

        logger.debug("---Output dir contents:----")
        _list_files(os.path.abspath(self.args.out_dir))
        logger.debug("--------------------------")
        scf.generate_stats(self.repeat_stats, self.polished_stats, scaffolds,
                           self.out_files["stats"])
        logger.info("Final assembly: %s", self.out_files["assembly"])


class JobConsensus(Job):
    def __init__(self, args, work_dir, in_contigs):
        super(JobConsensus, self).__init__()

        self.args = args
        self.in_contigs = in_contigs
        self.consensus_dir = os.path.join(work_dir, "10-consensus")
        self.out_consensus = os.path.join(self.consensus_dir, "consensus.fasta")
        self.name = "consensus"
        self.out_files["consensus"] = self.out_consensus

    def run(self):
        super(JobConsensus, self).run()
        if not os.path.isdir(self.consensus_dir):
            os.mkdir(self.consensus_dir)

        logger.info("Running Minimap2")
        out_alignment = os.path.join(self.consensus_dir, "minimap.bam")
        aln.make_alignment(self.in_contigs, self.args.reads, self.args.threads,
                           self.args.platform, self.args.read_type, out_alignment)

        contigs_info = aln.get_contigs_info(self.in_contigs)
        logger.info("Computing consensus")
        consensus_fasta = cons.get_consensus(out_alignment, self.in_contigs,
                                             contigs_info, self.args.threads,
                                             self.args.platform)

        fp.write_fasta_dict(consensus_fasta, self.out_consensus)
        os.remove(out_alignment)


class JobPolishing(Job):
    def __init__(self, args, work_dir, log_file, in_contigs, in_graph_edges,
                 in_graph_gfa):
        super(JobPolishing, self).__init__()

        self.args = args
        self.log_file = log_file
        self.in_contigs = in_contigs
        self.in_graph_edges = in_graph_edges
        self.in_graph_gfa = in_graph_gfa
        self.polishing_dir = os.path.join(work_dir, "40-polishing")

        self.name = "polishing"
        final_contigs = os.path.join(self.polishing_dir,
                                     "filtered_contigs.fasta")
        self.out_files["contigs"] = final_contigs
        self.out_files["stats"] = os.path.join(self.polishing_dir,
                                               "filtered_stats.txt")
        self.out_files["polished_gfa"] = os.path.join(self.polishing_dir,
                                                      "polished_edges.gfa")

    def run(self):
        super(JobPolishing, self).run()
        if not os.path.isdir(self.polishing_dir):
            os.mkdir(self.polishing_dir)

        contigs, stats = \
            pol.polish(self.in_contigs, self.args.reads, self.polishing_dir,
                       self.args.num_iters, self.args.threads, self.args.platform,
                       self.args.read_type, output_progress=True)
        #contigs = os.path.join(self.polishing_dir, "polished_1.fasta")
        #stats = os.path.join(self.polishing_dir, "contigs_stats.txt")
        pol.filter_by_coverage(self.args, stats, contigs,
                               self.out_files["stats"], self.out_files["contigs"])
        pol.generate_polished_edges(self.in_graph_edges, self.in_graph_gfa,
                                    self.out_files["contigs"],
                                    self.polishing_dir, self.args.platform, self.args.read_type,
                                    stats, self.args.threads)
        os.remove(contigs)
        if os.path.getsize(self.out_files["contigs"]) == 0:
            raise asm.AssembleException("No contigs were assembled - "
                                        "pipeline stopped")


def _create_job_list(args, work_dir, log_file):
    """
    Build pipeline as a list of consecutive jobs
    """
    jobs = []

    #Run configuration
    jobs.append(JobConfigure(args, work_dir))

    #Assembly job
    jobs.append(JobAssembly(args, work_dir, log_file))
    disjointigs = jobs[-1].out_files["assembly"]

    #Consensus
    if args.read_type != "subasm":
        jobs.append(JobConsensus(args, work_dir, disjointigs))
        disjointigs = jobs[-1].out_files["consensus"]

    #Repeat analysis
    jobs.append(JobRepeat(args, work_dir, log_file, disjointigs))
    #repeats_dump = jobs[-1].out_files["repeats_dump"]
    repeat_graph_edges = jobs[-1].out_files["repeat_graph_edges"]
    repeat_graph = jobs[-1].out_files["repeat_graph"]
    reads_alignment = jobs[-1].out_files["reads_alignment"]

    #Trestle: Resolve Unbridged Repeats
    #if not args.no_trestle and not args.meta and args.read_type == "raw":
    #if args.trestle:
    #    jobs.append(JobTrestle(args, work_dir, log_file,
    #                repeat_graph, repeat_graph_edges,
    #                reads_alignment))
    #    repeat_graph_edges = jobs[-1].out_files["repeat_graph_edges"]
    #    repeat_graph = jobs[-1].out_files["repeat_graph"]

    #Contigger
    jobs.append(JobContigger(args, work_dir, log_file, repeat_graph_edges,
                             repeat_graph, reads_alignment))
    raw_contigs = jobs[-1].out_files["contigs"]
    scaffold_links = jobs[-1].out_files["scaffold_links"]
    graph_file = jobs[-1].out_files["assembly_graph"]
    gfa_file = jobs[-1].out_files["gfa_graph"]
    final_graph_edges = jobs[-1].out_files["edges_sequences"]
    repeat_stats = jobs[-1].out_files["stats"]

    #Polishing
    contigs_file = raw_contigs
    polished_stats = None
    polished_gfa = gfa_file
    if args.num_iters > 0:
        jobs.append(JobPolishing(args, work_dir, log_file, raw_contigs,
                                 final_graph_edges, gfa_file))
        contigs_file = jobs[-1].out_files["contigs"]
        polished_stats = jobs[-1].out_files["stats"]
        polished_gfa = jobs[-1].out_files["polished_gfa"]

    #Report results
    jobs.append(JobFinalize(args, work_dir, log_file, contigs_file,
                            graph_file, repeat_stats, polished_stats,
                            polished_gfa, scaffold_links))

    return jobs


def _set_genome_size(args):
    """
    Select k-mer size based on the target genome size
    """
    if args.genome_size.isdigit():
        args.genome_size = int(args.genome_size)
    else:
        args.genome_size = human2bytes(args.genome_size.upper())


def _run_polisher_only(args, output_progress=True):
    """
    Runs standalone polisher
    """
    if output_progress:
        logger.info("Running Flye polisher")
        logger.debug("Cmd: %s", " ".join(sys.argv))
    bam_input = False

    for read_file in args.reads:
        if not os.path.exists(read_file):
            raise ResumeException("Can't open " + read_file)
        if " " in read_file:
            raise ResumeException("Path to reads contain spaces: " + read_file)

        without_gz = read_file.rstrip(".gz")
        if not any([without_gz.endswith(x) for x in ["fasta", "fa", "fastq", "fq", "bam"]]):
            raise ResumeException("Unsupported input. Supported types: fasta/fastq/bam")
        if without_gz.endswith("bam"):
            bam_input = True

    if bam_input and len(args.reads) > 1:
        raise ResumeException("Only single bam input supported")

    if bam_input and args.num_iters > 1:
        raise ResumeException("Bam input only supports single iteration. For multiple iterations, provide fastq instead")

    pol.polish(args.polish_target, args.reads, args.out_dir,
               args.num_iters, args.threads, args.platform,
               args.read_type, output_progress)
    if output_progress:
        logger.info("Done!")


def _run(args):
    """
    Runs the pipeline
    """
    logger.info("Starting Flye " + _version())
    logger.debug("Cmd: %s", " ".join(sys.argv))
    logger.debug("Python version: " + sys.version)

    if args.genome_size:
        _set_genome_size(args)

    for read_file in args.reads:
        if not os.path.exists(read_file):
            raise ResumeException("Can't open " + read_file)
        if " " in read_file:
            raise ResumeException("Path to reads contain spaces: " + read_file)

    save_file = os.path.join(args.out_dir, "params.json")
    jobs = _create_job_list(args, args.out_dir, args.log_file)

    if args.stop_after and not args.stop_after in [j.name for j in jobs]:
        raise ResumeException("Stop after: unknown stage '{0}'"
                                .format(args.stop_after))

    current_job = 0
    if args.resume or args.resume_from:
        if not os.path.exists(save_file):
            raise ResumeException("Can't find save file")

        logger.info("Resuming previous run")
        if args.resume_from:
            job_to_resume = args.resume_from
        else:
            job_to_resume = json.load(open(save_file, "r"))["stage_name"]

        can_resume = False
        for i in range(len(jobs)):
            if jobs[i].name == job_to_resume:
                jobs[i].load(save_file)
                current_job = i
                if not jobs[i - 1].completed(save_file):
                    raise ResumeException("Can't resume: stage '{0}' incomplete"
                                          .format(jobs[i - 1].name))
                can_resume = True
                break

        if not can_resume:
            raise ResumeException("Can't resume: stage {0} does not exist"
                                  .format(job_to_resume))

    for i in range(current_job, len(jobs)):
        jobs[i].save(save_file)
        jobs[i].run()
        if args.stop_after == jobs[i].name:
            if i + 1 < len(jobs):
                jobs[i + 1].save(save_file)
            logger.info("Pipeline stopped as requested by --stop-after")
            break


def _enable_logging(log_file, debug, overwrite):
    """
    Turns on logging, sets debug levels and assigns a log file
    """
    log_formatter = logging.Formatter("[%(asctime)s] %(name)s: %(levelname)s: "
                                      "%(message)s", "%Y-%m-%d %H:%M:%S")
    console_formatter = logging.Formatter("[%(asctime)s] %(levelname)s: "
                                          "%(message)s", "%Y-%m-%d %H:%M:%S")
    console_log = logging.StreamHandler()
    console_log.setFormatter(console_formatter)
    if not debug:
        console_log.setLevel(logging.INFO)

    if overwrite:
        open(log_file, "w").close()
    file_handler = logging.FileHandler(log_file, mode="a")
    file_handler.setFormatter(log_formatter)

    logger.setLevel(logging.DEBUG)
    logger.addHandler(console_log)
    logger.addHandler(file_handler)


def _usage():
    return ("flye (--pacbio-raw | --pacbio-corr | --pacbio-hifi | --nano-raw |\n"
            "\t     --nano-corr | --nano-hq ) file1 [file_2 ...]\n"
            "\t     --out-dir PATH\n\n"
            "\t     [--genome-size SIZE] [--threads int] [--iterations int]\n"
            "\t     [--meta] [--polish-target] [--min-overlap SIZE]\n"
            "\t     [--keep-haplotypes] [--debug] [--version] [--help] \n"
            "\t     [--scaffold] [--resume] [--resume-from] [--stop-after] \n"
            "\t     [--read-error float] [--extra-params] \n"
            "\t     [--deterministic]")


def _epilog():
    return ("Input reads can be in FASTA or FASTQ format, uncompressed\n"
            "or compressed with gz. Currently, PacBio (CLR, HiFi, corrected)\n"
            "and ONT reads (regular, HQ, corrected) are supported. Expected error rates are\n"
            "<15% for PB CLR/regular ONT; <5% for ONT HQ, <3% for corrected, and <1% for HiFi. Note that Flye\n"
            "was primarily developed to run on uncorrected reads. You may specify multiple\n"
            "files with reads (separated by spaces). Mixing different read\n"
            "types is not yet supported. The --meta option enables the mode\n"
            "for metagenome/uneven coverage assembly.\n\n"
            "To reduce memory consumption for large genome assemblies,\n"
            "you can use a subset of the longest reads for initial disjointig\n"
            "assembly by specifying --asm-coverage and --genome-size options. Typically,\n"
            "40x coverage is enough to produce good disjointigs.\n\n"
            "You can run Flye polisher as a standalone tool using\n"
            "--polish-target option.")


def _version():
    return __version__ + "-b" + str(__build__)
    """
    repo_root = os.path.dirname((os.path.dirname(__file__)))
    try:
        git_label = subprocess.check_output(["git", "-C", repo_root, "describe"],
                                            stderr=open(os.devnull, "w")).decode()
        commit_id = git_label.strip("\n").rsplit("-", 1)[-1]
        return __version__ + "-" + commit_id
    except (subprocess.CalledProcessError, OSError):
        pass
    return __version__ + "-release"
    """


def main():
    def check_int_range(value, min_val, max_val, require_odd=False):
        ival = int(value)
        if ival < min_val or ival > max_val:
            raise argparse.ArgumentTypeError("value should be in the "
                            "range [{0}, {1}]".format(min_val, max_val))
        if require_odd and ival % 2 == 0:
            raise argparse.ArgumentTypeError("should be an odd number")
        return ival

    parser = argparse.ArgumentParser \
        (description="Assembly of long reads with repeat graphs",
         formatter_class=argparse.RawDescriptionHelpFormatter,
         usage=_usage(), epilog=_epilog())

    read_group = parser.add_mutually_exclusive_group(required=True)
    read_group.add_argument("--pacbio-raw", dest="pacbio_raw",
                        default=None, metavar="path", nargs="+",
                        help="PacBio regular CLR reads (<20%% error)")
    read_group.add_argument("--pacbio-corr", dest="pacbio_corrected",
                        default=None, metavar="path", nargs="+",
                        help="PacBio reads that were corrected with other methods (<3%% error)")
    read_group.add_argument("--pacbio-hifi", dest="pacbio_hifi",
                        default=None, metavar="path", nargs="+",
                        help="PacBio HiFi reads (<1%% error)")
    read_group.add_argument("--nano-raw", dest="nano_raw", nargs="+",
                        default=None, metavar="path",
                        help="ONT regular reads, pre-Guppy5 (<20%% error)")
    read_group.add_argument("--nano-corr", dest="nano_corrected", nargs="+",
                        default=None, metavar="path",
                        help="ONT reads that were corrected with other methods (<3%% error)")
    read_group.add_argument("--nano-hq", dest="nano_hq", nargs="+",
                        default=None, metavar="path",
                        help="ONT high-quality reads: Guppy5+ SUP or Q20 (<5%% error)")
    read_group.add_argument("--subassemblies", dest="subassemblies", nargs="+",
                        default=None, metavar="path",
                        help="[deprecated] high-quality contigs input")
    parser.add_argument("-g", "--genome-size", dest="genome_size",
                        metavar="size", required=False, default=None,
                        help="estimated genome size (for example, 5m or 2.6g)")
    parser.add_argument("-o", "--out-dir", dest="out_dir",
                        default=None, required=True,
                        metavar="path", help="Output directory")
    parser.add_argument("-t", "--threads", dest="threads",
                        type=lambda v: check_int_range(v, 1, 128),
                        default=1, metavar="int", help="number of parallel threads [1]")
    parser.add_argument("-i", "--iterations", dest="num_iters",
                        type=lambda v: check_int_range(v, 0, 10),
                        default=1, help="number of polishing iterations [1]",
                        metavar="int")
    parser.add_argument("-m", "--min-overlap", dest="min_overlap", metavar="int",
                        type=lambda v: check_int_range(v, 1000, 10000),
                        default=None, help="minimum overlap between reads [auto]")
    parser.add_argument("--asm-coverage", dest="asm_coverage", metavar="int",
                        default=None, help="reduced coverage for initial "
                        "disjointig assembly [not set]", type=int)
    parser.add_argument("--hifi-error", dest="hifi_error", metavar="float",
                        default=None, help="[deprecated] same as --read-error", type=float)
    parser.add_argument("--read-error", dest="read_error", metavar="float",
                        default=None, help="adjust parameters for given read error rate (as fraction e.g. 0.03)", type=float)
    parser.add_argument("--extra-params", dest="extra_params",
                        metavar="extra_params", required=False, default=None,
                        help="extra configuration parameters list (comma-separated)")
    parser.add_argument("--plasmids", action="store_true",
                        dest="plasmids", default=False,
                        help="unused (retained for backward compatibility)")
    parser.add_argument("--meta", action="store_true",
                        dest="meta", default=False,
                        help="metagenome / uneven coverage mode")
    parser.add_argument("--keep-haplotypes", action="store_true",
                        dest="keep_haplotypes", default=False,
                        help="do not collapse alternative haplotypes")
    parser.add_argument("--no-alt-contigs", action="store_true",
                        dest="no_alt_contigs", default=False,
                        help="do not output contigs representing alternative haplotypes")
    parser.add_argument("--scaffold", action="store_true",
                        dest="scaffold", default=False,
                        help="enable scaffolding using graph [disabled by default]")
    parser.add_argument("--trestle", action="store_true",
                        dest="trestle", default=False,
                        help="[deprecated] enable Trestle [disabled by default]")
    parser.add_argument("--polish-target", dest="polish_target",
                        metavar="path", required=False,
                        help="run polisher on the target sequence")
    parser.add_argument("--resume", action="store_true",
                        dest="resume", default=False,
                        help="resume from the last completed stage")
    parser.add_argument("--resume-from", dest="resume_from", metavar="stage_name",
                        default=None, help="resume from a custom stage")
    parser.add_argument("--stop-after", dest="stop_after", metavar="stage_name",
                        default=None, help="stop after the specified stage completed")
    #parser.add_argument("--kmer-size", dest="kmer_size",
    #                    type=lambda v: check_int_range(v, 11, 31, require_odd=True),
    #                    default=None, help="kmer size (default: auto)")
    parser.add_argument("--debug", action="store_true",
                        dest="debug", default=False,
                        help="enable debug output")
    parser.add_argument("-v", "--version", action="version", version=_version())
    parser.add_argument("--deterministic", action="store_true",
                        dest="deterministic", default=False,
                        help="perform disjointig assembly single-threaded")
    args = parser.parse_args()

    if args.asm_coverage and (args.genome_size is None):
        parser.error("--asm-coverage option requires genome size estimate (--genome-size)")

    if args.asm_coverage and args.meta:
        parser.error("--asm-coverage is incompatible with --meta")

    if args.hifi_error and not args.read_error:
        args.read_error = args.hifi_error
    if args.read_error and (args.pacbio_raw or args.nano_raw):
        parser.error("--read-error can only be used with corr/hq/hifi modes")
    if args.read_error and args.read_error > 1:
        parser.error("--read-error expressed as a decimal fraction, e.g. 0.01 or 0.03")

    def _add_extra_param(param):
        if args.extra_params:
            args.extra_params += "," + param
        else:
            args.extra_params = param

    if args.read_error:
        hifi_str = "assemble_ovlp_divergence={0},repeat_graph_ovlp_divergence={0}".format(args.read_error)
        _add_extra_param(hifi_str)

    if args.no_alt_contigs:
        _add_extra_param("remove_alt_edges=1")

    if args.keep_haplotypes:
        _add_extra_param("aggressive_dup_filter=0")

    if args.pacbio_raw:
        args.reads = args.pacbio_raw
        args.platform = "pacbio"
        args.read_type = "raw"
    if args.pacbio_corrected:
        args.reads = args.pacbio_corrected
        args.platform = "pacbio"
        args.read_type = "corrected"
    if args.pacbio_hifi:
        args.reads = args.pacbio_hifi
        args.platform = "pacbio"
        args.read_type = "hifi"
    if args.nano_raw:
        args.reads = args.nano_raw
        args.platform = "nano"
        args.read_type = "raw"
    if args.nano_corrected:
        args.reads = args.nano_corrected
        args.platform = "nano"
        args.read_type = "corrected"
    if args.nano_hq:
        args.reads = args.nano_hq
        args.platform = "nano"
        args.read_type = "nano_hq"
    if args.subassemblies:
        args.reads = args.subassemblies
        args.platform = "pacbio"    #arbitrary
        args.read_type = "subasm"

    if not os.path.isdir(args.out_dir):
        os.mkdir(args.out_dir)
    args.out_dir = os.path.abspath(args.out_dir)

    args.reads = [os.path.abspath(r) for r in args.reads]

    args.log_file = os.path.join(args.out_dir, "flye.log")
    _enable_logging(args.log_file, args.debug,
                    overwrite=False)

    args.asm_config = os.path.join(cfg.vals["pkg_root"],
                                   cfg.vals["bin_cfg"][args.read_type])

    if args.plasmids:
        logger.warning("--plasmids mode is no longer available. Command line option will be removed in the future versions")
    if args.trestle:
        logger.warning("--trestle mode is being deprecated. It will be removed in the future versions.")
    if args.subassemblies:
        logger.warning("--subassemblies mode is being deprecated. It will be removed in the future versions.")

    try:
        aln.check_binaries()
        pol.check_binaries()
        asm.check_binaries()
        repeat.check_binaries()

        if not args.polish_target:
            _run(args)
        else:
            _run_polisher_only(args)

    except (AlignmentException, pol.PolishException,
            asm.AssembleException, repeat.RepeatException,
            ResumeException, fp.FastaError, ConfigException) as e:
        logger.error(e)
        logger.error("Pipeline aborted")
        return 1

    return 0
