# run_sample_microhaps_caller.py - Microhaps command
# [https://github.com/vivaxgen/Microhaps]

__author__ = "Hidayat Trimarsanto"
__copyright__ = "(C) 2024, Hidayat Trimarsanto"
__email__ = "trimarsanto@gmail.com,hidayat.trimarsanto@menzies.edu.au"
__license__ = "MIT"

# to improve the responsiveness during bash autocomplete, do not import heavy
# modules (such as numpy, pandas, etc) here, but instead import them within the
# functions that require the respective heavy modules

import sys
import os
import pathlib
from ngs_pipeline import cerr, cexit, check_multiplexer
from ngs_pipeline.cmds import run_snakefile
from glob import glob

basedir = os.environ.get("MICROHAPS_BASEDIR", None)
avail_panels = [os.path.basename(panel).replace(".yaml","") for panel in
                    glob(pathlib.posixpath.join(basedir, "configs", "*.yaml"))]

def init_argparser():
    p = run_snakefile.init_argparser("run microhaplotype caller per sample")

    m = p.add_mutually_exclusive_group()
    m.add_argument(
        "--single",
        default=False,
        action="store_true",
        help="fastq files are single (non-paired) such as ONT reads",
    )
    m.add_argument(
        "--paired",
        default=False,
        action="store_true",
        help="fastq files are paired such as Illumina paired-end reads",
    )

    p.add_argument(
        "-u",
        "--underscore",
        default=0,
        type=int,
        help="number of underscore character to be stripped, counted in reverse (see docs)",
    )

    p.add_argument(
        "--skip",
        default=["Undetermined"],
        action="append",
        help="skip samples with the given name (can be used multiple times)",
    )

    p.add_argument(
        "--no-skip", default=False, action="store_true", help="do not skip any samples"
    )

    p.arg_dict["panel"].help = f"the panel for this run. Available panels: {', '.join(avail_panels)}"

    p.add_argument(
        "--run-discovery",
        default=False,
        action="store_true",
        help="run discovery mode (default: False)",
    )


    p.add_argument(
        "--illumina-2-dye",
        default=False,
        action="store_true",
        help="data is from Illumina 2 dye instruments: NovaSeq, NextSeq, MiniSeq",
    )
    p.add_argument("-o", "--outdir", default="output-dir", help="outdir")
    p.add_argument("infiles", nargs="+")

    return p


def run_microhaps_caller(args):

    # check panel
    if not args.panel:
        cexit("Please provide a panel to use using --panel argument")

    # check multiplexer
    check_multiplexer(args)

    # check input files
    for infile in args.infiles:
        if not os.path.exists(infile):
            cexit(f"error: input file {infile} not found")

    # check skip
    if args.no_skip:
        args.skip = []

    args.snakefile = run_snakefile.get_snakefile_path(
        "microhaps_pipeline::msf_microhaps_caller.smk"
    )
    args.no_config_cascade = True
    args.force = True

    args.paired = True if not args.single else False

    config = dict(
        infiles=args.infiles,
        singleton=args.single,
        paired=args.paired,
        underscore=args.underscore,
        outdir=pathlib.Path(args.outdir).absolute().as_posix(),
        skip_list=args.skip if not args.no_skip else [],
        # use generic 2-dye instrument
        instrument="nextseq" if args.illumina_2_dye else "generic",
        # run discovery mode
        joint_discovery=args.run_discovery,
        gatk_drag_haplotypecaller="gatk_drag_haplotypecaller",
        sample_variant_caller_target="all_no_qc",
    )

    status, elapsed_time = run_snakefile.run_snakefile(
        args, config=config, show_status=False
    )

    if not status:
        cerr("[WARNING: run-full-analysis did not successfully complete]")
    cerr(f"[Finish run-full-analysis (time: {elapsed_time})]")


def main(args):
    run_microhaps_caller(args)


#  EOF
