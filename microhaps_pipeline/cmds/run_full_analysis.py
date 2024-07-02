
import os
import pathlib
from ngs_pipeline import cerr, cexit, get_snakefile_path, snakeutils
from ngs_pipeline.cmds import run_snakefile

# this is a wrapper to run initialize.smk

def init_argparser():

    p = snakeutils.init_argparser('run microhaps full analysis')


    m = p.add_mutually_exclusive_group()
    m.add_argument('--single', default=False, action='store_true',
                   help='fastq files are single (non-paired) such as ONT reads')
    m.add_argument('--paired', default=False, action='store_true',
                   help='fastq files are paired such as Illumina paired-end reads')

    p.add_argument('-u', '--underscore', default=0, type=int,
                   help='number of undercore character to be stripped, counted in reverse')
    
    p.add_argument('-o', '--outdir',
                   help='outdir')
    p.add_argument('infiles', nargs='+')
    return p


def run_full_analysis(args):

    import microhaps_pipeline

    args.snakefile = get_snakefile_path('Amplicon_Pipeline.smk', from_module=microhaps_pipeline)
    args.no_config_cascade = True
    args.force = True

    config = dict(
        infiles=args.infiles,
        singleton=args.single,
        paired_end=args.paired,
        underscore=args.underscore,
        outdir=args.outdir
    )

    status, elapsed_time = snakeutils.run_snakefile(args, config=config)

    if not status:
        cerr('[WARNING: run-full-analysis did not successfully complete]')
    cerr(f'[Finish run-full-analysis (time: {elapsed_time})]')

def main(args):
    run_full_analysis(args)


# EOF

