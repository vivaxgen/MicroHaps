#!/usr/bin/env python

import argparse
import subprocess
import threading
import multiprocessing
import sys
import os
import subprocess as sp
#from fastq2matrix import run_cmd

MICROHAPS_BASEDIR = os.path.expanduser("~/tools/MicroHaps/scripts")
#MICROHAPS_BASEDIR = os.environ['MICROHAPS_BASEDIR']

def run_cmd(cmd):
    sys.stderr.write("Running command:\n%s\n\n" % cmd)
    with open("/dev/null","w") as O:
        res = sp.call(cmd,shell=True,stderr=O,stdout=O)
    if res!=0:
        sys.exit("Error running last command, please check!\n")

def trim_primer(sampleid, fileF, fileR, pr1, pr2):
    if os.path.isfile(fileF) and os.path.isfile(fileR):
        cmd = [
            'cutadapt', '-g', ('file:' + pr1), '-G', ('file:' + pr2),
            '-o', os.path.join(run_dir, "prim_fq", sampleid + "_prim_1.fq.gz"),
            '-p', os.path.join(run_dir, "prim_fq", sampleid + "_prim_2.fq.gz"),
            '--pair-adapters', '--discard-untrimmed', '--action=trim',
            fileF, fileR
        ]
        print("Primer Removal Command:", " ".join(cmd))
        proc = subprocess.Popen(cmd, stdout=sys.stdout, stderr=sys.stderr)
        proc.wait()
    else:
        sys.exit('Pre-process halted: one or both of the fastq files not found! Exiting..')


def main():
    global run_dir
    global path
    path = os.path.dirname(__file__)

    parser = argparse.ArgumentParser()
    parser.add_argument('--path_to_meta', help="Path to input fastq files")
    parser.add_argument('--keep_primers', action="store_true", help="Skip primer removal step")
    parser.add_argument('--pr1', help="Path to forward primers FASTA file")
    parser.add_argument('--pr2', help="Path to reverse primers FASTA file")
    parser.add_argument('--Class', help="Specify Analysis class. Accepts one of two: parasite/vector")
    parser.add_argument('--maxEE', help="Maximum Expected errors (dada2 filtering argument)")
    parser.add_argument('--trimRight', help="Hard trim number of bases at 5` end (dada2 filtering argument)")
    parser.add_argument('--minLen', help="Minimum length filter (dada2 filtering argument)")
    parser.add_argument('--truncQ', help="Soft trim bases based on quality (dada2 filtering argument)")
    parser.add_argument('--max_consist', help="Number of cycles for consistency in error model (dada2 argument)")
    parser.add_argument('--omegaA', help="p-value for the partitioning algorithm (dada2 argument)")
    parser.add_argument('--justConcatenate', help="whether reads should be concatenated with N's during merge (dada2 argument)")
    #parser.add_argument('--saveRdata', help="Optionally save dada2 part of this run as Rdata object")

    args = parser.parse_args()

    if args.path_to_meta is not None:
        print('NOTE: --path_to_meta argument found. Using provided path.')
        path_to_meta = args.path_to_meta
    else:
        sys.exit('Execution halted: --path_to_meta not found! Exiting..')

    if os.path.isfile(path_to_meta):
        run_dir = os.path.abspath(os.path.join(path_to_meta, os.pardir))
        sys.stdout = open((run_dir + "/stdout.txt"), "a")
        sys.stderr = open((run_dir + "/stderr.txt"), "a")
    else:
        sys.exit('Execution halted: Metafile with sample list not found! Exiting..')

    # ... (rest of the argument parsing logic)

    if args.keep_primers:
        print("skipping Primer removal step..")
        pass
    else:
        if args.pr1 is not None:
            print('NOTE: --pr1 argument found. Using provided path.')
            pr1 = args.pr1
            mark = True
        else:
            print('NOTE: --pr1 argument not found. Skipping primer removal..')
            mark = False

        if args.pr2 is not None:
            print('NOTE: --pr2 argument found. Using provided path.')
            pr2 = args.pr2
            mark = True
        else:
            print('NOTE: --pr2 argument not found. Skipping primer removal..')
            mark = False

        if mark:
            # Primer removal steps here
            if not os.path.exists(os.path.join(run_dir, "prim_fq")):
                os.mkdir(os.path.join(run_dir, "prim_fq"))
            else:
                print("Directory %s already exists.." % (os.path.join(run_dir, "prim_fq")))

            print("Now running Primer removal..")
            meta = open(path_to_meta, 'r')
            samples = meta.readlines()
            p = multiprocessing.Pool()
            for sample in samples:
                slist = sample.split()
                p.apply_async(trim_primer, args=(slist[0], slist[1], slist[2], pr1, pr2))
            p.close()
            p.join()
            meta.close()

            #create_meta(os.path.join(run_dir, "prim_fq"), os.path.join(run_dir, "prim_meta.txt"),
            #            pattern_fw="*_prim_1.fq.gz", pattern_rv="*_prim_2.fq.gz")
            #path_to_meta = os.path.join(run_dir, "prim_meta.txt")

            run_cmd('create_meta.py --path_to_fq prim_fq --output_file prim_meta.txt --pattern_fw "*_prim_1.fq.gz" --pattern_rv "*_prim_2.fq.gz"' % vars(args))
            

    # Steps after Primer removal
    if not os.path.exists(os.path.join(run_dir, "run_dada2")):
        os.mkdir(os.path.join(run_dir, "run_dada2"))
    else:
        print("Directory %s already exists.." % (os.path.join(run_dir, "run_dada2")))

    print("Now running DADA2..")
    dada2_cmd = [
    'Rscript', f'runDADA2.R',
    '-p', 'prim_meta.txt',
    '-d', os.path.join(run_dir, 'run_dada2'),
    '-o', 'seqtab.tsv',
    '-c', args.Class,
    '-ee', str(args.maxEE),
    '-tR', str(args.trimRight),
    '-mL', str(args.minLen),
    '-tQ', str(args.truncQ),
    '-mC', str(args.max_consist),
    '-wA', str(args.omegaA),
    '-jC', str(args.justConcatenate),
#    '--saveRdata', args.saveRdata,
    '--bimera'
    ]
    dada2_cmd_str = ' '.join(map(str, dada2_cmd))
    run_cmd(dada2_cmd_str)
    print('DADA2 step complete!')

    return

if __name__ == "__main__":
    main()
