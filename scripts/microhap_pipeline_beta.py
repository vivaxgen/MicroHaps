#! /usr/bin/env python
import sys
import argparse
import subprocess as sp
import csv
import fastq2matrix as fm
from fastq2matrix import run_cmd
from collections import defaultdict
import gzip
import os
import shutil
import threading
import multiprocessing
import glob

def run_cmd(cmd):
    sys.stderr.write("Running command:\n%s\n\n" % cmd)
    with open("/dev/null","w") as O:
        res = sp.call(cmd,shell=True,stderr=O,stdout=O)
    if res!=0:
        sys.exit("Error running last command, please check!\n")
      

def main(args):

    samples = []
    reader = csv.DictReader(open(args.index_file))
    if "sample" not in reader.fieldnames:
        reader = csv.DictReader(open(args.index_file,encoding='utf-8-sig'))
    for row in reader:
        if row["sample"]=="": continue
        samples.append(row["sample"])

    fm.bwa_index(args.ref)
    fm.create_seq_dict(args.ref)
    fm.faidx(args.ref)
    
    run_cmd("mkdir FASTQC_results")
    run_cmd("mkdir bam_files")
    run_cmd("mkdir cov_stats")
    run_cmd("mkdir untrimmed_fastq")
    run_cmd("mkdir trimmed_fastq")
  
    for sample in samples:
        args.sample = sample
        run_cmd("fastqc -t 6 %(sample)s_R1.fastq.gz %(sample)s_R2.fastq.gz -o FASTQC_results" % vars(args))

        if args.trim:
            run_cmd("trimmomatic PE %(sample)s_R1.fastq.gz %(sample)s_R2.fastq.gz %(sample)s_R1.trimmed.fastq.gz %(sample)s_R1.untrimmed.fastq.gz %(sample)s_R2.trimmed.fastq.gz %(sample)s_R2.untrimmed.fastq.gz LEADING:3 TRAILING:3 SLIDINGWINDOW:4:%(trim_qv)s MINLEN:20 2> %(sample)s.trimlog" % vars(args))
            run_cmd("bwa mem -t 6 -R \"@RG\\tID:M00859\\tSM:%(sample)s\\tLB:MicroHap\\tPU:L6WVN:1\\tPL:Illumina\" %(ref)s %(sample)s_R1.trimmed.fastq.gz %(sample)s_R2.trimmed.fastq.gz | samclip --ref %(ref)s --max 50 | samtools sort -o %(sample)s.bam -" % vars(args))
        else:
            run_cmd("bwa mem -t 6 -R \"@RG\\tID:M00859\\tSM:%(sample)s\\tLB:MicroHap\\tPU:L6WVN:1\\tPL:Illumina\" %(ref)s %(sample)s_R1.fastq.gz %(sample)s_R2.fastq.gz | samclip --ref %(ref)s --max 50 | samtools sort -o %(sample)s.bam -" % vars(args))

        run_cmd("samtools index %(sample)s.bam" % vars(args))
        run_cmd("samtools flagstat %(sample)s.bam > %(sample)s.quickstats.txt" % vars(args))
        run_cmd("bedtools coverage -a %(bed)s -b %(sample)s.bam -mean > %(sample)s_region_coverage.txt" % vars(args))
        run_cmd("sambamba depth base %(sample)s.bam > %(sample)s.position_coverage.txt" % vars(args))

    run_cmd("multiqc FASTQC_results")
    
    destination_directory = 'trimmed_fastq'
    os.makedirs(destination_directory, exist_ok=True)
    for filename in os.listdir(os.getcwd()):
        if filename.endswith(".trimmed.fastq.gz") or filename.endswith(".trimlog"):
            source_path = os.path.join(os.getcwd(), filename)
            destination_path = os.path.join(destination_directory, filename)

            # Move the file
            shutil.move(source_path, destination_path)
    print("Trimmed FASTQ files moved successfully.")

    run_cmd('create_meta.py --path_to_fq trimmed_fastq --output_file %(output_file)s --pattern_fw "%(pattern_fw)s" --pattern_rv "%(pattern_rv)s"' % vars(args))
    
    destination_directory = 'untrimmed_fastq'
    os.makedirs(destination_directory, exist_ok=True)
    for filename in os.listdir(os.getcwd()):
        if filename.endswith(".untrimmed.fastq.gz"):
            source_path = os.path.join(os.getcwd(), filename)
            destination_path = os.path.join(destination_directory, filename)

            # Move the file
            shutil.move(source_path, destination_path)
    print("Untrimmed FASTQ files moved successfully.")
    
    destination_directory = 'bam_files'
    os.makedirs(destination_directory, exist_ok=True)
    for filename in os.listdir(os.getcwd()):
        if filename.endswith(".bam") or filename.endswith(".bam.bai"):
            source_path = os.path.join(os.getcwd(), filename)
            destination_path = os.path.join(destination_directory, filename)

            # Move the file
            shutil.move(source_path, destination_path)
    print("Bam and index files moved successfully.")

    destination_directory = 'cov_stats'
    os.makedirs(destination_directory, exist_ok=True)
    for filename in os.listdir(os.getcwd()):
        if filename.endswith(".txt"):
            source_path = os.path.join(os.getcwd(), filename)
            destination_path = os.path.join(destination_directory, filename)

            # Move the file
            shutil.move(source_path, destination_path)
    print("Coverage stats moved successfully.")

    # Run Amplicon Pipeline and DADA 2
    run_cmd('AmpliconPipeline.py --path_to_meta %(output_file)s --pr1 %(pr1)s --pr2 %(pr2)s --Class %(Class)s --maxEE %(maxEE)s --trimRight %(trimRight)s --minLen %(minLen)s --truncQ %(truncQ)s --max_consist %(max_consist)s --omegaA %(omegaA)s --justConcatenate %(justConcatenate)s --saveRdata %(saveRdata)s' % vars(args))

# Set up the parser
parser = argparse.ArgumentParser(description='MicroHaplotype Pipeline',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--index-file',type=str,help='CSV file containing field "Sample"',required=True)
parser.add_argument('--ref',type=str,help='Reference fasta',required=True)
parser.add_argument('--bed',type=str,help='BED file with MicroHaplotype locations',required=True)
parser.add_argument('--trim',action="store_true",help='Perform triming')
parser.add_argument('--trim-qv',default=5,type=int,help='Quality value to use in the sliding window analysis')
parser.add_argument('--output_file', type=str, help='Output meta file; to be used in path to meta', required=True)
parser.add_argument('--pattern_fw', type=str, help='Pattern for forward reads, e.g. "*_R1.fastq.gz"', required=True)
parser.add_argument('--pattern_rv', type=str, help='Pattern for reverse reads, e.g. "*_R2.fastq.gz"', required=True)
parser.add_argument('--path_to_meta', help="Path to input fastq files", required=True)
#parser.add_argument('--keep_primers', action="store_true",default=1, help="Skip primer removal step")
parser.add_argument('--pr1', help="Path to forward primers FASTA file", required=True)
parser.add_argument('--pr2', help="Path to reverse primers FASTA file", required=True)
parser.add_argument('--Class', default="parasite", help="Specify Analysis class. Accepts one of two: parasite/vector")
parser.add_argument('--maxEE', default="5,5", help="Maximum Expected errors (dada2 filtering argument)")
parser.add_argument('--trimRight', default="10,10", help="Hard trim number of bases at 5` end (dada2 filtering argument)")
parser.add_argument('--minLen', default=30, help="Minimum length filter (dada2 filtering argument)")
parser.add_argument('--truncQ', default="5,5", help="Soft trim bases based on quality (dada2 filtering argument)")
parser.add_argument('--max_consist', default=10, help="Number of cycles for consistency in error model (dada2 argument)")
parser.add_argument('--omegaA', default=1e-120, help="p-value for the partitioning algorithm (dada2 argument)")
parser.add_argument('--justConcatenate', default=0, help="whether reads should be concatenated with N's during merge (dada2 argument)")
parser.add_argument('--saveRdata',default="", help="Optionally save dada2 part of this run as Rdata object")
parser.add_argument('--version', action='version', version='%(prog)s 1.0')
parser.set_defaults(func=main)

args = parser.parse_args()
args.func(args)
