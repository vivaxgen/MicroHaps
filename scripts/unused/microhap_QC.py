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
  
    for sample in samples:
        args.sample = sample
        run_cmd("fastqc -t 6 %(sample)s_R1.fastq.gz %(sample)s_R2.fastq.gz -o FASTQC_results" % vars(args))
        run_cmd("bwa mem -t 6 -R \"@RG\\tID:M00859\\tSM:%(sample)s\\tLB:MicroHap\\tPU:L6WVN:1\\tPL:Illumina\" %(ref)s %(sample)s_R1.fastq.gz %(sample)s_R2.fastq.gz | samclip --ref %(ref)s --max 50 | samtools sort -o %(sample)s.bam -" % vars(args))
        run_cmd("samtools index %(sample)s.bam" % vars(args))
        run_cmd("samtools flagstat %(sample)s.bam > %(sample)s.flagstat.txt" % vars(args))
        #run_cmd("mosdepth -x -b %(bed)s %(sample)s --thresholds 1,10,20,30  %(sample)s.bam" % vars(args))
        run_cmd("bedtools coverage -a %(bed)s -b %(sample)s.bam -mean > %(sample)s_region_coverage.txt" % vars(args))
        run_cmd("sambamba depth base %(sample)s.bam > %(sample)s.coverage.txt" % vars(args))

    run_cmd("multiqc FASTQC_results")

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
    
#    with open("bam_list.txt","w") as O:
#        for s in samples:
#            O.write("%s.bam\n" % (s))

#    run_cmd("freebayes -f %(ref)s -t %(bed)s -L bam_list.txt --haplotype-length -1 --min-coverage 50 --min-base-quality %(min_base_qual)s --gvcf --gvcf-dont-use-chunk true | bcftools view -T %(bed)s | bcftools norm -f %(ref)s | bcftools sort -Oz -o combined.genotyped.vcf.gz" % vars(args))

# Set up the parser
parser = argparse.ArgumentParser(description='MicroHaplotype Quality Control script',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--index-file',type=str,help='CSV file containing field "Sample"',required=True)
parser.add_argument('--ref',type=str,help='Reference fasta',required=True)
#parser.add_argument('--gff',type=str,help='GFF file',required=True)
parser.add_argument('--bed',type=str,help='BED file with MicroHaplotype locations',required=True)
#parser.add_argument('--min-base-qual',default=30,type=int,help='Minimum base quality to use by freebayes')
#parser.add_argument('--min-adf',type=float,help='Set a minimum frequency for a mixed call')
#parser.add_argument('--min-variant-qual',default=30,type=int,help='Quality value to use in the sliding window analysis')
#parser.add_argument('--min-sample-af',default=0.05,type=float,help='Quality value to use in the sliding window analysis')
parser.add_argument('--version', action='version', version='%(prog)s 1.0')
parser.set_defaults(func=main)

args = parser.parse_args()
args.func(args)
