# Snakefile: Amplicon_Pipeline.smk
#
# (C) 2023 Mariana Barnes, Ashley Osborne, Hidayat Trimarsanto
#
#this is the snakefile for a dada2 based analysis of microhap fastq data. This is based on the broads institute's malaria amplicon pipeline: https://github.com/broadinstitute/malaria-amplicon-pipeline/blob/main/inputs.json


import os
import pathlib
from ngs_pipeline import cerr, fileutils

microhaps_basedir = os.environ['MICROHAPS_BASEDIR']

fasta = microhaps_basedir + '/microhaps_pipeline/' + config['fasta']
primer_fw = microhaps_basedir + '/microhaps_pipeline/' + config['primer_fw']
primer_rev = microhaps_basedir + '/microhaps_pipeline/' + config['primer_rev']


# define all output files 

out_dir = config['outdir']
#in_dir = config['indir']
in_dir=''

read_files = fileutils.ReadFileDict(config['infiles'], config['underscore'])
IDs = read_files.keys()

rule all:
    input:
        f"{out_dir}/malamp/meta",
        f"{out_dir}/malamp/dada2/seqtab.tsv",
        #"ASVTable.txt",
        #"outputCIGAR.tsv",
        expand(f'{out_dir}/trimmed/{{sample}}-0_R1.trimmed.fastq.gz', sample=IDs)
 

#create input file list
#rule create_input:
#    localrule: True
#    input:
#        f"{in_dir}",
#    output:
#        f"{out_dir}/mal_amp/samples_file.csv",
#    run:
#        import pathlib
#        import pandas as pd
#
#        filenames = [fn.name.split('_')[0] for fn in pathlib.Path(input[0]).glob("*_R1.fastq.gz")]
#        pd.DataFrame({'sample': filenames}).to_csv(output[0], index=False)


#trim reads using trimmomatic - adapters and primers 

rule trim:
    input:
        R1 = f"{in_dir}/{{sample}}_R1.fastq.gz",
        R2 = f"{in_dir}/{{sample}}_R2.fastq.gz"
    output:
        R1 = f"{out_dir}/trimmed/{{sample}}_R1.trimmed.fastq.gza",
        R2 = f"{out_dir}/trimmed/{{sample}}_R2.trimmed.fastq.gza"
    params:
        platform = config['platform'],
    shell: 
        "trimmomatic {params.platform} {input.R1} {input.R2} {output.R1} {output.R2} LEADING:3 TRAILING:3 SLIDINGWINDOW:4:%(trim_qv)s MINLEN:20 2> %(sample)s.trimlog"

rule trim_1:
    localrule: True
    input:
        lambda w: read_files.get_read_file(w),
        prim_fw = primer_fw,
        prim_rv = primer_rev
    output:
        #f"{outdir}/{{sample}}/reads/raw-{{idx}}_R1.fastq.gz",
        #f"{outdir}/{{sample}}/reads/raw-{{idx}}_R2.fastq.gz"
        R1 = f"{out_dir}/trimmed/{{sample}}-{{idx}}_R1.trimmed.fastq.gz",
        R2 = f"{out_dir}/trimmed/{{sample}}-{{idx}}_R2.trimmed.fastq.gz",
    params:
        platform = config['platform'],
        trim_qv = config['trimqv']
    shell: 
        "cutadapt -g file:{input.prim_fw} -G file:{input.prim_rv} -o {output.R1} -p {output.R2} --pair-adapters --discard-untrimmed --action=trim {input[0]} {input[1]}"

#run inividual preprocessing for dada2 in R TO DO

rule create_meta:
    localrule: True
    input:
        expand(f'{out_dir}/trimmed/{{sample}}-0_R1.trimmed.fastq.gz', sample=IDs)
    output:
        f"{out_dir}/malamp/meta"
    log:
        f"{out_dir}/mylog.txt"
    shell: 
        """
        mkdir -p {out_dir}/malamp
        python {microhaps_basedir}/scripts/create_meta.py --path_to_fq {out_dir}/trimmed/ --output_file {output} --pattern_fw *R1.trimmed.fastq.gz --pattern_rv *R2.trimmed.fastq.gz 1> {log} 2> {log}
        """


rule run_dada2R:
    localrule: True
    input:
        meta = f"{out_dir}/malamp/meta",
    output:
        f"{out_dir}/malamp/dada2/seqtab.tsv"
    params:
        dir = f"{out_dir}/malamp/dada2",
        output_filename = "seqtab.tsv",
        maxEE = config['maxEE'],
        trimRight = config['trim_right'],
        minLen = config['min_length'],
        truncQ = config['truncQ'],
        max_consist = config['max_consist'],
        omega_a = config['omegaA'],
        justConcatenate = config['justconcat'],
        _class = config['class']

    shell:
        """
        Rscript {microhaps_basedir}/scripts/runDADA2.R  \
            --path_to_meta {input.meta} \
            --class {params._class} \
            --dir {params.dir} \
            --output_filename {params.output_filename} \
            --maxEE {params.maxEE} \
            --trimRight {params.trimRight} \
            --minLen {params.minLen} \
            --truncQ {params.truncQ} \
            --max_consist {params.max_consist} \
            --omega_a {params.omega_a} \
            --justConcatenate {params.justConcatenate} \
        """