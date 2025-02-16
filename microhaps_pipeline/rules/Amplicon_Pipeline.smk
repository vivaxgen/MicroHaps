# Snakefile: Amplicon_Pipeline.smk
#
# (C) 2023 Mariana Barnes, Ashley Osborne, Hidayat Trimarsanto
# (C) 2024 Ludwig Hoon 
#this is the snakefile for a dada2 based analysis of microhap fastq data. This is based on the broads institute's malaria amplicon pipeline: https://github.com/broadinstitute/malaria-amplicon-pipeline/blob/main/inputs.json


import os
import pathlib
from ngs_pipeline import cerr, fileutils

include: "params_mhap.smk"

rule all:
    input:
        f"{outdir}/malamp/meta",
        f"{outdir}/malamp/dada2/seqtab.tsv",
        #"ASVTable.txt",
        #"outputCIGAR.tsv",
        expand(f'{outdir}/trimmed/{{sample}}-0_R1.trimmed.fastq.gz', sample=IDs),
        f"{outdir}/malamp/ASVTable.txt",
        f"{outdir}/malamp/ASVSeqs.fasta",
        f"{outdir}/malamp/outputCIGAR.tsv"
 

#create input file list
#rule create_input:
#    localrule: True
#    input:
#        f"{indir}",
#    output:
#        f"{outdir}/malamp/samples_file.csv",
#    run:
#        import pathlib
#        import pandas as pd
#
#        filenames = [fn.name.split('_')[0] for fn in pathlib.Path(input[0]).glob("*_R1.fastq.gz")]
#        pd.DataFrame({'sample': filenames}).to_csv(output[0], index=False)


#trim reads using trimmomatic - adapters and primers 

rule trim:
    input:
        R1 = f"{indir}/{{sample}}_R1.fastq.gz",
        R2 = f"{indir}/{{sample}}_R2.fastq.gz"
    output:
        R1 = f"{outdir}/trimmed/{{sample}}_R1.trimmed.fastq.gza",
        R2 = f"{outdir}/trimmed/{{sample}}_R2.trimmed.fastq.gza"
    params:
        platform = config['platform'],
    shell: 
        "trimmomatic {params.platform} {input.R1} {input.R2} {output.R1} {output.R2} LEADING:3 TRAILING:3 SLIDINGWINDOW:4:%(trim_qv)s MINLEN:20 2> %(sample)s.trimlog"

rule trim_1:
    input:
        unpack(read_files.get_read_file_as_dict),
        prim_fw = primer_fw_file,
        prim_rv = primer_rev_file,
    output:
        #f"{outdir}/{{sample}}/reads/raw-{{idx}}_R1.fastq.gz",
        #f"{outdir}/{{sample}}/reads/raw-{{idx}}_R2.fastq.gz"
        R1 = f"{outdir}/trimmed/{{sample}}-{{idx}}_R1.trimmed.fastq.gz",
        R2 = f"{outdir}/trimmed/{{sample}}-{{idx}}_R2.trimmed.fastq.gz",
    params:
        platform = config['platform'],
        trim_qv = config['trimqv'],
        additional_params = f"--nextseq-trim={config.get('trimqv', 20)}" if is_nextseq_or_novaseq() else ""
    shell: 
        "cutadapt -g file:{input.prim_fw} -G file:{input.prim_rv} -o {output.R1} -p {output.R2} --pair-adapters --discard-untrimmed {params.additional_params} --action=trim {input.read1} {input.read2}"

#run inividual preprocessing for dada2 in R TO DO

rule create_meta:
    localrule: True
    input:
        expand(f'{outdir}/trimmed/{{sample}}-0_R1.trimmed.fastq.gz', sample=IDs)
    output:
        f"{outdir}/malamp/meta"
    log:
        f"{outdir}/mylog.txt"
    shell: 
        """
        mkdir -p {outdir}/malamp
        python {microhaps_basedir}/scripts/create_meta.py --path_to_fq {outdir}/trimmed/ --output_file {output} --pattern_fw *-0_R1.trimmed.fastq.gz --pattern_rv *-0_R2.trimmed.fastq.gz 1> {log} 2> {log}
        """


rule run_dada2R:
    threads: 16
    input:
        meta = f"{outdir}/malamp/meta",
    output:
        f"{outdir}/malamp/dada2/seqtab.tsv"
    params:
        dir = f"{outdir}/malamp/dada2",
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


rule post_process:
    input:
        f"{outdir}/malamp/dada2/seqtab.tsv"
    output:
        Table = f"{outdir}/malamp/ASVTable.txt",
        Seqs = f"{outdir}/malamp/ASVSeqs.fasta"
    shell:
        """
        Rscript {microhaps_basedir}/scripts/postProc_dada2.R \
            -s {input} \
            --strain PvP01 \
            -ref {insertseq_file} \
            -o {output.Table} \
            --fasta \
            --parallel \
        """


rule asv_to_cigar:
    input:
        Table = f"{outdir}/malamp/ASVTable.txt",
        Seqs = f"{outdir}/malamp/ASVSeqs.fasta",
        seqtab = f"{outdir}/malamp/dada2/seqtab.tsv"
    output:
        cigar = f"{outdir}/malamp/outputCIGAR.tsv",
        asv_to = f"{outdir}/malamp/asv_to_cigar"
    shell:
        """
        python {microhaps_basedir}/scripts/ASV_to_CIGAR.py {input.Seqs} {input.Table} {input.seqtab} {output.cigar} --asv_to_cigar {output.asv_to} \
        -a {outdir}/alignments  --amp_db {insertseq_file} \
        """

# EOF

