__author__ = "Hidayat Trimarsanto"
__copyright__ = "(C) 2024, Hidayat Trimarsanto"
__email__ = "trimarsanto@gmail.com,hidayat.trimarsanto@menzies.edu.au"
__license__ = "MIT"

import ngs_pipeline.rules

include: ngs_pipeline.rules.path("general_params.smk")
include: ngs_pipeline.rules.path("params_region.smk")
include: ngs_pipeline.rules.path("utilities.smk")

include: "params_mhap.smk"


# 0 - link read files to individual sample directories
include: ngs_pipeline.rules.path("msf_prepare_sample_files.smk")

include: ngs_pipeline.rules.path("msf_trimmer_null.smk")

# 1 - map to the reference sequence
#     need to use bwa-mem2 since minimap2 (and mm2plus) does not set flags for
#     READ1 and READ2 in the output file, which is necessary for samtools fastq
#     to generate paired-end fastq files
include: ngs_pipeline.rules.path("msf_mapper_bwa-mem2.smk")

# 2 - merge multiple bams to a single bam and then generate new bam for reads 
#     that map to the specific regions of the reference sequence
#     the bam file {sample}-{idx}.bam is suitable for uploading to public databases
include: ngs_pipeline.rules.path("msf_merge_map.smk")

# 3 - include all statistics utitlities
include: ngs_pipeline.rules.path("msf_stats.smk")

include: "msf_final_bam_to_fastq.smk"


# 3 - perform microhaps calling

include: "msf_trim_merge_denoise_dada2.smk"

# 4 - hard-trim primers from each of merged reads

# 5 - generate a FASTQ file from the trimmmed reads

# 6 - call dada2 for denoising FASTQ reads

# 7 - generate ASV table

rule all:
    input:
        f"{outdir}/stats.tsv",
        f"{outdir}/final.depths.tsv",
        f"{outdir}/final.coverages.tsv",
        f"{outdir}/malamp/dada2/seqtab.tsv",
        f"{outdir}/malamp/outputCIGAR.tsv",
        f"{outdir}/malamp/depths.tsv"

rule seqtab:
    input:
        f"{outdir}/malamp/dada2/seqtab.tsv",

rule asvtable:
    input:
        f"{outdir}/malamp/ASVTable.txt",

# EOF
