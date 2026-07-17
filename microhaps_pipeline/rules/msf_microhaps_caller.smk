__author__ = "Hidayat Trimarsanto"
__copyright__ = "(C) 2024-2026 Hidayat Trimarsanto"
__email__ = "trimarsanto@gmail.com,hidayat.trimarsanto@menzies.edu.au"
__license__ = "MIT"

from ngs_pipeline import cerr
from ngs_pipeline.rules import pkg

#include: inc("ngs_pipeline::general_params.smk")
#include: inc("ngs_pipeline::params_region.smk")
# include: inc("ngs_pipeline::helper/utilities.smk")

include: pkg("ngs_pipeline::msf/init.smk")
include: pkg("ngs_pipeline::params_region.smk")

include: pkg("params_mhap.smk")


# 0 - link read files to individual sample directories
#include: inc("ngs_pipeline::msf_prepare_sample_files.smk")

include: pkg("ngs_pipeline::trimmer/null.smk")

# 1 - map to the reference sequence
#     need to use bwa-mem2 since minimap2 (and mm2plus) does not set flags for
#     READ1 and READ2 in the output file, which is necessary for samtools fastq
#     to generate paired-end fastq files
include: pkg("ngs_pipeline::mapper/bwa-mem2.smk")

# 2 - merge multiple bams to a single bam and then generate new bam for reads 
#     that map to the specific regions of the reference sequence
#     the bam file {sample}-{idx}.bam is suitable for uploading to public databases
#include: inc("ngs_pipeline::legacy/msf_merge_map.smk")
include: pkg("ngs_pipeline::helper/map_handler.smk")

# 3 - include all statistics utitlities
#include: inc("ngs_pipeline::legacy/msf_stats.smk")
include: pkg("ngs_pipeline::helper/stats.smk")


include: pkg("msf_final_bam_to_fastq.smk")


# 3 - perform microhaps calling

# 4 - hard-trim primers from each of merged reads
include: pkg("msf_trim_dedup.smk")

# 5 - generate a FASTQ file from the trimmmed reads

# 6 - call dada2 for denoising FASTQ reads
if config.get("merge_map") == "dada2":
    include: pkg("msf_merge_denoise_dada2.smk")
elif config.get("merge_map") == "fastp":
    include: pkg("msf_fastp_merge.smk")
elif config.get("merge_map") == "fastp_dada2":
    include: pkg("msf_fastp_merge.smk")
    include: pkg("msf_dada2_denoise.smk")
else:
    raise ValueError(f"Unknown merge_map option: {config.get('merge_map')}")


# 7 - generate Haplotype table
include: pkg("msf_post_process_merged.smk")

# 8 - qc Haplotypes
include: pkg("msf_qc_haplotype.smk")

include: pkg("msf_discovery_calling.smk")

new_postprocess = config.get('post_process', "old")

merging_output = f"{outdir}/malamp/{config.get('merge_map')}/seqtab.tsv"

drug_resistance_output = []
if config.get("drugs_resistance_aa_pos", None) is not None:
    include: pkg("msf_drug_resistance.smk")
    drug_resistance_output.append(f"{outdir}/malamp/drug_resistance.tsv")
    drug_resistance_output.append(f"{outdir}/malamp/drug_resistance_flagged.tsv")

presence_absence_output = []
if config.get("presence_absence_markers", None) is not None:
    include: pkg("msf_presence_absence.smk")
    presence_absence_output.append(f"{outdir}/malamp/presence_absence.tsv")


rule all_microhaps:
    input:
        f"{outdir}/stats.tsv",
        f"{outdir}/depths-mapped.png",
        f"{outdir}/coverages-mapped.tsv",
        f"{outdir}/.__discovery__",
        merging_output,
        *([f"{outdir}/malamp/outputHaplotypes.tsv", f"{outdir}/malamp/outputHaplotypes_rm_ins.tsv"] if new_postprocess != "old" else [f"{outdir}/malamp/outputCIGAR.tsv"]),
        f"{outdir}/malamp/depths-microhaps.png",
        f"{outdir}/malamp/depth-ratio-markers.png",
        *drug_resistance_output,
        *presence_absence_output,

rule seqtab:
    input:
        f"{outdir}/malamp/dada2/seqtab.tsv",

rule asvtable:
    input:
        f"{outdir}/malamp/ASVTable.txt",

# EOF
