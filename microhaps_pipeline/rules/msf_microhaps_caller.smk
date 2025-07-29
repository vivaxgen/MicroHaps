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

include: ngs_pipeline.rules.path("utilities.smk")

include: "msf_final_bam_to_fastq.smk"


# 3 - perform microhaps calling

# 4 - hard-trim primers from each of merged reads
include: "msf_trim.smk"

# 5 - generate a FASTQ file from the trimmmed reads

# 6 - call dada2 for denoising FASTQ reads
if config.get("merge_map") == "dada2":
    include: "msf_merge_denoise_dada2.smk"
elif config.get("merge_map") == "bbmap_merge" or config.get("merge_map") == "bbmerge":
    include: "msf_bbmap_bbmerge_vsearch.smk"
else:
    raise ValueError(f"Unknown merge_map option: {config.get('merge_map')}")


# 7 - generate Haplotype table
include: "msf_post_process_merged.smk"

# 8 - qc Haplotypes
include: "msf_qc_haplotype.smk"

include: "msf_discovery_calling.smk"

new_postprocess = config.get('post_process', "old")

dada2_output = f"{outdir}/malamp/dada2/seqtab.tsv" # default
bbmap_merge_output = f"{outdir}/malamp/bbmap_merge/seqtab.tsv"
bbmerge_output = f"{outdir}/malamp/bbmerge/seqtab.tsv"

match config.get("merge_map", "dada2"):
    case "dada2": merging_output = dada2_output
    case "bbmap_merge": merging_output = bbmap_merge_output
    case "bbmerge": merging_output = bbmerge_output

rule all_microhaps:
    input:
        f"{outdir}/stats.tsv",
        f"{outdir}/depths-mapped.png",
        f"{outdir}/coverages-mapped.tsv",
        f"{outdir}/.__discovery__",
        merging_output,
        f"{outdir}/malamp/outputHaplotypes.tsv" if new_postprocess != "old" else f"{outdir}/malamp/outputCIGAR.tsv",
        f"{outdir}/malamp/depths-microhaps.png",
        f"{outdir}/malamp/depth-ratio-markers.png",
        f"{outdir}/malamp/{bbmap_or_bbmerge}/dereplicated_counted.fasta"

rule seqtab:
    input:
        f"{outdir}/malamp/dada2/seqtab.tsv",

rule asvtable:
    input:
        f"{outdir}/malamp/ASVTable.txt",

# EOF
