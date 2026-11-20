from ngs_pipeline import cerr
from ngs_pipeline.rules import pkg

include: pkg("ngs_pipeline::msf/init.smk")
include: pkg("ngs_pipeline::params_region.smk")

include: pkg("params_mhap.smk")

include: pkg("ngs_pipeline::trimmer/null_lr.smk")


use rule reads_trimming_lr as reads_trimming_lr_no_model with:
    input:
        read = "{anypath}reads/raw-{idx}_R0.fastq.gz",

ruleorder: reads_trimming_lr_no_model > reads_trimming_lr

# 1 - map to the reference sequence
include: pkg(config.get("mapper_wf", "ngs_pipeline::mapper/minimap2_lr.smk"))

# 2 - merge multiple bams to a single bam (maps/mapped-final.bam)
include: pkg("ngs_pipeline::helper/map_handler.smk")

# 3 - include all statistics utitlities
include: pkg("ngs_pipeline::helper/stats.smk")

# depth & coverage per amplicon
include: pkg("msf_final_bam_to_fastq.smk")

# 4 - drug resistance report
include: pkg("msf_drug_resistance.smk")


depths_mapped_plot = [f"{outdir}/depths-mapped.png"]
if len(Plates) > 1:
    depths_mapped_plot += [f"{outdir}/depths-mapped_{plate}.png" for plate in Plates]


rule all_drug_resistance_lr:
    input:
        f"{outdir}/stats.tsv",
        *depths_mapped_plot,
        f"{outdir}/coverages-mapped.tsv",
        f"{outdir}/depths-mapped.tsv",
        f"{outdir}/malamp/drug_resistance.tsv",
        f"{outdir}/malamp/drug_resistance_flagged.tsv",
        expand(f"{outdir}/samples/{{sample}}/drugs/drug_resistance_stats.tsv", sample=IDs),

rule mapping_stats:
    localrule: True
    input:
        expand(f"{outdir}/samples/{{sample}}/logs/stats.tsv", sample=read_files.samples()),
    output:
        f"{outdir}/stats.tsv"
    shell:
        'ngs-pl gather-stats -o {output} {outdir}/samples'

# EOF
