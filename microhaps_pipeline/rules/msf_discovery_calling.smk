

import ngs_pipeline.rules

include: ngs_pipeline.rules.path("multistep_variant_calling.smk")

if config.get("joint_discovery"):
    final_output = f"{outdir}/joint/concatenated.vcf.gz.tbi"
else:
    final_output = []


ruleorder: mark_prepared > run_prepare_sample_directory
ruleorder: gather_stats_sample > gather_stats
ruleorder: prepare_shallow_dir > index_bai

rule mark_prepared:
    localrule: True
    input:
        expand(f"{outdir}/analysis/{{sample}}/maps/mapped-final.bam", sample=IDs),
        expand(f"{outdir}/analysis/{{sample}}/maps/mapped-final.bam.bai", sample=IDs),
        expand(f"{outdir}/analysis/{{sample}}/reads", sample=IDs),
    output:
        f"{outdir}/analysis/._prepared_",
    shell:
        """
        touch {output}
        """

rule prepare_shallow_dir:
    localrule: True
    input:
        f"{outdir}/samples/{{sample}}/maps/final.bam",
        f"{outdir}/samples/{{sample}}/maps/final.bam.bai",
    output:
        f"{outdir}/analysis/{{sample}}/maps/mapped-final.bam",
        f"{outdir}/analysis/{{sample}}/maps/mapped-final.bam.bai",
        directory(f"{outdir}/analysis/{{sample}}/reads"),
    params:
        outdir = f"{outdir}/analysis/{{sample}}/maps",
    shell:
        """
        mkdir -p {params.outdir};
        ln -s {input[0]} {output[0]};
        ln -s {input[1]} {output[1]};
        mkdir -p {output[2]};
        """

rule mark_discovery_done:
    localrule: True
    input:
        final_output
    output:
        f"{outdir}/.__discovery__"
    shell:
        """
        touch {output}
        """

