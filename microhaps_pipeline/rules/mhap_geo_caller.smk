__author__ = "Hidayat Trimarsanto"
__copyright__ = "(C) 2024, Hidayat Trimarsanto"
__email__ = "trimarsanto@gmail.com,hidayat.trimarsanto@menzies.edu.au"
__license__ = "MIT"

import ngs_pipeline.rules

include: ngs_pipeline.rules.path("general_params.smk")
include: ngs_pipeline.rules.path("params_region.smk")
# include: ngs_pipeline.rules.path("utilities.smk")

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


rule collate_bam:
    group: "{sample}_hap_step1"
    resources:
        runtime="2h"
    threads: 4
    input:
        bam = f"{outdir}/samples/{{sample}}/maps/final.bam",
        refseq = refseq,
    output:
        collated_bam = f"{outdir}/samples/{{sample}}/maps/final.collate.bam"
    shell:
        """samtools collate -O {input.bam} \
           | samtools calmd -u - {input.refseq} \
           | samtools view -o {output}
        """

rule generate_haplotype_table:
    group: "{sample}_hap_step1"
    threads: 4
    resources:
        runtime="9h"
    input:
        bam = f"{outdir}/samples/{{sample}}/maps/final.collate.bam",
        variants_list = get_abspath(config['geo_variant_list'], microhaps_basedir),
    output:
        haplotype = f"{outdir}/samples/{{sample}}/geo/partial_haps.tsv",
        snps = f"{outdir}/samples/{{sample}}/geo/snps_table.tsv",
    params:
        min_qual = 20,
        min_mapq = 30,
        GOOD_QUAL = 35,
        BAD_QUAL = 20,
    shell:
        """
        ngs-pl construct-pseudo-haplotypes --min_qual {params.min_qual} --min_mapq {params.min_mapq} \
        --GOOD_QUAL {params.GOOD_QUAL} --BAD_QUAL {params.BAD_QUAL} --variants_list {input.variants_list} \
        --pileup --pileout {output.snps} --sample {wildcards.sample} -o {output.haplotype} -j {threads} \
        {input.bam}
        """

rule filter_partial_pseudohaplotype:
    group: "{sample}_hap_step2"
    resources:
        runtime="1h"
    input:
        f"{outdir}/samples/{{sample}}/geo/partial_haps.tsv"
    output:
        f"{outdir}/samples/{{sample}}/geo/filtered_partial_haps.tsv"
    run:
        import pandas as pd
        df = pd.read_table(input[0], dtype={"marker": str})
        df["n_known"] = df["haplotype"].apply(lambda x: sum([1 for a in x if a in "ACGT"]))
        df["completeness"] = df["n_known"] / df["haplotype"].apply(len)
        filtered_df = df.query("completeness >= 0.25 and count > 5").copy() # and n_known >= 2 
        filtered_df.to_csv(output[0], sep="\t", index=False)

rule assemble_partial_pseudohaplotye:
    group: "{sample}_hap_step2"
    resources:
        runtime="2h"
    input:
        f"{outdir}/samples/{{sample}}/geo/filtered_partial_haps.tsv"
    output:
        f"{outdir}/samples/{{sample}}/geo/assembled_{{sample}}.tsv"
    params:
        max_consider = 15,
        prioritise = "completeness", # or "count"
        strict = "--strict"
    shell:
        "ngs-pl assemble-partial-haplotype -o {output} {params.strict} -m {params.max_consider} -p {params.prioritise} {input}"

min_depth_per_base_MHAPS = 25
haplotype_pass_min_count_MHAPS = 25
haplotype_pass_min_ratio_MHAPS = 0.01

error_rate = 0.1
confidence = 0.99

rule filter_assembled_pseudohaplotype:
    threads: 16
    resources:
        runtime="1h"
    input:
        files = expand(f"{outdir}/samples/{{sample}}/geo/assembled_{{sample}}.tsv", sample=IDs),
        marker_df = get_abspath(config['geo_variant_list'], microhaps_basedir),
    output:
        long_df = f"{outdir}/geo/long.tsv",
        stats = f"{outdir}/geo/stats.tsv",
        wide_hap_df = f"{outdir}/geo/wide_haplotypes.tsv",
        marker_stats = f"{outdir}/geo/marker_stats.tsv",
    params:
        min_depth_per_base = min_depth_per_base_MHAPS,
        haplotype_pass_min_count= haplotype_pass_min_count_MHAPS,
        haplotype_pass_min_ratio = haplotype_pass_min_ratio_MHAPS,
        use_poisson_threshold=True,
        error_rate=error_rate,
        confidence=confidence,
    run:
        from ngs_pipeline.pseudohap_utils import process_individual_sample_haplotype, long_to_wide_ordered, wide_sample_to_wide_haplotype, marker_missingness
        import pandas as pd
        from concurrent.futures import ProcessPoolExecutor, as_completed
        from tqdm import tqdm
        from os.path import basename
        marker_df = pd.read_table(input.marker_df, header=None, names=["chr", "pos0", "pos", "marker"], dtype={"marker": str})
        append_list = []
        stat_list = []
        with ProcessPoolExecutor(max_workers=threads) as executor:
            futures = {}
            for f in tqdm(input.files, total=len(input.files), desc="Submitting samples for processing"):
                df = pd.read_table(f, dtype={"marker": str})
                sample_id = basename(f).replace(".tsv", "").replace("assembled_", "")

                futures[executor.submit(
                    process_individual_sample_haplotype, df, marker_df, params.min_depth_per_base, params.haplotype_pass_min_count,
                    params.haplotype_pass_min_ratio, sample_id, params.use_poisson_threshold, params.error_rate, params.confidence)
                ] = sample_id

            for future in tqdm(as_completed(futures), total=len(input.files), desc="Processing completed samples"):
                sample_id = futures[future]
                try:
                    res = future.result()
                except AssertionError as e:
                    print(f"Error processing file {sample_id}: {e}")
                    continue
                append_list.append(res["df"])
                stat_list.append(res["stat"])

        long_df = pd.concat(append_list, axis=0)
        stats_df = pd.DataFrame(stat_list)
        long_df.to_csv(output.long_df, sep="\t", index=False)
        stats_df.to_csv(output.stats, sep="\t", index=False)

        wide_df = long_to_wide_ordered(long_df)
        marker_stats = marker_missingness(wide_df)
        marker_stats.to_csv(output.marker_stats, sep="\t", index=False)
        final_wide_df = wide_sample_to_wide_haplotype(wide_df)
        final_wide_df.to_csv(output.wide_hap_df, sep="\t", index=False)


rule geo_mhaps:
    input:
        long_df = f"{outdir}/geo/long.tsv",
        stats = f"{outdir}/geo/stats.tsv",
        wide_hap_df = f"{outdir}/geo/wide_haplotypes.tsv",
        marker_stats = f"{outdir}/geo/marker_stats.tsv",
