# /home/data/malaria/work/ldwg/LATEST_NGSPL/vvg-MicroHaps/bin/activate
# snakemake --snakefile drug_resistance_workflow.smk


drugs_resistance_variant_list = get_abspath(config.get("drugs_resistance_variant_list"))
aa_pos_file = get_abspath(config.get("drugs_resistance_aa_pos"))
gff = get_abspath(config.get("gff_file"))
target_bed = get_abspath(config.get("drugs_resistance_bed"))

rule final_drug_resistance_report:
    input:
        f"{outdir}/malamp/drug_resistance.tsv",
        f"{outdir}/samples/{{sample}}/drugs/drug_resistance_stats.tsv"

rule merge_drug_resistance_report:
    input:
        reports = expand(f"{outdir}/samples/{{sample}}/drugs/drug_resistance_report.tsv", sample=IDs),
        variant_list = drugs_resistance_variant_list,
    output:
        f"{outdir}/malamp/drug_resistance.tsv"
    run:
        import pandas as pd
        variants_ = pd.read_table(input.variant_list, header=None, names=["chrom", "pos0", "pos", "variant_id"])

        columns = ["sample"] + sorted(variants_["variant_id"].unique().tolist(), key=lambda x: (x.split(":")[0], int(x.split(":")[1])))
        all_reports = []
        for report_ in input.reports:
            report_df = pd.read_table(report_)
            result = [report_df["sample"].iloc[0]]
            for col in columns[1:]:
                if not col in report_df["marker"].values:
                    result.append(".|0")
                else:
                    res = report_df.query("marker == @col").sort_values("count", ascending=False)
                    aas = ','.join(res['aa'].values)
                    counts = ','.join([str(a) for a in res['count'].values])
                    result.append(f"{aas}|{counts}")
            all_reports.append(result)
        combined_report_df = pd.DataFrame(all_reports, columns = columns).to_csv(output[0], index=False, sep="\t")

rule predict_aa:
    input:
        reference = refseq,
        gff = gff,
        pseudohaplotypes =  f"{outdir}/samples/{{sample}}/drugs/haplotype_pseudohaplotypes.tsv",
        variant_list = drugs_resistance_variant_list
    output:
        protein_predictions = f"{outdir}/samples/{{sample}}/drugs/haplotype_protein_predictions.tsv",
        final_output = f"{outdir}/samples/{{sample}}/drugs/drug_resistance_report.tsv"
    params:
        min_support = 5,
    shell:
        """
        python {microhaps_basedir}/scripts/predict_aa.py --reference {input.reference} --gff_file {input.gff} --pseudohaplotypes {input.pseudohaplotypes} --variant_list {input.variant_list} --sample {wildcards.sample} --min_support {params.min_support} --protein_predictions {output.protein_predictions} --drug_resistance_report {output.final_output}
        """

rule generate_pseudohaplotype:
    input:
        filtered_bam = f"{outdir}/samples/{{sample}}/maps/named_sorted_drugs.bam", 
        variant_list = drugs_resistance_variant_list,
    output:
        details = f"{outdir}/samples/{{sample}}/drugs/haplotype_details.tsv",
        result = f"{outdir}/samples/{{sample}}/drugs/haplotype_pseudohaplotypes.tsv",
    shell:
        """
        ngs-pl construct-pseudo-haplotypes --variants_list {input.variant_list} --pileup --pileout {output.details} --sample {wildcards.sample} -o {output.result} {input.filtered_bam}
        """

rule generate_stats_for_drug_resistance:
    input:
        filtered_bam = f"{outdir}/samples/{{sample}}/maps/index_sorted_drugs.bam",
        target_bed = target_bed,
    output:
        stats = f"{outdir}/samples/{{sample}}/drugs/drug_resistance_stats.tsv",
    shell:
        """
        samtools bedcov -H -d 10 -c {input.target_bed} {input.filtered_bam} > {output.stats}
        """

rule name_sort_bam:
    input:
        bam = f"{outdir}/samples/{{sample}}/maps/index_sorted_drugs.bam",
    output:
        named_sorted = f"{outdir}/samples/{{sample}}/maps/named_sorted_drugs.bam",
    shell:
        """
        samtools sort -n -o {output.named_sorted} {input.bam}
        """

rule filter_and_sort_bam:
    input:
        merged_bam = f"{outdir}/samples/{{sample}}/maps/final.bam",
        bed = target_bed,
    output:
        index_sorted = f"{outdir}/samples/{{sample}}/maps/index_sorted_drugs.bam",
        bai = f"{outdir}/samples/{{sample}}/maps/index_sorted_drugs.bam.bai",
        rejects = f"{outdir}/samples/{{sample}}/maps/rejects_drugs.bam",
    log:
        f"{outdir}/samples/{{sample}}/logs/filter_and_sort_drugs.log"
    params:
        primers_bed = get_abspath(config.get("primers_bed")),
        samtools_ampliconclip_params = "--strand --clipped",
    shell:
        """
        samtools view -h -L {input.bed} {input.merged_bam} | 
        samtools ampliconclip -b {params.primers_bed} -f {log} {params.samtools_ampliconclip_params} --rejects-file {output.rejects} - |
        samtools sort -n |
        samtools fixmate -mMu -  - |
        samtools calmd -u - {refseq} - |
        samtools sort -o {output.index_sorted}  && samtools index {output.index_sorted}
        """

rule generate_variant_list:
    localrule: True
    input:
        aa_pos_file = aa_pos_file,
        reference = refseq,
        gff = gff,
    output:
        drugs_resistance_variant_list
    shell:
        """
        python {microhaps_basedir}/scripts/aa_to_genomic.py --aa_pos_file {input.aa_pos_file} --reference {input.reference} --gff_file {input.gff} --output {output}
        """