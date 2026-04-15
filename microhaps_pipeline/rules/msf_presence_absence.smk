
presence_absence_bed = get_abspath(config.get("presence_absence_markers"))

rule final_presence_absence_report:
    input:
        f"{outdir}/malamp/presence_absence.tsv"


rule merge_presence_absence_report:
    input:
        reports = expand(f"{outdir}/samples/{{sample}}/presence_absence/stats.tsv", sample=IDs),
    output:
        f"{outdir}/malamp/presence_absence.tsv"
    log:
        f"{outdir}/logs/merge_presence_absence_report.log"
    params:
        min_coverage = config.get("presence_absence_min_coverage", 90),
        min_numreads = config.get("presence_absence_min_numreads", 5),
    run:
        import pandas as pd
        all_result = [pd.read_table(f) for f in input.reports]
        full_result = pd.concat(all_result)
        full_result["presence"] = "absent"
        full_result.loc[(full_result["coverage"] >= params.min_coverage) &
                        (full_result["numreads"] >= params.min_numreads), "presence"] = "present"
        full_result["coverage|nread"] = full_result["coverage"].astype(str) + "|" + full_result["numreads"].astype(str)
        full_result.to_csv(log[0], index=False, sep="\t")
        final_result = full_result.pivot_table(index=["Chr", "Start", "End", "Amplicon_name"], columns="sample", values="presence", aggfunc="first", fill_value="absent").reset_index()
        final_result.to_csv(output[0], index=False, sep="\t")


use rule final_bam_depth_coverage_per_inserts as generate_stats_for_presence_absence with:
    input:
        bam = f"{outdir}/samples/{{sample}}/maps/index_sorted_presence_absence.bam",
        bai = f"{outdir}/samples/{{sample}}/maps/index_sorted_presence_absence.bam.bai",
        targetregion_file = presence_absence_bed,
    output:
        depth_coverage = f"{outdir}/samples/{{sample}}/presence_absence/stats.tsv",
    params:
        suffix_to_remove = "/maps/index_sorted_presence_absence.bam",
    

rule filter_and_index_bam_presence_absence:
    input:
        merged_bam = f"{outdir}/samples/{{sample}}/maps/final.bam",
        bed = presence_absence_bed
    output:
        index_sorted = f"{outdir}/samples/{{sample}}/maps/index_sorted_presence_absence.bam",
        rejects = f"{outdir}/samples/{{sample}}/maps/rejects_presence_absence.bam"
    log:
        f"{outdir}/samples/{{sample}}/logs/filter_and_sort_presence_absence.log"
    params:
        primers_bed = get_abspath(config.get("primers_bed")),
        samtools_ampliconclip_params = "--strand --clipped"
    shell:
        """
        samtools view -h -L {input.bed} {input.merged_bam} | 
        samtools ampliconclip -b {params.primers_bed} -f {log} {params.samtools_ampliconclip_params} --rejects-file {output.rejects} - |
        samtools fixmate -mMu -  - |
        samtools calmd -u - {refseq} - |
        samtools sort -o {output.index_sorted}
        """
