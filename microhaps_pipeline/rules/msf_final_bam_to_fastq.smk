
rule bam_to_fastq:
    threads: 4
    input:
        bam = "{pfx}/samples/{sample}/maps/final.bam"
    output:
        R1 = "{pfx}/samples/{sample}/mhaps-reads/target_R1.fastq.gz",
        R2 = "{pfx}/samples/{sample}/mhaps-reads/target_R2.fastq.gz"
    shell:
        "samtools collate -u -O {input.bam}"
        " | samtools fastq --thread 2 -1 {output.R1} -2 {output.R2} -0 /dev/null -s /dev/null -n"

rule final_bam_depth_coverage_per_inserts:
    input:
        bam = "{pfx}/samples/{sample}/maps/final.bam",
        bai = "{pfx}/samples/{sample}/maps/final.bam.bai",
        insertseq_file = targetregion_file,
    output:
        depth_coverage = "{pfx}/samples/{sample}/logs/final.depth_coverage.tsv",
    params:
        prefix_to_remove = "{pfx}/samples/",
    run:
        import pandas as pd
        from io import StringIO
        markers = pd.read_table(input.insertseq_file, header=None, names=["Chr", "Start", "End", "Amplicon_name"])
        all_results = []
        markers["region"] = markers["Chr"] + ":" + markers["Start"].astype(str) + "-" + markers["End"].astype(str)
        for marker in markers["region"]:
            temp = pd.read_table(StringIO(shell(f"samtools coverage -H -r {marker} {input.bam}", read= True)), header=None, names = ["rname", "startpos", "endpos", "numreads", "covbases", "coverage", "meandepth", "meanbaseq", "meanmapq"])
            temp["sample"] = input.bam.replace(params.prefix_to_remove, "").replace("/maps/final.bam", "")
            temp["region"] = marker
            all_results.append(temp)
        all_results = pd.concat(all_results)
        full_result = markers.merge(all_results, left_on="region", right_on="region", how="outer").drop("region", axis=1)
        full_result.to_csv(output.depth_coverage, sep="\t", index=False)

rule aggregate_depth_coverage:
    input:
        per_inserts_depth_coverage = expand("{{pfx}}/samples/{sample}/logs/final.depth_coverage.tsv", sample=IDs),
        insertseq_file = targetregion_file,
    output:
        depth = "{pfx}/final.depths.tsv",
        coverage = "{pfx}/final.coverages.tsv"
    run:
        import pandas as pd
        from io import StringIO
        all_results = []

        markers = pd.read_table(input.insertseq_file, header=None, names=["Chr", "Start", "End", "Amplicon_name"])
        markers["region"] = markers["Chr"] + ":" + markers["Start"].astype(str) + "-" + markers["End"].astype(str)
        all_result = [pd.read_table(f) for f in input.per_inserts_depth_coverage]
        full_result = pd.concat(all_result)

        # pivot full result to get count of reads per sample per marker
        full_result2 = full_result.pivot_table(index=["Chr", "Start", "End", "Amplicon_name"], columns="sample", values="numreads", fill_value=0).reset_index()
        full_result2[full_result2.columns[4:]] = full_result2[full_result2.columns[4:]].astype(int)
        full_result2.to_csv(output.depth, sep="\t", index=False)

        full_result3 = full_result.pivot_table(index=["Chr", "Start", "End", "Amplicon_name"], columns="sample", values="coverage", fill_value=0).reset_index()
        full_result3[full_result3.columns[4:]] = full_result3[full_result3.columns[4:]].astype(int)
        full_result3.to_csv(output.coverage, sep="\t", index=False)
# EOF
