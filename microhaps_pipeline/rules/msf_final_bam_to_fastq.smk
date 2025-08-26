

# see https://www.htslib.org/algorithms/duplicate.html
rule optical_dedup:
    group: "filter_dedup"
    threads: 4
    input:
        bam = f"{outdir}/samples/{{sample}}/maps/final.bam"
    output:
        marked = f"{outdir}/samples/{{sample}}/maps/final.temp.mark_dup.bam",
        deduped = f"{outdir}/samples/{{sample}}/maps/final.temp.op_dedup.bam",
        stats = f"{outdir}/samples/{{sample}}/maps/dedup_stat.txt"
    params:
        pixel_distance = config.get("optical_dedup_pixel_distance", 10) # 10 to be conservative, values dependent on platform
    shell:
        # mode s -  measure positions based on sequence start.
        """
        samtools collate -@ {threads} -O -u {input.bam} | samtools fixmate -@ {threads} -m -u - - | samtools sort -@ {threads} -u - | samtools markdup -@ {threads} -f {output.stats} -m s -d {params.pixel_distance} - {output.marked} &&
        samtools view -b -e '[dt]!="SQ"' {output.marked} | samtools sort -@ {threads} -o {output.deduped}
        """

rule filter_bam:
    group: "filter_dedup"
    threads: 4
    input:
        tempbam = f"{outdir}/samples/{{sample}}/maps/final.temp.op_dedup.bam",
    output:
        unsorted = f"{outdir}/samples/{{sample}}/maps/unsorted.final.filtered.bam",
        filtered = f"{outdir}/samples/{{sample}}/maps/final.filtered.bam",
    log:
        filter_stats = f"{outdir}/samples/{{sample}}/logs/bam_to_fastq_filter_stats.log"
    shell:
        "ngs-pl filter-reads-orientation --remove_FF --remove_RR --remove_RF --remove_trans --remove_unmapped --remove_secondary --remove_supplementary --max-insert-length 350 -o {output.unsorted} --outstat {log.filter_stats} {input.tempbam} && "
        "samtools sort -@ {threads} -o {output.filtered} {output.unsorted}"

rule bam_to_fastq:
    group: "filter_dedup"
    threads: 4
    input: 
        tempbam = f"{outdir}/samples/{{sample}}/maps/final.filtered.bam",
    output:
        R1 = f"{outdir}/samples/{{sample}}/mhaps-reads/target_R1.fastq.gz",
        R2 = f"{outdir}/samples/{{sample}}/mhaps-reads/target_R2.fastq.gz",
    shell:
        "samtools collate -u -O {input.tempbam} "
        " | samtools fastq --thread {threads} -1 {output.R1} -2 {output.R2} -0 /dev/null -s /dev/null -n"

rule bam_to_marker_fastq:
    group: "filter_dedup"
    threads: 4
    input:
        tempbam = f"{outdir}/samples/{{sample}}/maps/final.filtered.bam",
        tempbam_bai = f"{outdir}/samples/{{sample}}/maps/final.filtered.bam.bai",
        insertseq_file = targetregion_file,
    output:
        filtered = [f"{outdir}/samples/{{sample}}/markers-reads/{marker}/final.filtered.bam" for marker in Markers],
        collated = [f"{outdir}/samples/{{sample}}/markers-reads/{marker}/final.temp.op_dedup.collated.bam" for marker in Markers],
        markers_target = [f"{outdir}/samples/{{sample}}/markers-reads/{marker}/target_R{read1_2}.fastq.gz" for marker in Markers for read1_2 in [1,2]],
        # complete_flag = f"{outdir}/samples/{{sample}}/markers-reads/complete.flag",
    params:
        per_marker_directory = f"{outdir}/samples/{{sample}}/markers-reads"
    run:
        shell(f"mkdir -p {params.per_marker_directory}")
        print(f"mkdir -p {params.per_marker_directory}")
        import pandas as pd
        from concurrent.futures import ThreadPoolExecutor

        markers = pd.read_table(input.insertseq_file, header=None, names=["Chr", "Start", "End", "Amplicon_name"])
        markers["region"] = markers["Chr"] + ":" + markers["Start"].astype(str) + "-" + markers["End"].astype(str)

        args = []
        for _, row in markers.iterrows():
            marker = row["Amplicon_name"]
            args.append({
                "bamfile": input.tempbam,
                "region": row["region"],
                "marker": marker,
                "output_R1": f"{params.per_marker_directory}/{marker}/target_R1.fastq.gz",
                "output_R2": f"{params.per_marker_directory}/{marker}/target_R2.fastq.gz",
                "filtered": f"{params.per_marker_directory}/{marker}/final.filtered.bam",
                "collated": f"{params.per_marker_directory}/{marker}/final.temp.op_dedup.collated.bam",
                "dir_marker": f"{params.per_marker_directory}/{marker}"
            })

        def process_marker(args):
            bamfile = args["bamfile"]
            region = args["region"]
            output_R1 = args["output_R1"]
            output_R2 = args["output_R2"]
            filtered = args["filtered"]
            collated = args["collated"]
            dir_marker = args["dir_marker"]
            shell(f"mkdir -p {dir_marker} && samtools view -b {bamfile} \"{region}\" | ngs-pl filter-reads-orientation --remove_unmapped -o {filtered} - ")
            shell(f"samtools collate -o {collated} {filtered}")
            shell(f"samtools fastq -1 {output_R1} -2 {output_R2} {collated}")

        nworker = int(threads)
        print(nworker)
        with ThreadPoolExecutor(max_workers=nworker) as executor:
            executor.map(process_marker, args)
        print(f"Completed processing markers for sample {wildcards.sample}")


rule final_bam_depth_coverage_per_inserts:
    input:
        bam = f"{outdir}/samples/{{sample}}/maps/final.bam",
        bai = f"{outdir}/samples/{{sample}}/maps/final.bam.bai",
        insertseq_file = targetregion_file,
    output:
        depth_coverage = f"{outdir}/samples/{{sample}}/logs/depth_coverage-mapped.tsv",
    params:
        prefix_to_remove = f"{outdir}/samples/",
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
        per_inserts_depth_coverage = expand(f"{outdir}/samples/{{sample}}/logs/depth_coverage-mapped.tsv", sample=IDs),
        insertseq_file = targetregion_file,
    output:
        depth = f"{outdir}/depths-mapped.tsv",
        coverage = f"{outdir}/coverages-mapped.tsv"
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


rule depth_heatmap:
    input:
        depth = f"{outdir}/depths-mapped.tsv",
    output:
        depth = f"{outdir}/depths-mapped.png",
        marker = f"{outdir}/markers-mapped.png",
    shell:
        "ngs-pl tab-to-plots --index-column Amplicon_name --additional-title 'Reads'"
        "  --outheatmap {output.depth} --outmarker {output.marker} {input.depth}"


# EOF
