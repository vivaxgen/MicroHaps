
# see https://www.htslib.org/algorithms/duplicate.html
# mode s -  measure positions based on sequence start.
def optical_dedup(input_bam, output_stats, output_marked_bam, output_deduped_bam, pixel_distance, nthread):
    return shell(f"""
        samtools collate -@ {nthread} -O -u {input_bam} \
        | samtools fixmate -@ {nthread} -m -u - - \
        | samtools sort -@ {nthread} -u - \
        | samtools markdup -@ {nthread} -f {output_stats} -m s -d {pixel_distance} - {output_marked_bam} \
        && samtools view -b -e '[dt]!="SQ"' {output_marked_bam} | samtools sort -@ {nthread} -o {output_deduped_bam}
    """)

def filter_bam(input_bam, output_unsorted_bam, output_filtered_bam, log_file, nthread, max_insert = 350):
    return shell(f"""ngs-pl filter-reads-orientation --remove_FF --remove_RR --remove_RF --remove_trans \
    --remove_unmapped --remove_secondary --remove_supplementary --max-insert-length {max_insert} \
    -o {output_unsorted_bam} --outstat {log_file} {input_bam} \
    && samtools sort -@ {nthread} -o {output_filtered_bam} {output_unsorted_bam}
    """)

def bam_to_fastq(input_bam, output_R1, output_R2, nthread):
    return shell(f"""
        samtools collate -u -O {input_bam} \
        | samtools fastq --thread {nthread} -1 {output_R1} -2 {output_R2} -0 /dev/null -s /dev/null -n
    """)

def filter_bam_chrom(input_bam, chrom, output_filtered_bam, output_filtered_bam_inverse):
    cmds = f"""samtools view -e 'rname=="{chrom}"' -b {input_bam} -o {output_filtered_bam} && \
        samtools view -e 'rname!="{chrom}"' -b {input_bam} -o {output_filtered_bam_inverse}"""
    return shell(cmds)

rule optical_dedup_filter_sample:
    threads: 4
    input:
        bam = f"{outdir}/samples/{{sample}}/maps/final.bam"
    output:
        marked = temp(f"{outdir}/samples/{{sample}}/maps/final.temp.mark_dup.bam"),
        deduped = temp(f"{outdir}/samples/{{sample}}/maps/final.temp.op_dedup.bam"),
        stats = f"{outdir}/samples/{{sample}}/maps/dedup_stat.txt",
    params:
        pixel_distance = config.get("optical_dedup_pixel_distance", 10) # 10 to be conservative, values dependent on platform
    run:
        optical_dedup(input_bam=input.bam, output_stats=output.stats, output_marked_bam=output.marked, output_deduped_bam=output.deduped, pixel_distance=params.pixel_distance, nthread=threads)


rule bam_to_sample_fastq:
    threads: 4
    input:
        deduped = f"{outdir}/samples/{{sample}}/maps/final.temp.op_dedup.bam",
        deduped_index = f"{outdir}/samples/{{sample}}/maps/final.temp.op_dedup.bam.bai",
    output:
        filtered_unsorted = temp(f"{outdir}/samples/{{sample}}/maps/final.temp.filtered_unsorted.bam"),
        filtered = temp(f"{outdir}/samples/{{sample}}/maps/final.temp.filtered.bam"),
        R1 = f"{outdir}/samples/{{sample}}/mhaps-reads/target_R1.fastq.gz",
        R2 = f"{outdir}/samples/{{sample}}/mhaps-reads/target_R2.fastq.gz",
    log:
        filter_stats = f"{outdir}/logs/{{sample}}/bam_to_fastq_filter_stats.log"
    run:
        filter_bam(input_bam=input.deduped, output_unsorted_bam=output.filtered_unsorted, output_filtered_bam=output.filtered, log_file=log.filter_stats, nthread=threads, max_insert = 350)
        bam_to_fastq(input_bam=output.filtered, output_R1=output.R1, output_R2=output.R2, nthread=threads)

rule bam_to_marker_fastq:
    threads: 4
    input:
        tempbam = f"{outdir}/samples/{{sample}}/maps/final.temp.op_dedup.bam",
        tempbam_bai = f"{outdir}/samples/{{sample}}/maps/final.temp.op_dedup.bam.bai",
        insertseq_file = targetregion_file,
    output:
        marker_read_dir = directory(f"{outdir}/samples/{{sample}}/markers-reads"),
        completion_flag = f"{outdir}/samples/{{sample}}/markers-reads/.completed",
    log:
        filter_stats = f"{outdir}/logs/{{sample}}/bam_to_fastq_filter_stats.log"
    params:
        log_dir = f"{outdir}/logs/{{sample}}",
    run:
        shell(f"mkdir -p {output.marker_read_dir}")
        import pandas as pd
        from concurrent.futures import ThreadPoolExecutor
        import sys
        import os

        markers = pd.read_table(input.insertseq_file, header=None, names=["Chr", "Start", "End", "Amplicon_name"])
        markers["region"] = markers["Chr"] + ":" + markers["Start"].astype(str) + "-" + markers["End"].astype(str)

        args = []
        for _, row in markers.iterrows():
            marker = row["Amplicon_name"]
            region = row["region"]
            dir_marker = os.path.join(output.marker_read_dir, marker)
            region_temp_file = os.path.join(dir_marker, "region.temp.bam")
            filtered_region_unsorted_temp_file = os.path.join(dir_marker, "filtered.temp.bam")
            final_filtered_region_temp_file = os.path.join(dir_marker, "final.filtered.temp.bam")
            filter_log = os.path.join(params.log_dir, f"{marker}.filter.log")
            output_R1 = os.path.join(dir_marker, "target_R1.fastq.gz")
            output_R2 = os.path.join(dir_marker, "target_R2.fastq.gz")
            args.append({
                "bamfile": input.tempbam,
                "region": region,
                "dir_marker": f"{output.marker_read_dir}/{marker}",
                "region_temp_file": region_temp_file,
                "filtered_region_unsorted_temp_file": filtered_region_unsorted_temp_file,
                "final_filtered_region_temp_file": final_filtered_region_temp_file,
                "filter_log": filter_log,
                "output_R1": output_R1,
                "output_R2": output_R2,
                "max_insert ": 350
            })

        def process_per_samples(args):
            # 1. make marker subdir
            shell(f"mkdir -p {args['dir_marker']}")
            # 2. samtools view the region to a temp file
            shell(f"samtools view -b {args['bamfile']} \"{args['region']}\" -o {args['region_temp_file']}")
            # 3. filter_bam
            filter_bam(input_bam=args['region_temp_file'], output_unsorted_bam=args['filtered_region_unsorted_temp_file'],
                output_filtered_bam=args['final_filtered_region_temp_file'], log_file=args['filter_log'], nthread=1, max_insert = args['max_insert '])
            # 4. bam_to_fastq
            bam_to_fastq(input_bam=args['final_filtered_region_temp_file'], output_R1=args['output_R1'], output_R2=args['output_R2'], nthread=1)
            # 5. remove temp files
            shell(f"rm -f {args['region_temp_file']} {args['filtered_region_unsorted_temp_file']} {args['final_filtered_region_temp_file']}")
            return 1


        nworker = int(threads)
        with ThreadPoolExecutor(max_workers=nworker) as executor:
            result = list(executor.map(process_per_samples, args))
        
        if not all([res == 1 for res in result ]):
            raise Exception("Error in processing some markers, please check previous logs")
        # Combine all logs
        all_temp_log = " ".join([a["filter_log"] for a in args])
        shell(f"cat {all_temp_log} > {log.filter_stats} && rm -f {all_temp_log}")

        print(f"Completed processing markers for sample {wildcards.sample}", file = sys.stderr)
        shell(f"touch {output.completion_flag}")

rule final_bam_depth_coverage_per_inserts:
    input:
        bam = f"{outdir}/samples/{{sample}}/maps/final.bam",
        bai = f"{outdir}/samples/{{sample}}/maps/final.bam.bai",
        targetregion_file = targetregion_file,
    output:
        depth_coverage = f"{outdir}/samples/{{sample}}/logs/depth_coverage-mapped.tsv",
    params:
        prefix_to_remove = f"{outdir}/samples/",
    run:
        import pandas as pd
        from io import StringIO
        markers = pd.read_table(input.targetregion_file, header=None, names=["Chr", "Start", "End", "Amplicon_name"])
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
        targetregion_file = targetregion_file,
    output:
        depth = f"{outdir}/depths-mapped.tsv",
        coverage = f"{outdir}/coverages-mapped.tsv"
    run:
        import pandas as pd
        from io import StringIO
        all_results = []

        markers = pd.read_table(input.targetregion_file, header=None, names=["Chr", "Start", "End", "Amplicon_name"])
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
