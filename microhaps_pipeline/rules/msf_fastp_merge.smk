wildcard_constraints:
    marker=f"({'|'.join(Markers)})"

def fastp_merge_reads(input_R1, input_R2, output_merged, log_json, log_html, params_):
    return shell(f"""
        fastp -A -Q -L -m --in1 {input_R1} --in2 {input_R2} --merged_out {output_merged} -j {log_json} -h {log_html} {params_}
    """)

def vsearch_filter_insert_length(input_merged, output_merged_fasta, params_, log_ = None):
    maxEE = params_.get("maxEE", 5)
    min_max = params_.get("min_max", (0, 0))
    vsearch_filter_log = "/dev/null" if log_ is None or log_ == "" else log_
    return shell(f"""
        vsearch --fastx_filter {input_merged} --fastaout {output_merged_fasta} --fastq_maxee {maxEE} \
        --fastq_minlen {min_max[0]} --fastq_maxlen {min_max[1]} > {vsearch_filter_log} 2>&1
    """)

def vsearch_filter_unique(input_merged_fasta, output_uniq_fasta, sample):
    return shell(f"""
        vsearch --fastx_uniques {input_merged_fasta} --fastaout {output_uniq_fasta} --fasta_width 0 --sizeout --sample {sample}
    """)

def vsearch_unoise(input_merged_fasta, output_denoised_fasta, sample, params_, log_ = None):
    minsize_ratio = params_.get('minsize_ratio', 1) # final minsize if minsize_ratio < 1 = max(minsize_ratio * total size, minsize)
    minsize = params_.get("minsize", 8)
    unoise_log = "/dev/null" if log_ is None or log_ == "" else log_
    unoise_alpha = params_.get("unoise_alpha", 2)

    import math
    import sys
    def get_size_from_fasta(fasta):
        with open(fasta, 'r') as f:
            for line in f:
                if line.startswith(">"):
                    match = re.search(r"size=(\d+)", line)
                    if match:
                        yield int(match.group(1))
    if minsize_ratio < 1:
        total_read_in_fa = sum([read_size for read_size in get_size_from_fasta(input_merged_fasta)])
        unoise_minsize = max( math.floor(total_read_in_fa * minsize_ratio), minsize)
        print(f"Dynamic unoise minsize for {input_merged_fasta}: {unoise_minsize}", file=sys.stderr)
    else:
        unoise_minsize = params_.minsize
        
    return shell(f"""
        vsearch --cluster_unoise {input_merged_fasta} --centroids {output_denoised_fasta} \
        --fasta_width 0  --minsize {unoise_minsize} --unoise_alpha {unoise_alpha} \
        --threads 1 > {unoise_log} 2>&1
        """)


def vsearch_uchime(input_denoised_fasta, output_uchime_fasta):
    return shell(f"""
        vsearch --uchime3_denovo {input_denoised_fasta} --nonchimeras {output_uchime_fasta} --fasta_width 0 --qmask none --threads 1
    """)

rule merge_sample_marker_denoise:
    threads: 3
    input:
        f"{outdir}/samples/{{sample}}/markers-reads/.quality_trimmed"
    output:
        f"{outdir}/samples/{{sample}}/markers-reads/.merged_denoised"
    params:
        fastp_merge_params = config.get("fastp_merge_params", ""),
        log_dir = f"{outdir}/samples/{{sample}}/logs/",
        maxEE = config.get("maxEE", 5),
        insert_size_tolerant_threshold = config.get("insert_size_tolerant_threshold", 0.25),
        minsize = config.get("unoise_minsize", 8),  # minimum size of clusters to be kept
        minsize_ratio = config.get('unoise_minsize_ratio', 1), # final minsize if minsize_ratio < 1 = max(minsize_ratio * total size, minsize)
        unoise_alpha = config.get("unoise_alpha", 2),  # alpha parameter for unoise algorithm
    run:
        from concurrent.futures import ThreadPoolExecutor
        import os

        marker_read_dir = os.path.dirname(input[0])
        
        def get_min_max_len(marker):
            tolerant_threshold = params.insert_size_tolerant_threshold
            import pandas as pd
            all_markers = pd.read_table(targetregion_file, header=None, names=["Chr", "Start", "End", "Amplicon_name"])
            all_markers["length"] = all_markers["End"] - all_markers["Start"] + 1
            rel_marker = all_markers.query("Amplicon_name == @marker")["length"].values[0]
            min_len = int(rel_marker * (1 - tolerant_threshold))
            max_len = int(rel_marker * (1 + tolerant_threshold))
            return min_len, max_len

        all_markers = list(Markers)
        logging_dir_marker = [
            os.path.join(params.log_dir, marker) for marker in all_markers
        ]
        fastp_merge_args = [
            {
                "input_R1": os.path.join(marker_read_dir, marker, "trimmed-filtered_R1.fastq.gz"),
                "input_R2": os.path.join(marker_read_dir, marker, "trimmed-filtered_R2.fastq.gz"),
                "output_merged": os.path.join(marker_read_dir, marker, "merged.fastq.gz"),
                "log_json": os.path.join(params.log_dir, marker, "fastp_merge.json"),
                "log_html": os.path.join(params.log_dir, marker, "fastp_merge.html"),
                "params_": params.fastp_merge_params
            } for marker in all_markers
        ]
        min_max_inserts_length = [get_min_max_len(marker) for marker in all_markers]
        vsearch_filter_inserts_args = [
            {
                "input_merged": os.path.join(marker_read_dir, marker, "merged.fastq.gz"),
                "output_merged_fasta": os.path.join(marker_read_dir, marker, "merged.fasta"),
                "params_": {
                    "maxEE": params.maxEE,
                    "min_max": min_max_inserts_length[i]
                },
                "log_": os.path.join(params.log_dir, marker, "vsearch_filter.log")
            } for i, marker in enumerate(all_markers)
        ]
        vsearch_filter_unique_args = [
            {
                "input_merged_fasta": os.path.join(marker_read_dir, marker, "merged.fasta"),
                "output_uniq_fasta": os.path.join(marker_read_dir, marker, "merged_uniq.fasta"),
                "sample": wildcards.sample
            } for marker in all_markers
        ]
        vsearch_unoise_args = [
            {
                "input_merged_fasta": os.path.join(marker_read_dir, marker, "merged_uniq.fasta"),
                "output_denoised_fasta": os.path.join(marker_read_dir, marker, "merged_denoised.fasta"),
                "sample": wildcards.sample,
                "params_": params,
                "log_": os.path.join(params.log_dir, marker, "vsearch_unoise.log")
            } for marker in all_markers
        ]
        vsearch_uchime_args = [
            {
                "input_denoised_fasta": os.path.join(marker_read_dir, marker, "merged_denoised.fasta"),
                "output_uchime_fasta": os.path.join(marker_read_dir, marker, "final_merged.denoised_nochimera.fasta")
            } for marker in all_markers
        ]

        def process_per_marker(marker_log_dir, fastp_merge_args, vsearch_filter_inserts_args,
            vsearch_filter_unique_args, vsearch_unoise_args, vsearch_uchime_args):
            completion_step = 0

            # 0. prepare log dir
            shell(f"mkdir -p {marker_log_dir}")
            completion_step += 1

            # 1. fastp merge reads
            fastp_merge_reads(**fastp_merge_args)
            completion_step += 1

            # 2. (get_min_max_len) + vsearch filter insert lengths
            vsearch_filter_insert_length(**vsearch_filter_inserts_args)
            completion_step += 1

            # 3. vsearch filter unique count occurrence
            vsearch_filter_unique(**vsearch_filter_unique_args)
            completion_step += 1

            # 4. (determine_minsize_ratio) + vsearch unoise denoise
            vsearch_unoise(**vsearch_unoise_args)
            completion_step += 1

            # 5. vsearch uchime denovo
            vsearch_uchime(**vsearch_uchime_args)
            completion_step += 1

            return completion_step

        with ThreadPoolExecutor(max_workers=threads) as executor:
            result = list(executor.map(process_per_marker, logging_dir_marker, fastp_merge_args, vsearch_filter_inserts_args,
                vsearch_filter_unique_args, vsearch_unoise_args, vsearch_uchime_args))
        
        if not all([res == 6 for res in result]):
            raise Exception("Error in processing some markers, please check previous logs")

        shell(f"touch {output}")


def get_same_marker_all_sample(wildcards):
    all_markers = []
    for sample in IDs:
        output_dir = checkpoints.bam_to_marker_fastq.get(sample=sample).output[0]
        all_markers.extend(glob_wildcards(os.path.join(output_dir, "{marker}", "target_R1.fastq.gz")).marker)
    all_markers = list(set(all_markers))
    return expand(f"{outdir}/samples/{{sample}}/markers-reads/{wildcards.marker}/final_merged.denoised_nochimera.fasta",
                  sample=IDs)


rule merge_per_marker:
    priority: 10
    threads: 2
    input:
        # expand(f"{outdir}/samples/{{sample}}/markers-reads/complete.flag", sample=IDs),
        expand(f"{outdir}/samples/{{sample}}/markers-reads/.merged_denoised", sample=IDs),
    output:
        merged_denoise_nochim_dir =directory(f"{outdir}/malamp/merged/"),
        merged_flag = f"{outdir}/malamp/merged/.merged"
    params:
        sample_dir = f"{outdir}/samples"
    run:
        from concurrent.futures import ThreadPoolExecutor
        import os

        output_dir = os.path.dirname(output.merged_flag)
        # 1. making the output directory
        shell(f"mkdir -p {output.merged_denoise_nochim_dir}")

        def process_per_marker(input_, output_, temp_non_uniq, temp_unsorted):
            completed = 0

            # 2. merge every markers
            shell(f"cat {input_} > {temp_non_uniq}")
            completed += 1

            # 3. filter to unique sequences
            shell(f"vsearch --fastx_uniques {temp_non_uniq} --fastaout {temp_unsorted} --fasta_width 0 --sizein --sizeout --relabel_sha1")
            completed += 1

            # 4. sort by size
            shell(f"vsearch --sortbysize {temp_unsorted} --output {output_}")
            completed += 1

            # 5. clean up
            shell(f"rm {temp_non_uniq} {temp_unsorted}")
            completed += 1
            return completed
        
        
        all_args = []
        for marker in Markers:
            all_inputs = [os.path.join(params.sample_dir, sample_, "markers-reads", marker, "final_merged.denoised_nochimera.fasta") for sample_ in IDs]
            all_args.append(
            {
                "input_": " ".join(all_inputs),
                "output_": os.path.join(output_dir, f"final_merged_{marker}.denoised_nochimera.fasta"),
                "temp_non_uniq": os.path.join(output_dir, f"final_merged_{marker}.denoised_nochimera.temp.fasta"),
                "temp_unsorted": os.path.join(output_dir, f"final_merged_{marker}.denoised_nochimera.temp1.fasta")
            })
        
        with ThreadPoolExecutor(max_workers=threads) as executor:
            result = list(executor.map(lambda args: process_per_marker(**args), all_args))

        if not all([res == 4 for res in result]):
            raise Exception("Error in processing some markers, please check previous logs")

        shell(f"touch {output.merged_flag}")

rule create_otutab_per_sample:
    threads: 1
    input:
        marker_merged_flag = f"{outdir}/malamp/merged/.merged",
        sample_merged_flag = f"{outdir}/samples/{{sample}}/markers-reads/.merged_denoised",
    output:
        otutab = f"{outdir}/samples/{{sample}}/markers-reads/seqtab.tsv",
    params:
        sequence_to_denoise_criteria = config.get('sequence_to_denoise_criteria', "--id 0.98 --iddef 1 --target_cov 0.98 --query_cov 0.98")
    run:
        import os
        merged_marker_dir = os.path.dirname(input.marker_merged_flag)
        sample_dir = os.path.dirname(input.sample_merged_flag)

        per_marker_otu = []
        for marker in Markers:
            fasta = os.path.join(merged_marker_dir, f"final_merged_{marker}.denoised_nochimera.fasta")
            sample_fa = os.path.join(sample_dir, marker, "merged_uniq.fasta")
            marker_otutab = os.path.join(sample_dir, marker, "seqtab.tsv")
            shell(f"""vsearch --usearch_global {sample_fa} --db {fasta} {params.sequence_to_denoise_criteria} \
            --otutabout {marker_otutab} \
            --maxaccepts 0 --maxrejects 0 --maxhits 1 \
            --sizein --sizeout --fasta_width 0 --qmask none --dbmask none --threads {threads}
            """)
            per_marker_otu.append(marker_otutab)

        import pandas as pd
        all_dfs = [pd.read_table(f, sep="\t", header=0) for f in per_marker_otu]
        to_concat = [df for df in all_dfs if df.shape[0] > 0]
        if len(to_concat) == 0:
            # All empty
            otutab = pd.DataFrame(columns=["#OTU ID", wildcards.sample])
        else:
            otutab = pd.concat(to_concat, axis=0, ignore_index=True).fillna(0)
            if otutab.shape[1] == 2:
                otutab.columns = ["#OTU ID", wildcards.sample]
            else:
                otutab.loc[:, wildcards.sample] = 0
            otutab.iloc[:, 1] = otutab.iloc[:, 1].astype(int)
        otutab.to_csv(output.otutab, sep="\t", header=True, index=None)

rule merge_otutab_joint:
    input:
        otutabs = expand(f"{outdir}/samples/{{sample}}/markers-reads/seqtab.tsv", sample=IDs),
        merged_flag = f"{outdir}/malamp/merged/.merged"
    output:
        otutab = f"{outdir}/malamp/fastp/seqtab.tsv",
        filtered_fasta = f"{outdir}/malamp/fastp/dereplicated_counted_no_chimeras-unique.fa",
    run:
        denoise_fa = " ".join([f"{outdir}/malamp/merged/final_merged_{marker}.denoised_nochimera.fasta" for marker in Markers])
        shell(f"""
        cat {denoise_fa} > {output.filtered_fasta} &&
        python {microhaps_basedir}/scripts/merge_otutab.py \
         --denoised_fasta {output.filtered_fasta} \
         --out_seqtab {output.otutab} \
         --indv_seqtabs {input.otutabs}
        """)
