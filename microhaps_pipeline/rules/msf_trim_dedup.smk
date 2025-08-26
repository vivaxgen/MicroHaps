new_postprocess = config.get('post_process', "old")
primers_trimmed = config.get("primers_trimmed", False)

if not primers_trimmed:
    # trim primers from reads per index sample
    rule trim_primers:
        group: "trimming_primer_qual"
        threads: 3
        input:
            R1 = "{pfx}/target_R1.fastq.gz",
            R2 = "{pfx}/target_R2.fastq.gz",
            prim_fw = primer_fw_file,
            prim_rv = primer_rev_file,
        output:
            R1 = "{pfx}/primer-trimmed_R1.fastq.gz",
            R2 = "{pfx}/primer-trimmed_R2.fastq.gz",
        log:
            "{pfx}/logs/trim_before_dada2.log",
        params:
            platform = config['platform'],
            trim_qv = config['trimqv'],
            additional_params = f"--nextseq-trim={config.get('trimqv', 20)}" if is_nextseq_or_novaseq() else ""
        shell: 
            "cutadapt -g file:{input.prim_fw} -G file:{input.prim_rv} -o {output.R1} -p {output.R2}"
            " --pair-adapters --discard-untrimmed {params.additional_params} --action=trim {input.R1} {input.R2} > {log} 2>&1"
else:
    rule pseudo_trim_primers:
        group: "trimming_primer_qual"
        localrule: True
        input:
            R1 = "{pfx}/target_R1.fastq.gz",
            R2 = "{pfx}/target_R2.fastq.gz",
        output:
            R1 = "{pfx}/primer-trimmed_R1.fastq.gz",
            R2 = "{pfx}/primer-trimmed_R2.fastq.gz",
        shell:
            "ln -s {input.R1} {output.R1} && ln -s {input.R2} {output.R2}"

# Dedup now performed prior to getting the reads back from mapping
# rule optical_dedup:
#     input:
#         R1 = "{pfx}/primer-trimmed_R1.fastq.gz",
#         R2 = "{pfx}/primer-trimmed_R2.fastq.gz",
#     output:
#         R1 = "{pfx}/Odedup-primer-trimmed_R1.fastq.gz",
#         R2 = "{pfx}/Odedup-primer-trimmed_R2.fastq.gz",
#     log:
#         clump_log = "{pfx}/logs/clumpify.log",
#     params:
#         instrument_specific = "spany adjacent" if config.get("instrument", "generic") == "nextseq" else "",
#     shell:
#         """
#         clumpify.sh in1={input.R1} in2={input.R2} out1={output.R1} out2={output.R2} dedupe optical {params.instrument_specific} > {log.clump_log} 2>&1
#         """

if config.get('additional_trim_before_haplotype_generation', False) == True or config.get('merge_map', 'dada2') == 'dada2':
    rule trim_before_haplotype_generation:
        group: "trimming_primer_qual"
        threads: 3
        input:
            R1 = "{pfx}/primer-trimmed_R1.fastq.gz",
            R2 = "{pfx}/primer-trimmed_R2.fastq.gz",
        output:
            R_trimmed_R1 = temp("{pfx}/T-primer-trimmed_R1.fastq.gz"),
            R_trimmed_R2 = temp("{pfx}/T-primer-trimmed_R2.fastq.gz"),
            Filtered_final_R1 = "{pfx}/trimmed-filtered_R1.fastq.gz",
            Filtered_final_R2 = "{pfx}/trimmed-filtered_R2.fastq.gz",
        params:
            trim_right_R1 = config.get("trim_right_R1", 5),
            trim_right_R2 = config.get("trim_right_R2", 5),
            cut_window = config.get("cut_window", 1),
            window_minqual = config.get("truncQ", 20),
            min_length = config.get("min_length", 30),
            maxEE = config.get("maxEE", 5)
        log:
            fastp_json = "{pfx}/logs/fastp_tr_before_hapgen.json",
            fastp_html = "{pfx}/logs/fastp_tr_before_hapgen.html",
            vsearch_filter_log = "{pfx}/logs/fastx_mee.log",
        shell:
            """
            fastp -A -w {threads} -i {input.R1} -I {input.R2} -o {output.R_trimmed_R1} -O {output.R_trimmed_R2} \
                --trim_tail1 {params.trim_right_R1} --trim_tail2 {params.trim_right_R2} \
                --cut_tail --cut_tail_window_size {params.cut_window} --cut_tail_mean_quality {params.window_minqual} \
                -l {params.min_length} \
                -j {log.fastp_json} -h {log.fastp_html} &&
            vsearch --fastx_filter {output.R_trimmed_R1} --reverse {output.R_trimmed_R2} --fastqout >(gzip -c  > {output.Filtered_final_R1}) \
                --fastqout_rev >(gzip -c > {output.Filtered_final_R2}) --fastq_maxee {params.maxEE} > {log.vsearch_filter_log} 2>&1
            """
else:
    rule pseudo_trim_before_haplotype_generation:
        group: "trimming_primer_qual"
        localrule: True
        input:
            R1 = "{pfx}/primer-trimmed_R1.fastq.gz",
            R2 = "{pfx}/primer-trimmed_R2.fastq.gz",
        output:
            Filtered_final_R1 = "{pfx}/trimmed-filtered_R1.fastq.gz",
            Filtered_final_R2 = "{pfx}/trimmed-filtered_R2.fastq.gz",
        shell:
            "ln -s {input.R1} {output.Filtered_final_R1} && ln -s {input.R2} {output.Filtered_final_R2}"