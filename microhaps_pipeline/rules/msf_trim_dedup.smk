primers_trimmed = config.get("primers_trimmed", False)
quality_trim_required = config.get('additional_trim_before_haplotype_generation', False) == True or config.get('merge_map', 'dada2') == 'dada2'

def trim_primers(input_R1, input_R2, output_R1, output_R2, primer_fwd, primer_rv, additional_params, logfile):
    return shell(f"""cutadapt -g file:{primer_fwd} -G file:{primer_rv} -o {output_R1} -p {output_R2}
         --pair-adapters --discard-untrimmed {additional_params} --action=trim {input_R1} {input_R2} > {logfile} 2>&1
        """)

def pseudo_trim_primers(input_R1, input_R2, output_R1, output_R2):
    return shell(f"ln -s {input_R1} {output_R1} && ln -s {input_R2} {output_R2}")

rule trim_primers:
    threads: 1 if primers_trimmed else 3
    localrule: True if primers_trimmed else False
    input:
        R1 = f"{outdir}/samples/{{sample}}/mhaps-reads/target_R1.fastq.gz",
        R2 = f"{outdir}/samples/{{sample}}/mhaps-reads/target_R2.fastq.gz",
        prim_fw = primer_fw_file,
        prim_rv = primer_rev_file,
    output:
        R1 = f"{outdir}/samples/{{sample}}/mhaps-reads/primer-trimmed_R1.fastq.gz",
        R2 = f"{outdir}/samples/{{sample}}/mhaps-reads/primer-trimmed_R2.fastq.gz",
    log:
        f"{outdir}/samples/{{sample}}/logs/primers_trimming.log",
    params:
        platform = config['platform'],
        trim_qv = config['trimqv'],
        additional_params = f"--nextseq-trim={config.get('trimqv', 20)}" if is_nextseq_or_novaseq() else ""
    run:
        if primers_trimmed: 
            pseudo_trim_primers(input.R1, input.R2, output.R1, output.R2)
        else:
            trim_primers(input.R1, input.R2, output.R1, output.R2, prim_fw, prim_rv, params.additional_params, log)

rule per_marker_trim_primers:
    threads: 1 if primers_trimmed else 3
    localrule: True if primers_trimmed else False
    input:
        prior_completion = f"{outdir}/samples/{{sample}}/markers-reads/.completed",
        prim_fw = primer_fw_file,
        prim_rv = primer_rev_file,
    output:
        trim_completion_flag = f"{outdir}/samples/{{sample}}/markers-reads/.trimmed"
    log:
        f"{outdir}/samples/{{sample}}/logs/primers_trimming.log",
    params:
        platform = config['platform'],
        trim_qv = config['trimqv'],
        additional_params = f"--nextseq-trim={config.get('trimqv', 20)}" if is_nextseq_or_novaseq() else ""
    run:
        import os

        marker_read_dir = os.path.dirname(input.prior_completion)
        all_input_files = [(os.path.join(marker_read_dir, marker, "target_R1.fastq.gz"),
            os.path.join(marker_read_dir, marker, "target_R2.fastq.gz")) for marker in Markers]
        all_output_files = [(os.path.join(marker_read_dir, marker, "trimmed_R1.fastq.gz"),
            os.path.join(marker_read_dir, marker, "trimmed_R2.fastq.gz")) for marker in Markers]

        if primers_trimmed:
            for (inR1, inR2), (outR1, outR2) in zip(all_input_files, all_output_files):
                pseudo_trim_primers(inR1, inR2, outR1, outR2)
        else:
            for (inR1, inR2), (outR1, outR2) in zip(all_input_files, all_output_files):
                trim_primers(inR1, inR2, outR1, outR2, input.prim_fw, input.prim_rv, params.additional_params, log)

        shell(f"touch {output.trim_completion_flag}")

def quality_trim(input_R1, input_R2, output_R1, output_R2, filtered_R1, filtered_R2, nthread, log_vsearch, fastp_html, fastp_json, parameters):
    trim_right_R1 = parameters.get("trim_right_R1", 5)
    trim_right_R2 = parameters.get("trim_right_R2", 5)
    cut_window = parameters.get("cut_window", 1)
    window_minqual = parameters.get("truncQ", 20)
    min_length = parameters.get("min_length", 30)
    maxEE = parameters.get("maxEE", 5)
    return shell(f"""
        fastp -A -w 3 -i {input_R1} -I {input_R2} -o {output_R1} -O {output_R2} \
            --trim_tail1 {trim_right_R1} --trim_tail2 {trim_right_R2} \
            --cut_tail --cut_tail_window_size {cut_window} --cut_tail_mean_quality {window_minqual} \
            -l {min_length} \
            -j {fastp_json} -h {fastp_html} &&
        vsearch --fastx_filter {output_R1} --reverse {output_R2} --fastqout >(gzip -c  > {filtered_R1}) \
            --fastqout_rev >(gzip -c > {filtered_R2}) --fastq_maxee {maxEE} > {log_vsearch} 2>&1
        """)

def pseudo_quality_trim(input_R1, input_R2, filtered_R1, filtered_R2):
    return shell(f"ln -s {input_R1} {filtered_R1} && ln -s {input_R2} {filtered_R2}")

rule trim_before_haplotype_generation:
    threads: 3 if quality_trim_required else 1
    localrule: False if quality_trim_required else True
    input:
        R1 = f"{outdir}/samples/{{sample}}/mhaps-reads/primer-trimmed_R1.fastq.gz",
        R2 = f"{outdir}/samples/{{sample}}/mhaps-reads/primer-trimmed_R2.fastq.gz",
    output:
        R_trimmed_R1 = temp(f"{outdir}/samples/{{sample}}/mhaps-reads/T-primer-trimmed_R1.fastq.gz"),
        R_trimmed_R2 = temp(f"{outdir}/samples/{{sample}}/mhaps-reads/T-primer-trimmed_R2.fastq.gz"),
        Filtered_final_R1 = f"{outdir}/samples/{{sample}}/mhaps-reads/trimmed-filtered_R1.fastq.gz",
        Filtered_final_R2 = f"{outdir}/samples/{{sample}}/mhaps-reads/trimmed-filtered_R2.fastq.gz",
    params:
        trim_right_R1 = config.get("trim_right_R1", 5),
        trim_right_R2 = config.get("trim_right_R2", 5),
        cut_window = config.get("cut_window", 1),
        window_minqual = config.get("truncQ", 20),
        min_length = config.get("min_length", 30),
        maxEE = config.get("maxEE", 5)
    log:
        fastp_json = f"{outdir}/samples/{{sample}}/logs/fastp_tr_before_read_merging.json",
        fastp_html = f"{outdir}/samples/{{sample}}/logs/fastp_tr_before_read_merging.html",
        vsearch_filter_log = f"{outdir}/samples/{{sample}}/logs/fastx_mee_after_qtrim.log",
    run:
        if quality_trim_required:
            quality_trim(input.R1, input.R2, output.R_trimmed_R1, output.R_trimmed_R2,
                output.Filtered_final_R1, output.Filtered_final_R2, threads, log.vsearch_filter_log,
                log.fastp_html, log.fastp_json, params)
        else:
            pseudo_quality_trim(input.R1, input.R2, output.Filtered_final_R1, output.Filtered_final_R2)
            shell(f"touch {output.R_trimmed_R1} {output.R_trimmed_R2}")

rule per_marker_trim_before_haplotype_generation:
    threads: 3 if quality_trim_required else 1
    localrule: False if quality_trim_required else True
    input:
        trim_completion = f"{outdir}/samples/{{sample}}/markers-reads/.trimmed"
    output:
        qtrim_completion_flag = f"{outdir}/samples/{{sample}}/markers-reads/.quality_trimmed"
    params:
        trim_right_R1 = config.get("trim_right_R1", 5),
        trim_right_R2 = config.get("trim_right_R2", 5),
        cut_window = config.get("cut_window", 1),
        window_minqual = config.get("truncQ", 20),
        min_length = config.get("min_length", 30),
        maxEE = config.get("maxEE", 5)
    log:
        fastp_json = f"{outdir}/samples/{{sample}}/logs/fastp_tr_before_read_merging.json",
        fastp_html = f"{outdir}/samples/{{sample}}/logs/fastp_tr_before_read_merging.html",
        vsearch_filter_log = f"{outdir}/samples/{{sample}}/logs/fastx_mee_after_qtrim.log",
    params:
        log_path = f"{outdir}/samples/{{sample}}/logs/"
    run:
        import os

        marker_read_dir = os.path.dirname(input.trim_completion)
        all_input_files = [(os.path.join(marker_read_dir, marker, "trimmed_R1.fastq.gz"),
            os.path.join(marker_read_dir, marker, "trimmed_R2.fastq.gz")) for marker in Markers]
        all_trimmed_files = [(os.path.join(marker_read_dir, marker, "T-primer-trimmed_R1.fastq.gz"),
            os.path.join(marker_read_dir, marker, "T-primer-trimmed_R2.fastq.gz")) for marker in Markers]
        all_trimmed_filtered_files = [(os.path.join(marker_read_dir, marker, "trimmed-filtered_R1.fastq.gz"),
            os.path.join(marker_read_dir, marker, "trimmed-filtered_R2.fastq.gz")) for marker in Markers]
        all_log_files = [(os.path.join(params.log_path, marker, "fastp_tr_before_read_merging.json"),
            os.path.join(params.log_path, marker, "fastp_tr_before_read_merging.html"),
            os.path.join(params.log_path, marker, "fastx_mee_after_qtrim.log")) for marker in Markers]

        if quality_trim_required:
            for (inR1, inR2), (outTR1, outTR2), (outR1, outR2), (log_json, log_html, log_vsearch) in zip(
                all_input_files, all_trimmed_files, all_trimmed_filtered_files, all_log_files):
                quality_trim(inR1, inR2, outTR1, outTR2, outR1, outR2, threads, log_vsearch, log_html, log_json, params)
                shell("rm {outTR1} {outTR2}")
        else:
            for (inR1, inR2), (outR1, outR2) in zip(all_input_files, all_trimmed_filtered_files):
                pseudo_quality_trim(inR1, inR2, outR1, outR2)