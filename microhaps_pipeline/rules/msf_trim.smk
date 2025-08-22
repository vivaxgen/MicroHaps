new_postprocess = config.get('post_process', "old")
primers_trimmed = config.get("primers_trimmed", False)

if not primers_trimmed:
    # trim primers from reads per index sample
    rule trim:
        input:
            R1 = f"{outdir}/samples/{{sample}}/mhaps-reads/target_R1.fastq.gz",
            R2 = f"{outdir}/samples/{{sample}}/mhaps-reads/target_R2.fastq.gz",
            prim_fw = primer_fw_file,
            prim_rv = primer_rev_file,
        output:
            R1 = f"{outdir}/samples/{{sample}}/mhaps-reads/primer-trimmed_R1.fastq.gz",
            R2 = f"{outdir}/samples/{{sample}}/mhaps-reads/primer-trimmed_R2.fastq.gz",
        log:
            f"{outdir}/samples/{{sample}}/logs/trim_before_dada2.log",
        params:
            platform = config['platform'],
            trim_qv = config['trimqv'],
            additional_params = f"--nextseq-trim={config.get('trimqv', 20)}" if is_nextseq_or_novaseq() else ""
        shell: 
            "cutadapt -g file:{input.prim_fw} -G file:{input.prim_rv} -o {output.R1} -p {output.R2}"
            " --pair-adapters --discard-untrimmed {params.additional_params} --action=trim {input.R1} {input.R2} > {log} 2>&1"
else:
    rule trim:
        input:
            R1 = f"{outdir}/samples/{{sample}}/mhaps-reads/target_R1.fastq.gz",
            R2 = f"{outdir}/samples/{{sample}}/mhaps-reads/target_R2.fastq.gz",
        output:
            R1 = f"{outdir}/samples/{{sample}}/mhaps-reads/primer-trimmed_R1.fastq.gz",
            R2 = f"{outdir}/samples/{{sample}}/mhaps-reads/primer-trimmed_R2.fastq.gz",
        shell:
            "ln -s {input.R1} {output.R1} && ln -s {input.R2} {output.R2}"

rule optical_dedup:
    input:
        R1 = f"{outdir}/samples/{{sample}}/mhaps-reads/primer-trimmed_R1.fastq.gz",
        R2 = f"{outdir}/samples/{{sample}}/mhaps-reads/primer-trimmed_R2.fastq.gz",
    output:
        R1 = f"{outdir}/samples/{{sample}}/mhaps-reads/Odedup-primer-trimmed_R1.fastq.gz",
        R2 = f"{outdir}/samples/{{sample}}/mhaps-reads/Odedup-primer-trimmed_R2.fastq.gz",
    log:
        clump_log = f"{outdir}/samples/{{sample}}/logs/clumpify.log",
    params:
        instrument_specific = "spany adjacent" if config.get("instrument", "generic") == "nextseq" else "",
    shell:
        """
        clumpify.sh in1={input.R1} in2={input.R2} out1={output.R1} out2={output.R2} dedupe optical {params.instrument_specific} > {log.clump_log} 2>&1
        """

rule map_dedup_trimmed_seq:
    threads: 2
    input:
        R1 = f"{outdir}/samples/{{sample}}/mhaps-reads/Odedup-primer-trimmed_R1.fastq.gz",
        R2 = f"{outdir}/samples/{{sample}}/mhaps-reads/Odedup-primer-trimmed_R2.fastq.gz",
    output:
        ori_bam =  temp(f"{outdir}/samples/{{sample}}/mhaps-reads/mapped.bam"),
        bam = f"{outdir}/samples/{{sample}}/mhaps-reads/mapped.cigar_fixed.bam"
    params:
        rg = lambda w: f"-R '@RG\tID:{w.sample}\tSM:{w.sample}\tLB:LIB-{w.sample}\tPL:generic'",
    log:
        map_dedup = f"{outdir}/samples/{{sample}}/logs/map_dedup.log",
        sort_dedup = f"{outdir}/samples/{{sample}}/logs/sort_dedup.log",
        reformat = f"{outdir}/samples/{{sample}}/logs/reformat.log"
    shell:
        "bwa-mem2 mem -M -t {threads} {params.rg} {refseq} {input.R1} {input.R2} 2> {log.map_dedup} | samtools sort -o {output.ori_bam} 2> {log.sort_dedup} && "
        "reformat.sh sam=1.4 in={output.ori_bam} out={output.bam} > {log.reformat} 2>&1"

### This improve bbmerge <--- potential 
if config.get('recal_qual', False) == True and config.get('merge_map', 'dada2') != 'dada2': 
    rule calctruequality: 
        threads: 4
        input:
            bams = expand(f"{outdir}/samples/{{sample}}/mhaps-reads/mapped.cigar_fixed.bam", sample=IDs),
        output:
            expected_dir = directory(f"{outdir}/ref/qual"),
            done_signal = f"{outdir}/ref/qual/.done"
        log:
            f"{outdir}/logs/calctruequality.log"
        params:
            outdir = outdir,
            recal_param = config.get("recal_param", "")
        run:
            shell(f"calctruequality.sh in={','.join(input.bams)} passes=1 path={params.outdir} t={threads} {params.recal_param} > {log} 2>&1")
            shell("touch {output.done_signal}")

    rule recalibrate_qual:
        threads: 1
        input:
            R1 = f"{outdir}/samples/{{sample}}/mhaps-reads/Odedup-primer-trimmed_R1.fastq.gz",
            R2 = f"{outdir}/samples/{{sample}}/mhaps-reads/Odedup-primer-trimmed_R2.fastq.gz",
            qual = f"{outdir}/ref/qual/.done",
        output:
            recal_non_standard_1 = temp(f"{outdir}/samples/{{sample}}/mhaps-reads/recalibrated_R1.temp.fastq.gz"),
            recal_non_standard_2 = temp(f"{outdir}/samples/{{sample}}/mhaps-reads/recalibrated_R2.temp.fastq.gz"),
            R1 = f"{outdir}/samples/{{sample}}/mhaps-reads/recalibrated_R1.fastq.gz",
            R2 = f"{outdir}/samples/{{sample}}/mhaps-reads/recalibrated_R2.fastq.gz",
        params:
            outdir = outdir,
            recal_param = config.get("recal_param", "")
        shell:
            "bbduk.sh in1={input.R1} in2={input.R2} out1={output.recal_non_standard_1} out2={output.recal_non_standard_2} t={threads} {params.recal_param} recalibrate path={params.outdir} && "
            "reformat.sh in1={output.recal_non_standard_1} in2={output.recal_non_standard_2} out1={output.R1} out2={output.R2} t={threads} maxcalledquality=41"
else:
    rule pseudo_recalibrate_qual:
        localrule: True
        input:
            R1 = f"{outdir}/samples/{{sample}}/mhaps-reads/Odedup-primer-trimmed_R1.fastq.gz",
            R2 = f"{outdir}/samples/{{sample}}/mhaps-reads/Odedup-primer-trimmed_R2.fastq.gz",
        output:
            R1 = f"{outdir}/samples/{{sample}}/mhaps-reads/recalibrated_R1.fastq.gz",
            R2 = f"{outdir}/samples/{{sample}}/mhaps-reads/recalibrated_R2.fastq.gz",
        shell:
            "ln -s {input.R1} {output.R1} && ln -s {input.R2} {output.R2}"

if config.get('additional_trim_before_haplotype_generation', False) == True or config.get('merge_map', 'dada2') == 'dada2':
    rule trim_before_haplotype_generation:
        threads: 3
        input:
            R1 = f"{outdir}/samples/{{sample}}/mhaps-reads/recalibrated_R1.fastq.gz",
            R2 = f"{outdir}/samples/{{sample}}/mhaps-reads/recalibrated_R2.fastq.gz",
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
            fastp_json = f"{outdir}/samples/{{sample}}/logs/fastp_tr_before_hapgen.json",
            fastp_html = f"{outdir}/samples/{{sample}}/logs/fastp_tr_before_hapgen.html",
            vsearch_filter_log = f"{outdir}/samples/{{sample}}/logs/fastx_mee.log",
        shell:
            """
            fastp -A -Q -w {threads} -i {input.R1} -I {input.R2} -o {output.R_trimmed_R1} -O {output.R_trimmed_R2} \
                --trim_tail1 {params.trim_right_R1} --trim_tail2 {params.trim_right_R2} \
                --cut_tail --cut_tail_window_size {params.cut_window} --cut_tail_mean_quality {params.window_minqual} \
                -l {params.min_length} \
                -j {log.fastp_json} -h {log.fastp_html} &&
            vsearch --fastx_filter {output.R_trimmed_R1} --reverse {output.R_trimmed_R2} --fastqout >(gzip -c  > {output.Filtered_final_R1}) \
                --fastqout_rev >(gzip -c > {output.Filtered_final_R2}) --fastq_maxee {params.maxEE} > {log.vsearch_filter_log} 2>&1
            """
else:
    rule pseudo_trim_before_haplotype_generation:
        localrule: True
        input:
            R1 = f"{outdir}/samples/{{sample}}/mhaps-reads/recalibrated_R1.fastq.gz",
            R2 = f"{outdir}/samples/{{sample}}/mhaps-reads/recalibrated_R2.fastq.gz",
        output:
            Filtered_final_R1 = f"{outdir}/samples/{{sample}}/mhaps-reads/trimmed-filtered_R1.fastq.gz",
            Filtered_final_R2 = f"{outdir}/samples/{{sample}}/mhaps-reads/trimmed-filtered_R2.fastq.gz",
        shell:
            "ln -s {input.R1} {output.Filtered_final_R1} && ln -s {input.R2} {output.Filtered_final_R2}"