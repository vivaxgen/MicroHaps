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
