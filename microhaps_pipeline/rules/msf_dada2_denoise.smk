rule merge_fastp_merged_markers_to_samples:
    input:
        f"{outdir}/samples/{{sample}}/markers-reads/.merged_denoised"
    output:
        f"{outdir}/samples/{{sample}}/markers-reads/merged_filtered.fastq.gz"
    params:
        merged_markers = lambda wildcards: " ".join([
            f"{outdir}/samples/{wildcards.sample}/markers-reads/{marker}/merged_filtered.fastq.gz" for marker in Markers
        ]),
    shell:
        "cat {params.merged_markers} > {output}"

rule create_denoise_meta:
    localrule: True
    input:
        expand(f'{outdir}/samples/{{sample}}/markers-reads/merged_filtered.fastq.gz', sample=IDs)
    output:
        f"{outdir}/malamp/meta_denoise"
    log:
        f"{outdir}/mylog_denoise.txt"
    params:
        samples = IDs,
        outdir = outdir
    run:
        import pandas as pd

        read_list = []
        for sample in params.samples:
            read_list.append(f"{params.outdir}/samples/{sample}/markers-reads/merged_filtered.fastq.gz")

        pd.DataFrame({
            'id': params.samples,
            'ip': read_list
        }).to_csv(output[0], sep='\t', index=False, header=False)

rule run_dada2R_denoise:
    threads: 16
    input:
        meta = f"{outdir}/malamp/meta_denoise",
    output:
        f"{outdir}/malamp/fastp_dada2/seqtab.tsv"
    params:
        dir = f"{outdir}/malamp/fastp_dada2",
        output_filename = "seqtab.tsv",
        ####### read-filtering performed prior #######
        maxEE = {config.get("maxEE", 5)},          # Can recheck but should not perform additional filter, since already pre-filtered previously
        trimRight = 0,                                                          # Already trimmed, should not redo
        minLen = config['min_length'], # Can recheck but should not perform additional filter, since already pre-filtered previously
        truncQ = {config.get("dada2_truncQ_first", 5)}, 
        ##############################################
        max_consist = config['max_consist'],
        omega_a = config['omegaA'],
    shell:
        """
        Rscript {microhaps_basedir}/scripts/runDADA2_denoise_only.R  \
            --path_to_meta {input.meta} \
            --dir {params.dir} \
            --output_filename {params.output_filename} \
            --maxEE {params.maxEE} \
            --trimRight {params.trimRight} \
            --minLen {params.minLen} \
            --truncQ {params.truncQ} \
            --max_consist {params.max_consist} \
            --omega_a {params.omega_a}
        """
    