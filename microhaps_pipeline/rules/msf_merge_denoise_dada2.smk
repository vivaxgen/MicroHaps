
rule create_meta:
    # create a list
    localrule: True
    input:
        expand(f'{outdir}/samples/{{sample}}/mhaps-reads/trimmed-filtered_R1.fastq.gz', sample=IDs)
    output:
        f"{outdir}/malamp/meta"
    log:
        f"{outdir}/mylog.txt"
    params:
        samples = IDs,
        outdir = outdir
    run:
        import pandas as pd

        read1_list = []
        read2_list = []
        for sample in params.samples:
            read1_list.append(f"{params.outdir}/samples/{sample}/mhaps-reads/trimmed-filtered_R1.fastq.gz")
            read2_list.append(f"{params.outdir}/samples/{sample}/mhaps-reads/trimmed-filtered_R2.fastq.gz")

        pd.DataFrame({
            'id': params.samples,
            'ip1': read1_list, 'ip2': read2_list
        }).to_csv(output[0], sep='\t', index=False, header=False)

rule run_dada2R:
    threads: 16
    input:
        meta = f"{outdir}/malamp/meta",
    output:
        f"{outdir}/malamp/dada2/seqtab.tsv"
    params:
        dir = f"{outdir}/malamp/dada2",
        output_filename = "seqtab.tsv",
        ####### read-filtering performed prior #######
        maxEE = f"{config.get("maxEE", 5)},{config.get("maxEE", 5)}",          # Can recheck but should not perform additional filter, since already pre-filtered previously
        trimRight = "0,0",                                                          # Already trimmed, should not redo
        minLen = config['min_length'], # Can recheck but should not perform additional filter, since already pre-filtered previously
        truncQ = f"{config.get("truncQ", 20)},{config.get("truncQ", 20)}",   # Can recheck but should not perform additional trim, since already pre-trimmed previously
        ##############################################
        max_consist = config['max_consist'],
        omega_a = config['omegaA'],
        justConcatenate = config['justconcat'],
        _class = config['class'],

    shell:
        """
        Rscript {microhaps_basedir}/scripts/runDADA2.R  \
            --path_to_meta {input.meta} \
            --class {params._class} \
            --dir {params.dir} \
            --output_filename {params.output_filename} \
            --maxEE {params.maxEE} \
            --trimRight {params.trimRight} \
            --minLen {params.minLen} \
            --truncQ {params.truncQ} \
            --max_consist {params.max_consist} \
            --omega_a {params.omega_a} \
            --justConcatenate {params.justConcatenate} \
            --threads {threads} 
        """



# EOF
