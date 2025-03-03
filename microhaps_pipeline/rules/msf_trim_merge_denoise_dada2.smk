
# trim primers from reads per index sample
rule trim:
    input:
        #unpack(read_files.get_read_file_as_dict),
        R1 = f"{outdir}/samples/{{sample}}/mhaps-reads/target_R1.fastq.gz",
        R2 = f"{outdir}/samples/{{sample}}/mhaps-reads/target_R2.fastq.gz",
        prim_fw = primer_fw_file,
        prim_rv = primer_rev_file,
    output:
        #f"{outdir}/{{sample}}/reads/raw-{{idx}}_R1.fastq.gz",
        #f"{outdir}/{{sample}}/reads/raw-{{idx}}_R2.fastq.gz"
        R1 = f"{outdir}/samples/{{sample}}/mhaps-reads/primer-trimmed_R1.fastq.gz",
        R2 = f"{outdir}/samples/{{sample}}/mhaps-reads/primer-trimmed_R2.fastq.gz",
    params:
        platform = config['platform'],
        trim_qv = config['trimqv'],
        additional_params = f"--nextseq-trim={config.get('trimqv', 20)}" if is_nextseq_or_novaseq() else ""
    shell: 
        "cutadapt -g file:{input.prim_fw} -G file:{input.prim_rv} -o {output.R1} -p {output.R2}"
        " --pair-adapters --discard-untrimmed {params.additional_params} --action=trim {input.R1} {input.R2}"


rule create_meta:
    # create a list
    localrule: True
    input:
        expand(f'{outdir}/samples/{{sample}}/mhaps-reads/primer-trimmed_R1.fastq.gz', sample=IDs)
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
            read1_list.append(f"{params.outdir}/samples/{sample}/mhaps-reads/primer-trimmed_R1.fastq.gz")
            read2_list.append(f"{params.outdir}/samples/{sample}/mhaps-reads/primer-trimmed_R2.fastq.gz")

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
        maxEE = config['maxEE'],
        trimRight = config['trim_right'],
        minLen = config['min_length'],
        truncQ = config['truncQ'],
        max_consist = config['max_consist'],
        omega_a = config['omegaA'],
        justConcatenate = config['justconcat'],
        _class = config['class']

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
        """


rule post_process:
    threads: 4
    input:
        f"{outdir}/malamp/dada2/seqtab.tsv"
    output:
        Table = f"{outdir}/malamp/ASVTable.txt",
        Seqs = f"{outdir}/malamp/ASVSeqs.fasta"
    shell:
        "Rscript {microhaps_basedir}/scripts/postProc_dada2.R"
        " -s {input}"
        " --strain PvP01"
        " -ref {insertseq_file}"
        " -o {output.Table}"
        " --fasta --parallel"


rule asv_to_cigar:
    threads: 4
    input:
        Table = f"{outdir}/malamp/ASVTable.txt",
        Seqs = f"{outdir}/malamp/ASVSeqs.fasta",
        seqtab = f"{outdir}/malamp/dada2/seqtab.tsv"
    output:
        cigar = f"{outdir}/malamp/outputCIGAR.tsv",
        asv_to = f"{outdir}/malamp/asv_to_cigar"
    shell:
        "python {microhaps_basedir}/scripts/ASV_to_CIGAR.py"
        " {input.Seqs} {input.Table} {input.seqtab} {output.cigar}"
        " --asv_to_cigar {output.asv_to}"
        " -a {outdir}/alignments"
        " --amp_db {insertseq_file}"


rule qc_outputCIGAR:
    localrule: True
    input:
        tab = f"{outdir}/malamp/outputCIGAR.tsv"
    output:
        depths = f"{outdir}/malamp/depths.tsv"
    params:
        mindepth = 10,
        outdir = lambda w, output: output.depths.rsplit('/', 1)[0]
    shell:
        "ngs-pl tab-to-QC -d {params.mindepth} --outdir {params.outdir} {input}"


# EOF
