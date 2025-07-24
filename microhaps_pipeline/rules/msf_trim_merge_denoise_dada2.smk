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
            --threads {threads} 
        """

rule extract_sequences_from_dada2_seqtab:
    input:
        f"{outdir}/malamp/dada2/seqtab.tsv"
    output:
        fasta = f"{outdir}/malamp/haplotypes.fasta",
    shell:
        "python {microhaps_basedir}/scripts/seqtab_to_fasta.py --table {input[0]} --output_fasta {output.fasta}"

rule align_haplotypes_to_reference:
    threads: 3
    input:
        fasta = f"{outdir}/malamp/haplotypes.fasta",
    output:
        paf = f"{outdir}/malamp/haplotypes.paf",
    params:
        reference = insertseq_file,
        cs_style = "long" if new_postprocess == "cs_long" else "short"
    shell:
        "minimap2 -x sr -t {threads} --secondary=no --cs={params.cs_style} {params.reference} {input.fasta} --paf-no-hit  -o {output.paf}"

rule post_process_dada2:
    threads: 1
    input:
        fasta = f"{outdir}/malamp/haplotypes.fasta",
        paf = f"{outdir}/malamp/haplotypes.paf",
        seqtab = f"{outdir}/malamp/dada2/seqtab.tsv",
    output:
        output_table = f"{outdir}/malamp/outputHaplotypes.tsv",
    shell:
        "python {microhaps_basedir}/scripts/post_process_dada2.py"
        " --paf {input.paf}"
        " --fasta {input.fasta}"
        " --seqtab {input.seqtab}"
        " --output {output.output_table}"

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
        tab = f"{outdir}/malamp/outputHaplotypes.tsv" if new_postprocess != "old" else f"{outdir}/malamp/outputCIGAR.tsv",
    output:
        depths = f"{outdir}/malamp/depths.tsv"
    params:
        mindepth = 10,
        outdir = lambda w, output: output.depths.rsplit('/', 1)[0]
    shell:
        "ngs-pl tab-to-QC -d {params.mindepth} --outdir {params.outdir} {input}"


rule qc_heatmap:
    input:
        depth = "{pfx}/malamp/depths.tsv",
    output:
        depth = "{pfx}/malamp/depths-microhaps.png",
        marker = "{pfx}/malamp/markers-microhaps.png",
    shell:
        "ngs-pl tab-to-plots --index-column Amplicon_name --additional-title 'Microhaps'"
        "  --outheatmap {output.depth} --outmarker {output.marker} {input.depth}"


rule depth_ratio:
    localrule: True
    input:
        depth_mapped = "{pfx}/depths-mapped.tsv",
        depth_microhaps = "{pfx}/malamp/depths.tsv",
    output:
        ratio = "{pfx}/malamp/depth-ratio.tsv"
    run:

        import pandas as pd

        # Read the depth files
        mapped_df = pd.read_csv(input.depth_mapped, sep='\t')
        microhaps_df = pd.read_csv(input.depth_microhaps, sep='\t')

        depth_mapped = mapped_df.iloc[:, 4:].set_index(mapped_df.Amplicon_name)
        depth_microhaps = microhaps_df.iloc[:, 4:].set_index(microhaps_df.Amplicon_name)

        ratio_df = (depth_microhaps / depth_mapped).fillna(0).reindex(depth_mapped.index).reset_index()

        ratio_df.to_csv(output.ratio, sep='\t', index=False)


rule depth_ratio_heatmap:
    localrule: True
    input:
        ratio = "{pfx}/depth-ratio.tsv",
    output:
        heatmap = "{pfx}/depth-ratio-heatmap.png",
        marker = "{pfx}/depth-ratio-markers.png",
    shell:
        "ngs-pl tab-to-plots --value-is-ratio --index-column Amplicon_name --additional-title 'Microhaps Ratio'"
        "  --outheatmap {output.heatmap} --outmarker {output.marker} {input.ratio}"


# EOF
