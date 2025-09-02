dada2_bbmapmerge_bbmap = config.get('merge_map') # dada2, bbmap_merge or bbmerge

rule extract_sequences_from_seqtab:
    input:
        f"{outdir}/malamp/{dada2_bbmapmerge_bbmap}/seqtab.tsv"
    output:
        fasta = f"{outdir}/malamp/haplotypes.fasta",
    shell:
        "python {microhaps_basedir}/scripts/seqtab_to_fasta.py --table {input[0]} --output_fasta {output.fasta}"

rule align_haplotypes_to_reference:
    threads: 1
    input:
        fasta = f"{outdir}/malamp/haplotypes.fasta",
    output:
        paf = f"{outdir}/malamp/haplotypes.paf",
    params:
        reference = insertseq_file,
        minimap_params = config.get('minimap2_params', ''),
        cs_style = "long" if config.get("post_process", "cs_long") == "cs_long" else "short"
    shell:
        "minimap2 -x sr -t {threads} {params.minimap_params} --secondary=no --cs={params.cs_style} {params.reference} {input.fasta} --paf-no-hit  -o {output.paf}"

rule post_process_haplotype:
    threads: 1
    input:
        fasta = f"{outdir}/malamp/haplotypes.fasta",
        paf = f"{outdir}/malamp/haplotypes.paf",
        seqtab = f"{outdir}/malamp/{dada2_bbmapmerge_bbmap}/seqtab.tsv",
    output:
        output_table = f"{outdir}/malamp/outputHaplotypes.tsv",
    shell:
        "python {microhaps_basedir}/scripts/post_process_dada2.py"
        " --paf {input.paf}"
        " --fasta {input.fasta}"
        " --seqtab {input.seqtab}"
        " --output {output.output_table}"
        " --insert {insertseq_file}"

rule post_process_CIGAR:
    threads: 4
    input:
        f"{outdir}/malamp/{dada2_bbmapmerge_bbmap}/seqtab.tsv"
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
        seqtab = f"{outdir}/malamp/{dada2_bbmapmerge_bbmap}/seqtab.tsv"
    output:
        cigar = f"{outdir}/malamp/outputCIGAR.tsv",
        asv_to = f"{outdir}/malamp/asv_to_cigar"
    shell:
        "python {microhaps_basedir}/scripts/ASV_to_CIGAR.py"
        " {input.Seqs} {input.Table} {input.seqtab} {output.cigar}"
        " --asv_to_cigar {output.asv_to}"
        " -a {outdir}/alignments"
        " --amp_db {insertseq_file}"


