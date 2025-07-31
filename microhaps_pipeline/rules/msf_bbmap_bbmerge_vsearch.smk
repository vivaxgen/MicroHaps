
bbmap_or_bbmerge = config.get('merge_map') # should be either "bbmap" or "bbmerge"

rule index_ref:
    input:
        ref = refseq_file,
    output:
        directory(f"{outdir}/ref"),
    shell:
        f"bbmap.sh ref={{input.ref}} path={outdir}"

# reference: https://www.seqanswers.com/forum/bioinformatics/bioinformatics-aa/40046-non-overlapping-read-assembly?t=45477
rule bbmap_to_fq:
    threads: config.get('bbmap_merge_threads', 4)
    input:
        R1 = f"{outdir}/samples/{{sample}}/mhaps-reads/primer-trimmed_R1.fastq.gz",
        R2 = f"{outdir}/samples/{{sample}}/mhaps-reads/primer-trimmed_R2.fastq.gz",
        ref = f"{outdir}/ref",
    output:
        mapped = f"{outdir}/samples/{{sample}}/bbmap_merge/bbmapped.fastq.gz",
    log:
        f"{outdir}/samples/{{sample}}/logs/bbmap_merge_map.log",
    params:
        outdir = outdir,
        bbmap_parameters = config.get('bbmap_parameters', ''),
    shell:
        "bbmap.sh in1={input.R1} in2={input.R2} outm={output.mapped} t={threads} path={params.outdir} po rbm don {params.bbmap_parameters} > {log} 2>&1"
        # po = paired only, "rbm" and "don" (renamebymapping and deleteoldname)

rule bbmap_merge_mapped:
    threads: config.get('bbmap_merge_threads', 4)
    input:
        merged = f"{outdir}/samples/{{sample}}/bbmap_merge/bbmapped.fastq.gz",
    output:
        merged = f"{outdir}/samples/{{sample}}/bbmap_merge/final_merged.fastq.gz",
    log:
        f"{outdir}/samples/{{sample}}/logs/bbmap_merge_merge.log",
    params:
        bbmerge_parameters = config.get('bbmerge_parameters', ''),
    shell:
        "bbmerge.sh in={input.merged} out={output.merged} t={threads} int usemapping parsecustom {params.bbmerge_parameters} > {log} 2>&1"

rule bbmerge_to_fq:
    threads: config.get('bbmap_merge_threads', 4)
    input:
        R1 = f"{outdir}/samples/{{sample}}/mhaps-reads/primer-trimmed_R1.fastq.gz",
        R2 = f"{outdir}/samples/{{sample}}/mhaps-reads/primer-trimmed_R2.fastq.gz",
    output:
        merged = f"{outdir}/samples/{{sample}}/bbmerge/final_merged.fastq.gz",
    log:
        f"{outdir}/samples/{{sample}}/logs/bbmerge.log",
    params:
        outdir = outdir,
        bbmerge_parameters = config.get('bbmerge_parameters', ''),
    shell:
        "bbmerge.sh in1={input.R1} in2={input.R2} outm={output.merged} t={threads} path={params.outdir} {params.bbmerge_parameters} > {log} 2>&1"
        # po = paired only, "rbm" and "don" (renamebymapping and deleteoldname)


# Denoising with vsearch
# reference: https://learnmetabarcoding.github.io/LearnMetabarcoding/filtering/freq_filtering_and_denoising.html
rule indv_bbfastq_to_fasta:
    input:
        merged = f"{outdir}/samples/{{sample}}/{bbmap_or_bbmerge}/final_merged.fastq.gz",
    output:
        fasta = f"{outdir}/samples/{{sample}}/{bbmap_or_bbmerge}/final_merged.fasta",
    log:
        f"{outdir}/samples/{{sample}}/logs/fastx_uniques.log",
    shell:
        "vsearch --fastx_uniques {input.merged} --fastaout {output.fasta} --fasta_width 0 --sizeout --relabel {wildcards.sample}: > {log} 2>&1"

rule indv_fasta_to_fasta:
    input:
        fastas = expand(f"{outdir}/samples/{{sample}}/{bbmap_or_bbmerge}/final_merged.fasta", sample=IDs),
    output:
        tempfasta = temp(f"{outdir}/malamp/{bbmap_or_bbmerge}/combined.fasta"),
        fasta = f"{outdir}/malamp/{bbmap_or_bbmerge}/dereplicated_counted.fasta",
    log:
        f"{outdir}/logs/fastx_uniques.log",
    shell:
        "cat {input.fastas} > {output.tempfasta} && "
        "vsearch --fastx_uniques {output.tempfasta} --fastaout {output.fasta} --fasta_width 0 --sizein --sizeout --relabel_sha1 > {log} 2>&1"


# if the research questions are tolerant to the loss of rare real sequences and/or
# it is important that as high a proportion of ASVs be valid as possible,
# then stricter settings (higher --minsize, lower --unoise_alpha) would be sensible. 
# alpha - determines the threshold level of dissimilarity between frequent and infrequent reads
# for exclusion of infrequent reads.
rule cluster_fasta:
    threads: config.get('bbmap_merge_threads', 4)
    input:
        fasta = f"{outdir}/malamp/{bbmap_or_bbmerge}/dereplicated_counted.fasta",
    output:
        fasta = f"{outdir}/malamp/{bbmap_or_bbmerge}/dereplicated_counted_centroids.fasta",
    log:
        f"{outdir}/logs/vsearch_cluster.log",
    shell:
        "vsearch --cluster_unoise {input.fasta} --minsize 4 --fasta_width 0 --unoise_alpha 2 --centroids {output.fasta} > {log} 2>&1"

# chimeric removal
rule chimera_filtering:
    threads: config.get('vsearch_threads', 4)
    input:
        fasta = f"{outdir}/malamp/{bbmap_or_bbmerge}/dereplicated_counted_centroids.fasta",
    output:
        kept = f"{outdir}/malamp/{bbmap_or_bbmerge}/dereplicated_counted_no_chimeras.fasta",
        discarded = f"{outdir}/malamp/{bbmap_or_bbmerge}/dereplicated_counted_chimeras.fasta",
    log:
        f"{outdir}/logs/vsearch_chimera.log",
    shell:
        "vsearch --uchime3_denovo {input.fasta} --nonchimeras {output.kept} --chimeras {output.discarded} --fasta_width 0 --threads {threads} > {log} 2>&1"

rule indv_otutab:
    threads: config.get('vsearch_threads', 4)
    input:
        fasta = f"{outdir}/malamp/{bbmap_or_bbmerge}/dereplicated_counted_no_chimeras.fasta",
        sample_fa = f"{outdir}/samples/{{sample}}/{bbmap_or_bbmerge}/final_merged.fasta",
    output:
        otutab = f"{outdir}/samples/{{sample}}/{bbmap_or_bbmerge}/seqtab.tsv",
    shell:
        "vsearch --search_exact {input.sample_fa} --db {input.fasta} --otutabout {output.otutab} --threads {threads}"

rule merge_otutab:
    input:
        otutabs = expand(f"{outdir}/samples/{{sample}}/{bbmap_or_bbmerge}/seqtab.tsv", sample=IDs),
        filtered_fasta = f"{outdir}/malamp/{bbmap_or_bbmerge}/dereplicated_counted_no_chimeras.fasta",
    output:
        otutab = f"{outdir}/malamp/{bbmap_or_bbmerge}/seqtab.tsv",
    shell:
        "python {microhaps_basedir}/scripts/merge_otutab.py"
        " --denoised_fasta {input.filtered_fasta}"
        " --out_seqtab {output.otutab}"
        " --indv_seqtabs {input.otutabs}"

