bbmap_or_bbmerge_fastq_merge = config.get('merge_map') # should be either "bbmap" or "bbmerge" or "fastq_merge"
markers = [line.replace(">", "") for line in open(insertseq_file, 'r').read().split() if line.startswith('>')]

ruleorder:
    indv_bbfastq_to_fasta > bunzip2

rule index_ref:
    input:
        ref = refseq_file,
    output:
        directory(f"{outdir}/ref"),
    shell:
        f"bbmap.sh ref={{input.ref}} path={outdir}"

rule trim_R:
    input:
        R1 = f"{outdir}/samples/{{sample}}/mhaps-reads/primer-trimmed_R1.fastq.gz",
        R2 = f"{outdir}/samples/{{sample}}/mhaps-reads/primer-trimmed_R2.fastq.gz",
    output:
        R1 = temp(f"{outdir}/samples/{{sample}}/mhaps-reads/R-primer-trimmed_R1.fastq.gz"),
        R2 = temp(f"{outdir}/samples/{{sample}}/mhaps-reads/R-primer-trimmed_R2.fastq.gz"),
    shell:
        """
        vsearch --fastx_filter {input.R1} --reverse {input.R2} --fastqout >(gzip -c  > {output.R1}) \
        --fastqout_rev >(gzip -c > {output.R2}) --fastq_stripright 10 --fastq_maxee 1
        """

rule optical_dedup:
    input:
        R1 = f"{outdir}/samples/{{sample}}/mhaps-reads/R-primer-trimmed_R1.fastq.gz",
        R2 = f"{outdir}/samples/{{sample}}/mhaps-reads/R-primer-trimmed_R2.fastq.gz",
    output:
        R1 = temp(f"{outdir}/samples/{{sample}}/mhaps-reads/TR-primer-trimmed_R1.fastq.gz"),
        R2 = temp(f"{outdir}/samples/{{sample}}/mhaps-reads/TR-primer-trimmed_R2.fastq.gz"),
    params:
        instrument_specific = "spany adjacent" if config.get("instrument", "generic") == "nextseq" else "",
    shell:
        """
        clumpify.sh in1={input.R1} in2={input.R2} out1={output.R1} out2={output.R2} dedupe optical {params.instrument_specific}
        """

# reference: https://www.seqanswers.com/forum/bioinformatics/bioinformatics-aa/40046-non-overlapping-read-assembly?t=45477
rule bbmap_to_fq:
    threads: config.get('bbmap_merge_threads', 4)
    input:
        R1 = f"{outdir}/samples/{{sample}}/mhaps-reads/TR-primer-trimmed_R1.fastq.gz",
        R2 = f"{outdir}/samples/{{sample}}/mhaps-reads/TR-primer-trimmed_R2.fastq.gz",
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
        R1 = f"{outdir}/samples/{{sample}}/mhaps-reads/TR-primer-trimmed_R1.fastq.gz",
        R2 = f"{outdir}/samples/{{sample}}/mhaps-reads/TR-primer-trimmed_R2.fastq.gz",
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

rule fastq_merge_to_fq:
    threads: config.get('vsearch_threads', 4)
    input:
        R1 = f"{outdir}/samples/{{sample}}/mhaps-reads/TR-primer-trimmed_R1.fastq.gz",
        R2 = f"{outdir}/samples/{{sample}}/mhaps-reads/TR-primer-trimmed_R2.fastq.gz",
    output:
        merged = f"{outdir}/samples/{{sample}}/fastq_merge/final_merged.fastq.gz",
    log:
        f"{outdir}/samples/{{sample}}/logs/fastq_merge.log",
    params:
        fastq_merge_parameters = config.get('fastq_merge_parameters', ''),
    shell:
        "vsearch --fastq_mergepairs {input.R1} --reverse {input.R2} {params.fastq_merge_parameters} --fastqout {output.merged} > {log} 2>&1"

rule ngmerge_to_fq:
    threads: config.get('vsearch_threads', 4)
    input:
        R1 = f"{outdir}/samples/{{sample}}/mhaps-reads/TR-primer-trimmed_R1.fastq.gz",
        R2 = f"{outdir}/samples/{{sample}}/mhaps-reads/TR-primer-trimmed_R2.fastq.gz",
    output:
        merged = temp(f"{outdir}/samples/{{sample}}/ngmerge/merged.fastq.gz"),
    log:
        f"{outdir}/samples/{{sample}}/logs/ngmerge.log",
    params:
        ngmerge_parameters = config.get('ngmerge_parameters', ''),
    shell:
        "ngmerge -1 {input.R1} -2 {input.R2} -o {output.merged} {params.ngmerge_parameters} -v -l {log} -n {threads}"

rule filter_ngmerged_to_fq:
    threads: config.get('vsearch_threads', 4)
    input:
        merged = f"{outdir}/samples/{{sample}}/ngmerge/merged.fastq.gz",
    output:
        filtered = f"{outdir}/samples/{{sample}}/ngmerge/final_merged.fastq.gz",
    log:
        f"{outdir}/samples/{{sample}}/logs/ngmerge_filter.log",
    shell:
        "vsearch --fastq_filter {input.merged} --fastq_maxee 1 --fasta_width 0 --threads {threads} --fastqout {output.filtered} > {log} 2>&1"

# Denoising with vsearch
# reference: https://learnmetabarcoding.github.io/LearnMetabarcoding/filtering/freq_filtering_and_denoising.html
rule indv_bbfastq_to_fasta:
    input:
        merged = f"{outdir}/samples/{{sample}}/{bbmap_or_bbmerge_fastq_merge}/final_merged.fastq.gz",
    output:
        fasta = f"{outdir}/samples/{{sample}}/{bbmap_or_bbmerge_fastq_merge}/final_merged.fasta",
    log:
        f"{outdir}/samples/{{sample}}/logs/fastx_uniques.log",
    shell:
        "vsearch --fastx_uniques {input.merged} --fastaout {output.fasta} --fasta_width 0 --sizeout --relabel {wildcards.sample}: > {log} 2>&1"

rule dereplicated_fa_to_uc:
    threads: config.get('vsearch_threads', 4)
    input:
        fasta = f"{outdir}/samples/{{sample}}/{bbmap_or_bbmerge_fastq_merge}/final_merged.fasta",
        db = insertseq_file,
    output:
        uc = f"{outdir}/samples/{{sample}}/{bbmap_or_bbmerge_fastq_merge}/dereplicated_counted.uc",
    log:
        f"{outdir}/samples/{{sample}}/logs/vsearch_usearch_global.log",
    shell:
        "vsearch --usearch_global {input.fasta} --db {input.db} --id 0.8 --top_hits_only --uc {output.uc} --threads {threads} > {log} 2>&1"

checkpoint split_dereplicated_fa_uc_to_marker_fa:
    input:
        uc = f"{outdir}/samples/{{sample}}/{bbmap_or_bbmerge_fastq_merge}/dereplicated_counted.uc",
        fasta = f"{outdir}/samples/{{sample}}/{bbmap_or_bbmerge_fastq_merge}/final_merged.fasta",
        insertseq = insertseq_file,
    output:
        output_dir = directory(f"{outdir}/samples/{{sample}}/{bbmap_or_bbmerge_fastq_merge}/denoise"),
    shell:
        """
        mkdir -p {output.output_dir} &&
        python {microhaps_basedir}/scripts/split_dereplicated.py \
         -u {input.uc} \
         -f {input.fasta} \
         -i {input.insertseq} \
         -o {output.output_dir}
        """

# if the research questions are tolerant to the loss of rare real sequences and/or
# it is important that as high a proportion of ASVs be valid as possible,
# then stricter settings (higher --minsize, lower --unoise_alpha) would be sensible. 
# alpha - determines the threshold level of dissimilarity between frequent and infrequent reads
# for exclusion of infrequent reads.

# https://drive5.com/usearch/manual/cmd_unoise3.html <- need to denoise per sample
# https://drive5.com/usearch/manual/faq_combine_samples.html <- per this, if same protocols, etc. pool and no need to denoise per sample
# https://drive5.com/usearch/manual/pipe_otutab.html <- from denoised fasta to otutab <- mapping all available sequences to zotu (denoised)

rule cluster_fasta:
    threads: config.get('bbmap_merge_threads', 4)
    input:
        fasta = f"{outdir}/samples/{{sample}}/{bbmap_or_bbmerge_fastq_merge}/denoise/{{marker}}-split.fa",
    output:
        fasta = f"{outdir}/samples/{{sample}}/{bbmap_or_bbmerge_fastq_merge}/denoise/{{marker}}-split-denoised.fa",
        uc = f"{outdir}/samples/{{sample}}/{bbmap_or_bbmerge_fastq_merge}/denoise/{{marker}}-split-denoised.uc",
    group: "cluster_chimera_removal"
    log:
       f"{outdir}/samples/{{sample}}/logs/vsearch_cluster-{{marker}}.log",
    params:
        minsize = config.get("unoise_minsize", 4),  # minimum size of clusters to be kept
        unoise_alpha = config.get("unoise_alpha", 2),  # alpha parameter for unoise algorithm
    shell:
        "vsearch --cluster_unoise {input.fasta} --minsize {params.minsize} --fasta_width 0 --unoise_alpha {params.unoise_alpha} --centroids {output.fasta} --uc {output.uc} > {log} 2>&1"

# chimeric removal
rule chimera_filtering:
    threads: config.get('vsearch_threads', 4)
    input:
        fasta = f"{outdir}/samples/{{sample}}/{bbmap_or_bbmerge_fastq_merge}/denoise/{{marker}}-split-denoised.fa",
    output:
        kept = f"{outdir}/samples/{{sample}}/{bbmap_or_bbmerge_fastq_merge}/denoise/{{marker}}-split-denoised-no_chimera.fa",
        discarded = f"{outdir}/samples/{{sample}}/{bbmap_or_bbmerge_fastq_merge}/denoise/{{marker}}-split-denoised-chimera.fa",
    group: "cluster_chimera_removal"
    log:
        f"{outdir}/samples/{{sample}}/logs/vsearch_chimera-{{marker}}.log",
    shell:
        "vsearch --uchime3_denovo {input.fasta} --nonchimeras {output.kept} --chimeras {output.discarded} --fasta_width 0 --threads {threads} > {log} 2>&1"

def get_split_files(wildcards):
    checkpoint_output = checkpoints.split_dereplicated_fa_uc_to_marker_fa.get(**wildcards).output[0]
    return expand(f"{outdir}/samples/{wildcards.sample}/{bbmap_or_bbmerge_fastq_merge}/denoise/{{marker}}-split-denoised-no_chimera.fa", marker=markers)

rule indv_merge_denoise_fasta:
    input:
        fasta = get_split_files,
    output:
        merged = temp(f"{outdir}/samples/{{sample}}/{bbmap_or_bbmerge_fastq_merge}/merged-denoised-no_chimera.fa"),
    shell:
        "cat {input.fasta} > {output.merged} "

rule merge_denoised_fasta:
    input:
        fasta = expand(f"{outdir}/samples/{{sample}}/{bbmap_or_bbmerge_fastq_merge}/merged-denoised-no_chimera.fa", sample = IDs),
    output:
        merged = temp(f"{outdir}/malamp/{bbmap_or_bbmerge_fastq_merge}/dereplicated_counted_no_chimeras.fa"),
        unique = f"{outdir}/malamp/{bbmap_or_bbmerge_fastq_merge}/dereplicated_counted_no_chimeras-unique.fa",
    shell:
        "cat {input.fasta} > {output.merged} && "
        "vsearch --fastx_uniques {output.merged} --fastaout {output.unique} --fasta_width 0 --sizein --sizeout --relabel_sha1 "

rule indv_otutab:
    threads: config.get('vsearch_threads', 4)
    input:
        fasta = f"{outdir}/malamp/{bbmap_or_bbmerge_fastq_merge}/dereplicated_counted_no_chimeras-unique.fa",
        sample_fa = f"{outdir}/samples/{{sample}}/{bbmap_or_bbmerge_fastq_merge}/final_merged.fasta",
    output:
        otutab = f"{outdir}/samples/{{sample}}/{bbmap_or_bbmerge_fastq_merge}/seqtab.tsv",
    shell:
        "vsearch --search_exact {input.sample_fa} --db {input.fasta} --otutabout {output.otutab} --threads {threads}"

rule merge_otutab:
    input:
        otutabs = expand(f"{outdir}/samples/{{sample}}/{bbmap_or_bbmerge_fastq_merge}/seqtab.tsv", sample=IDs),
        filtered_fasta = f"{outdir}/malamp/{bbmap_or_bbmerge_fastq_merge}/dereplicated_counted_no_chimeras-unique.fa",
    output:
        otutab = f"{outdir}/malamp/{bbmap_or_bbmerge_fastq_merge}/seqtab.tsv",
    shell:
        "python {microhaps_basedir}/scripts/merge_otutab.py"
        " --denoised_fasta {input.filtered_fasta}"
        " --out_seqtab {output.otutab}"
        " --indv_seqtabs {input.otutabs}"

