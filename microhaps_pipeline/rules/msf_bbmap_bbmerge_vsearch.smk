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

rule trim_R_optical_dedup:
    input:
        R1 = f"{outdir}/samples/{{sample}}/mhaps-reads/primer-trimmed_R1.fastq.gz",
        R2 = f"{outdir}/samples/{{sample}}/mhaps-reads/primer-trimmed_R2.fastq.gz",
    output:
        Trimmed_R1 = temp(f"{outdir}/samples/{{sample}}/mhaps-reads/R-primer-trimmed_R1.fastq.gz"),
        Trimmed_R2 = temp(f"{outdir}/samples/{{sample}}/mhaps-reads/R-primer-trimmed_R2.fastq.gz"),
        Op_dedup_R1 = temp(f"{outdir}/samples/{{sample}}/mhaps-reads/TR-primer-trimmed_R1.fastq.gz"),
        Op_dedup_R2 = temp(f"{outdir}/samples/{{sample}}/mhaps-reads/TR-primer-trimmed_R2.fastq.gz"),
    params:
        instrument_specific = "spany adjacent" if config.get("instrument", "generic") == "nextseq" else "",
    shell:
        """
        vsearch --fastx_filter {input.R1} --reverse {input.R2} --fastqout >(gzip -c  > {output.Trimmed_R1}) \
        --fastqout_rev >(gzip -c > {output.Trimmed_R2}) --fastq_stripright 10 --fastq_maxee 1 &&
        clumpify.sh in1={output.Trimmed_R1} in2={output.Trimmed_R2} out1={output.Op_dedup_R1} out2={output.Op_dedup_R2} dedupe optical {params.instrument_specific}
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
        filtered = f"{outdir}/samples/{{sample}}/ngmerge/final_merged.fastq.gz",
    log:
        ngmerge = f"{outdir}/samples/{{sample}}/logs/ngmerge.log",
        ngmerge_filter = f"{outdir}/samples/{{sample}}/logs/ngmerge_filter.log",
    params:
        ngmerge_parameters = config.get('ngmerge_parameters', ''),
    shell:
        """
        NGmerge -1 {input.R1} -2 {input.R2} -o {output.merged} {params.ngmerge_parameters} -v -l {log.ngmerge} -n {threads} && 
        vsearch --fastq_filter {output.merged} --fastq_maxee 1 --fasta_width 0 --threads {threads} --fastqout {output.filtered} > {log.ngmerge_filter} 2>&1
        """

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
    params:
        min_id_amplicon = config.get('min_id_amplicon', 0.8), 
    shell:
        "vsearch --usearch_global {input.fasta} --db {input.db} --id {params.min_id_amplicon} --top_hits_only --uc {output.uc} --threads {threads} > {log} 2>&1"

rule split_dereplicated_fa_uc_to_marker_fa:
    threads: 4
    input:
        uc = expand(f"{outdir}/samples/{{sample}}/{bbmap_or_bbmerge_fastq_merge}/dereplicated_counted.uc", sample=IDs),
        fasta = expand(f"{outdir}/samples/{{sample}}/{bbmap_or_bbmerge_fastq_merge}/final_merged.fasta", sample=IDs),
        insertseq = insertseq_file,
    output:
        sample_output_dir = directory(expand(f"{outdir}/samples/{{sample}}/{bbmap_or_bbmerge_fastq_merge}/denoise", sample=IDs)),
        output_dir = directory(f"{outdir}/malamp/{bbmap_or_bbmerge_fastq_merge}/denoise"),
        indv_merged_denoised = f"{outdir}/malamp/{bbmap_or_bbmerge_fastq_merge}/dereplicated_counted_no_chimeras-unique.fa",
        temp_indv_merged_denoised = temp(f"{outdir}/malamp/{bbmap_or_bbmerge_fastq_merge}/dereplicated_counted_no_chimeras-unique.fa.temp"),
    params:
        minsize = config.get("unoise_minsize", 4),  # minimum size of clusters to be kept
        unoise_alpha = config.get("unoise_alpha", 2),  # alpha parameter for unoise algorithm
        log_dir = f"{outdir}/logs",
    run:
        from concurrent.futures import ThreadPoolExecutor
        from glob import glob
        import re

        def demux_per_sample(uc, fasta, sample_outdir):
            print(uc, fasta, sample_outdir)
            shell(f"mkdir -p {sample_outdir}")
            # Split the amplicon sequences from each samples
            shell(f"python {microhaps_basedir}/scripts/split_dereplicated.py \
            -u {uc} \
            -f {fasta} \
            -i {input.insertseq} \
            -o {sample_outdir}")


        with ThreadPoolExecutor(max_workers=threads) as executor:
            demux_result = executor.map(demux_per_sample, input.uc, input.fasta, output.sample_output_dir)


        # Join the sequences from multiple samples into a single fasta for each marker
        shell(f"mkdir -p {output.output_dir}")

        fastas = sum([glob(f"{sample_output_dir}/*-split.fa") for sample_output_dir in output.sample_output_dir], [])
        markers = [re.search(r".*/(.*?)-split.fa", f).group(1) for f in fastas]
        u_markers = list(set(markers))

        for marker in u_markers:
            index = [i for i, m in enumerate(markers) if m == marker]
            fasta_marker = [fastas[i] for i in index]
            sublist_fastas = " ".join(fasta_marker)
            shell(f"cat {sublist_fastas} > {output.output_dir}/{marker}-split-non_unique.fa")
            shell(f"vsearch --fastx_uniques {output.output_dir}/{marker}-split-non_unique.fa --fastaout {output.output_dir}/{marker}-split.fa --fasta_width 0 --sizein --sizeout --relabel_sha1 ")

        def run_vsearch_denoise_uchime(fasta, marker):
            cmd = f"""vsearch --cluster_unoise {fasta} --centroids {fasta.replace('-split.fa', '-split-denoised.fa')} \
                --fasta_width 0 \
                --minsize {params.minsize} --unoise_alpha {params.unoise_alpha} \
                """
                # > {params.log_dir}/vsearch_denoise-{marker}.log
            shell(cmd)
            cmd2 = f"""vsearch --uchime3_denovo {fasta.replace('-split.fa', '-split-denoised.fa')} \
                --nonchimeras {fasta.replace('-split.fa', '-split-denoised-no_chimera.fa')} \
                --chimeras {fasta.replace('-split.fa', '-split-denoised-chimera.fa')} \
                --fasta_width 0 --threads {threads} \
                """
                #> {params.log_dir}/vsearch_chimera-{marker}.log
            shell(cmd2)

        final_fasta = [f"{output.output_dir}/{marker}-split.fa" for marker in u_markers]
        with ThreadPoolExecutor(max_workers=threads) as executor:
            print("Denoising & Chimera filtering with vsearch...")
            denoise_results = executor.map(run_vsearch_denoise_uchime, final_fasta, u_markers)

        shell(f"cat {output.output_dir}/*-split-denoised-no_chimera.fa > {output.indv_merged_denoised}.temp")
        shell(f"vsearch --fastx_uniques {output.indv_merged_denoised}.temp --fastaout {output.indv_merged_denoised} --fasta_width 0 --sizein --sizeout --relabel_sha1 ")
        # swarm test:
        # swarm -d 1 -z -f --boundary 3 -i swarm-26982-in.tab -j swarm-26982-net.net -s swarm-26982-stats.tsv -w swarm-26982-centroid.fa malamp/bbmerge/denoise/PvP01_02_v2\:327547-327756\:26982-split.fa
        # normalise to 1000 per sample and unoise

# if the research questions are tolerant to the loss of rare real sequences and/or
# it is important that as high a proportion of ASVs be valid as possible,
# then stricter settings (higher --minsize, lower --unoise_alpha) would be sensible. 
# alpha - determines the threshold level of dissimilarity between frequent and infrequent reads
# for exclusion of infrequent reads.

# https://drive5.com/usearch/manual/cmd_unoise3.html <- need to denoise per sample
# https://drive5.com/usearch/manual/faq_combine_samples.html <- per this, if same protocols, etc. pool and no need to denoise per sample
# https://drive5.com/usearch/manual/pipe_otutab.html <- from denoised fasta to otutab <- mapping all available sequences to zotu (denoised)

# rule cluster_fasta:
#     threads: config.get('bbmap_merge_threads', 4)
#     input:
#         fasta = f"{outdir}/samples/{{sample}}/{bbmap_or_bbmerge_fastq_merge}/denoise/{{marker}}-split.fa",
#     output:
#         fasta = f"{outdir}/samples/{{sample}}/{bbmap_or_bbmerge_fastq_merge}/denoise/{{marker}}-split-denoised.fa",
#         uc = f"{outdir}/samples/{{sample}}/{bbmap_or_bbmerge_fastq_merge}/denoise/{{marker}}-split-denoised.uc",
#     resources:
#         runtime = "1h",
#     group: "cluster_chimera_removal"
#     log:
#        f"{outdir}/samples/{{sample}}/logs/vsearch_cluster-{{marker}}.log",
#     params:
#         minsize = config.get("unoise_minsize", 4),  # minimum size of clusters to be kept
#         unoise_alpha = config.get("unoise_alpha", 2),  # alpha parameter for unoise algorithm
#     shell:
#         "vsearch --cluster_unoise {input.fasta} --minsize {params.minsize} --fasta_width 0 --unoise_alpha {params.unoise_alpha} --centroids {output.fasta} --uc {output.uc} > {log} 2>&1"

# chimeric removal
# rule chimera_filtering:
#     threads: config.get('vsearch_threads', 4)
#     input:
#         fasta = f"{outdir}/samples/{{sample}}/{bbmap_or_bbmerge_fastq_merge}/denoise/{{marker}}-split-denoised.fa",
#     output:
#         kept = f"{outdir}/samples/{{sample}}/{bbmap_or_bbmerge_fastq_merge}/denoise/{{marker}}-split-denoised-no_chimera.fa",
#         discarded = f"{outdir}/samples/{{sample}}/{bbmap_or_bbmerge_fastq_merge}/denoise/{{marker}}-split-denoised-chimera.fa",
#     group: "cluster_chimera_removal"
#     resources:
#         runtime = "1h",
#     log:
#         f"{outdir}/samples/{{sample}}/logs/vsearch_chimera-{{marker}}.log",
#     shell:
#         "vsearch --uchime3_denovo {input.fasta} --nonchimeras {output.kept} --chimeras {output.discarded} --fasta_width 0 --threads {threads} > {log} 2>&1"

# def get_split_files(wildcards):
#     checkpoint_output = checkpoints.split_dereplicated_fa_uc_to_marker_fa.get(**wildcards).output[0]
#     return expand(f"{outdir}/samples/{wildcards.sample}/{bbmap_or_bbmerge_fastq_merge}/denoise/{{marker}}-split-denoised-no_chimera.fa", marker=markers)

# rule indv_merge_denoise_fasta:
#     input:
#         fasta = get_split_files,
#     output:
#         merged = temp(f"{outdir}/samples/{{sample}}/{bbmap_or_bbmerge_fastq_merge}/merged-denoised-no_chimera.fa"),
#     shell:
#         "cat {input.fasta} > {output.merged} "

# rule merge_denoised_fasta:
#     input:
#         fasta = expand(f"{outdir}/samples/{{sample}}/{bbmap_or_bbmerge_fastq_merge}/merged-denoised-no_chimera.fa", sample = IDs),
#     output:
#         merged = temp(f"{outdir}/malamp/{bbmap_or_bbmerge_fastq_merge}/dereplicated_counted_no_chimeras.fa"),
#         unique = f"{outdir}/malamp/{bbmap_or_bbmerge_fastq_merge}/dereplicated_counted_no_chimeras-unique.fa",
#     shell:
#         "cat {input.fasta} > {output.merged} && "
#         "vsearch --fastx_uniques {output.merged} --fastaout {output.unique} --fasta_width 0 --sizein --sizeout --relabel_sha1 "

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

