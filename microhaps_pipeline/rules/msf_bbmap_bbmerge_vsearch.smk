bbmap_or_bbmerge_fastq_merge = config.get('merge_map') # should be either "bbmap" or "bbmerge" or "fastq_merge"
markers = [line.replace(">", "") for line in open(insertseq_file, 'r').read().split() if line.startswith('>')]

ruleorder:
    indv_bbfastq_to_fasta > bunzip2


# bbmap -> bbmerge: tested, but no improvement over bbmerge only
# reference: https://www.seqanswers.com/forum/bioinformatics/bioinformatics-aa/40046-non-overlapping-read-assembly?t=45477

rule bbmerge_to_fq:
    threads: config.get('bbmmerge_threads', 4)
    input:
        R1 = f"{outdir}/samples/{{sample}}/mhaps-reads/trimmed-filtered_R1.fastq.gz",
        R2 = f"{outdir}/samples/{{sample}}/mhaps-reads/trimmed-filtered_R2.fastq.gz",
    output:
        merged = f"{outdir}/samples/{{sample}}/bbmerge/final_merged.corrected.fastq.gz", #final_merged.fastq.gz",
    log:
        f"{outdir}/samples/{{sample}}/logs/bbmerge.log",
    params:
        outdir = outdir,
        bbmerge_parameters = config.get('bbmerge_parameters', ''),
    shell:
        "bbmerge.sh in1={input.R1} in2={input.R2} outm={output.merged} t={threads} path={params.outdir} {params.bbmerge_parameters} > {log} 2>&1 "
        # po = paired only, "rbm" and "don" (renamebymapping and deleteoldname)
        # "bfc -t {threads} {output.merged} | gzip -1 > {output.corrected}"

rule fastq_merge_to_fq:
    threads: config.get('vsearch_threads', 2)
    input:
        R1 = f"{outdir}/samples/{{sample}}/mhaps-reads/trimmed-filtered_R1.fastq.gz",
        R2 = f"{outdir}/samples/{{sample}}/mhaps-reads/trimmed-filtered_R2.fastq.gz",
    output:
        merged = f"{outdir}/samples/{{sample}}/fastq_merge/final_merged.corrected.fastq.gz", # final_merged.fastq.gz",
    log:
        f"{outdir}/samples/{{sample}}/logs/fastq_merge.log",
    params:
        fastq_merge_parameters = config.get('fastq_merge_parameters', ''),
    shell:
        "vsearch --fastq_mergepairs {input.R1} --reverse {input.R2} {params.fastq_merge_parameters} --fastqout {output.merged} > {log} 2>&1 "


# rule bfc_count:
#     threads: 4
#     input:
#         expand(f"{outdir}/samples/{{sample}}/{bbmap_or_bbmerge_fastq_merge}/final_merged.fastq.gz", sample=IDs)
#     output:
#         all_merged = temp(f"{outdir}/samples/all_merged.fastq.gz"),
#         hashfile = temp(f"{outdir}/hashfile.bin"),
#     log:
#         f"{outdir}/logs/bfc_count.log",
#     shell:
#         "cat {input} > {output.all_merged} && "
#         "bfc -t {threads} -E -d {output.hashfile} {output.all_merged} 2> {log} "

# rule bfc_correct:
#     input:
#         hashfile = f"{outdir}/hashfile.bin",
#         merged = f"{outdir}/samples/{{sample}}/{bbmap_or_bbmerge_fastq_merge}/final_merged.fastq.gz",
#     output:
#         corrected = f"{outdir}/samples/{{sample}}/{bbmap_or_bbmerge_fastq_merge}/final_merged.corrected.fastq.gz",
#     log:
#         f"{outdir}/samples/{{sample}}/logs/bfc_correct.log"
#     shell:
#         "bfc -r {input.hashfile} {input.merged} 2> {log} | gzip -1 > {output.corrected}"

 #corrected = f"{outdir}/samples/{{sample}}/fastq_merge/final_merged.corrected.fastq.gz",


# Denoising with vsearch
# reference: https://learnmetabarcoding.github.io/LearnMetabarcoding/filtering/freq_filtering_and_denoising.html
rule indv_bbfastq_to_fasta:
    input:
        merged = f"{outdir}/samples/{{sample}}/{bbmap_or_bbmerge_fastq_merge}/final_merged.corrected.fastq.gz",
    output:
        fasta = f"{outdir}/samples/{{sample}}/{bbmap_or_bbmerge_fastq_merge}/final_merged.fasta",
    log:
        f"{outdir}/samples/{{sample}}/logs/fastx_uniques.log",
    shell:
        "vsearch --fastx_uniques {input.merged} --fastaout {output.fasta} --fasta_width 0 --sizeout --sample {wildcards.sample} > {log} 2>&1"

rule dereplicated_fa_to_uc:
    threads: config.get('vsearch_threads', 2)
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
        """
        vsearch --usearch_global {input.fasta} --db {input.db} --id {params.min_id_amplicon} \
        --maxaccepts 0 --maxrejects 0 --maxhits 1 \
        --uc {output.uc} --threads {threads} > {log} 2>&1
        """

rule split_dereplicated_fa_uc_to_marker_fa:
    threads: 2
    input:
        uc = f"{outdir}/samples/{{sample}}/{bbmap_or_bbmerge_fastq_merge}/dereplicated_counted.uc",
        fasta = f"{outdir}/samples/{{sample}}/{bbmap_or_bbmerge_fastq_merge}/final_merged.fasta",
        insertseq = insertseq_file,
    output:
        sample_output_dir = directory(f"{outdir}/samples/{{sample}}/{bbmap_or_bbmerge_fastq_merge}/denoise"),
        indv_merged_denoised = f"{outdir}/samples/{{sample}}/{bbmap_or_bbmerge_fastq_merge}/dereplicated_counted_no_chimeras-unique.fa",
        temp_indv_merged_denoised = temp(f"{outdir}/samples/{{sample}}/{bbmap_or_bbmerge_fastq_merge}/dereplicated_counted_no_chimeras-unique.fa.temp",),
    params:
        minsize = config.get("unoise_minsize", 4),  # minimum size of clusters to be kept
        minsize_ratio = config.get('unoise_minsize_ratio', 1), # final minsize if minsize_ratio < 1 = max(minsize_ratio * total size, minsize)
        unoise_alpha = config.get("unoise_alpha", 2),  # alpha parameter for unoise algorithm
        log_dir = f"{outdir}/samples/{{sample}}/logs",
    run:
        from concurrent.futures import ThreadPoolExecutor
        from glob import glob
        import re
        import math

        # This split the amplicon sequences from each samples into different markers (based on best global pairwise alignment)
        def demux_per_sample(uc, fasta, sample_outdir):
            print(uc, fasta, sample_outdir)
            shell(f"mkdir -p {sample_outdir}")
            shell(f"python {microhaps_basedir}/scripts/split_dereplicated.py \
            -u {uc} \
            -f {fasta} \
            -i {input.insertseq} \
            -o {sample_outdir}")

        demux_per_sample(input.uc, input.fasta, output.sample_output_dir)

        fastas = glob(f"{output.sample_output_dir}/*-split.fa")
        markers = [re.search(r".*/(.*?)-split.fa", f).group(1) for f in fastas]
        u_markers = list(set(markers))
        print(fastas)
        marker_fastas = []
        final_fasta = [f"{output.sample_output_dir}/{marker}-split.fa" for marker in u_markers]
        for marker in u_markers:
            index = [i for i, m in enumerate(markers) if m == marker]
            fasta_marker = [fastas[i] for i in index]
            marker_fastas.append((marker, fasta_marker))

        ### Now we have multiple Fasta consisting of unique merged reads with total read counts (for a sample)
        ### Rescue: only if <uchime_unoise>
        ##  already available: (1) unique sesequence size  (after uchime)
        ##  calculate: (2) total size of sequence from (1)
        ##  calculate: cluster_size2 = cluster_ratio * total size (3) 
        ##  put sequence into centroids (Zotus) if size (1) >= min of (2), (3)
        ##  ----------
        def run_vsearch_denoise_uchime(fasta, marker,
            order='unoise_uchime', rescue_cluster=False,
            rescue_params={"cluster_size": 10, "cluster_ratio": 0.01}):
            
            def get_size_from_fasta(fasta):
                with open(fasta, 'r') as f:
                    for line in f:
                        if line.startswith(">"):
                            match = re.search(r"size=(\d+)", line)
                            if match:
                                yield int(match.group(1))
            
            if params.minsize_ratio < 1:
                total_read_in_fa = sum([read_size for read_size in get_size_from_fasta(fasta)])
                unoise_minsize = max( math.floor(total_read_in_fa * params.minsize_ratio), params.minsize)
                print(f"Dynamic unoise minsize for {fasta}: {unoise_minsize}")
            else:
                unoise_minsize = params.minsize

            if order == 'unoise_uchime':
                # start -> denoise -> chimera_removal
                start_fasta = fasta
                intermediate_fasta = start_fasta.replace('-split.fa', '-split-denoised.fa')
                final_fasta = start_fasta.replace('-split.fa', '-split-denoised-no_chimera.fa')
                final_fasta_chimera = start_fasta.replace('-split.fa', '-split-denoised-chimera.fa')

                denoise = f"""vsearch --cluster_unoise {start_fasta} --centroids {intermediate_fasta} \
                    --fasta_width 0 \
                    --minsize {unoise_minsize} --unoise_alpha {params.unoise_alpha} \
                    --threads 1 > {params.log_dir}/vsearch_denoise-{marker}.log 2>&1
                    """
                shell(denoise)
                chmimera_removal = f"""vsearch --uchime3_denovo {intermediate_fasta} \
                    --qmask none \
                    --nonchimeras {final_fasta} \
                    --chimeras {final_fasta_chimera} \
                    --fasta_width 0 --threads 1 > {params.log_dir}/vsearch_chimera-{marker}.log 2>&1
                    """
                shell(chmimera_removal)

            if order == 'uchime_unoise':
                # start -> chimera_removal -> denoise
                start_fasta = fasta
                final_fasta_chimera = start_fasta.replace('-split.fa', '-split-chimera.fa')
                intermediate_fasta = start_fasta.replace('-split.fa', '-split-no_chimera.fa')
                final_fasta = start_fasta.replace('-split.fa', '-split-denoised-no_chimera.fa')
                
                chmimera_removal = f"""vsearch --uchime3_denovo {start_fasta} \
                    --qmask none \
                    --nonchimeras {intermediate_fasta} \
                    --chimeras {final_fasta_chimera} \
                    --fasta_width 0 --threads 1 > {params.log_dir}/vsearch_chimera-{marker}.log 2>&1
                    """
                shell(chmimera_removal)
                denoise = f"""vsearch --cluster_unoise {intermediate_fasta} --centroids {final_fasta} \
                    --fasta_width 0 \
                    --minsize {unoise_minsize} --unoise_alpha {params.unoise_alpha} \
                    --threads 1 > {params.log_dir}/vsearch_denoise-{marker}.log 2>&1
                    """
                shell(denoise)

        with ThreadPoolExecutor(max_workers=threads) as executor:
            print("Denoising & Chimera filtering with vsearch...")
            executor.map(run_vsearch_denoise_uchime, final_fasta, u_markers)

        all_outputs = glob(f"{output.sample_output_dir}/*-split-denoised-no_chimera.fa")

        if len(all_outputs) > 1:
            shell(f"cat {' '.join(all_outputs)} > {output.temp_indv_merged_denoised}")
            shell(f"vsearch --fastx_uniques {output.temp_indv_merged_denoised} --fastaout {output.indv_merged_denoised} --sample {wildcards.sample} --fasta_width 0 --sizein --sizeout --relabel_sha1 ")
        else:
            if len(fastas) == 0:
                print(f"{wildcards.sample} did not have any read-pairs")
            else:
                print(f"{wildcards.sample} did not have any read-pairs that pass the denoising & chimera removal")
            shell(f"touch {output.temp_indv_merged_denoised} {output.indv_merged_denoised}")

########################################## NOTES ##########################################
# if the research questions are tolerant to the loss of rare real sequences and/or
# it is important that as high a proportion of ASVs be valid as possible,
# then stricter settings (higher --minsize, lower --unoise_alpha) would be sensible. 
# alpha - determines the threshold level of dissimilarity between frequent and infrequent reads
# for exclusion of infrequent reads.

# https://drive5.com/usearch/manual/cmd_unoise3.html <- need to denoise per sample
# https://drive5.com/usearch/manual/faq_combine_samples.html <- per this, if same protocols, etc. pool and no need to denoise per sample
# https://drive5.com/usearch/manual/pipe_otutab.html <- from denoised fasta to otutab <- mapping all available sequences to zotu (denoised)
########################################## NOTES ##########################################


rule merge_denoised_fasta:
    input:
        fasta = expand(f"{outdir}/samples/{{sample}}/{bbmap_or_bbmerge_fastq_merge}/dereplicated_counted_no_chimeras-unique.fa", sample = IDs),
    output:
        merged = temp(f"{outdir}/malamp/{bbmap_or_bbmerge_fastq_merge}/dereplicated_counted_no_chimeras.fa"),
        unique = temp(f"{outdir}/malamp/{bbmap_or_bbmerge_fastq_merge}/dereplicated_counted_no_chimeras-unique.unsorted.fa"),
        sorted_unique = f"{outdir}/malamp/{bbmap_or_bbmerge_fastq_merge}/dereplicated_counted_no_chimeras-unique.fa",
    shell:
        "cat {input.fasta} > {output.merged} && "
        "vsearch --fastx_uniques {output.merged} --fastaout {output.unique} --fasta_width 0 --sizein --sizeout --relabel_sha1 &&"
        "vsearch --sortbysize {output.unique} --output {output.sorted_unique}"

rule indv_otutab:
    threads: config.get('vsearch_threads', 2)
    input:
        fasta = f"{outdir}/malamp/{bbmap_or_bbmerge_fastq_merge}/dereplicated_counted_no_chimeras-unique.fa",
        sample_fa = f"{outdir}/samples/{{sample}}/{bbmap_or_bbmerge_fastq_merge}/final_merged.fasta",
    output:
        otutab = f"{outdir}/samples/{{sample}}/{bbmap_or_bbmerge_fastq_merge}/seqtab.tsv",
    params:
        sequence_to_denoise_criteria = config.get('sequence_to_denoise_criteria', "--id 0.98 --iddef 1 --target_cov 0.98 --query_cov 0.98")
    shell:
        """
        vsearch --usearch_global {input.sample_fa} --db {input.fasta} {params.sequence_to_denoise_criteria} \
        --otutabout {output.otutab} \
        --maxaccepts 0 --maxrejects 0 --maxhits 1 \
        --sizein --sizeout --fasta_width 0 --qmask none --dbmask none --threads {threads}
        """

rule merge_otutab:
    input:
        otutabs = expand(f"{outdir}/samples/{{sample}}/{bbmap_or_bbmerge_fastq_merge}/seqtab.tsv", sample=IDs),
        filtered_fasta = f"{outdir}/malamp/{bbmap_or_bbmerge_fastq_merge}/dereplicated_counted_no_chimeras-unique.fa",
    output:
        otutab = f"{outdir}/malamp/{bbmap_or_bbmerge_fastq_merge}/seqtab.tsv",
    shell:
        """
        python {microhaps_basedir}/scripts/merge_otutab.py \
         --denoised_fasta {input.filtered_fasta} \
         --out_seqtab {output.otutab} \
         --indv_seqtabs {input.otutabs}
        """