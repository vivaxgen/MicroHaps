primers_trimmed = config.get("primers_trimmed", False)

wildcard_constraints:
    marker=f"({'|'.join(Markers)})"

def get_min_max_len(marker):
    tolerant_threshold = config.get("insert_size_tolerant_threshold", 0.25)
    import pandas as pd
    all_markers = pd.read_table(targetregion_file, header=None, names=["Chr", "Start", "End", "Amplicon_name"])
    all_markers["length"] = all_markers["End"] - all_markers["Start"] + 1
    rel_marker = all_markers.query("Amplicon_name == @marker")["length"].values[0]
    min_len = int(rel_marker * (1 - tolerant_threshold))
    max_len = int(rel_marker * (1 + tolerant_threshold))
    return min_len, max_len

rule merge_sample_marker:
    group: "merge_filter_unique_{sample}"
    input:
        R1 = f"{outdir}/samples/{{sample}}/markers-reads/{{marker}}/trimmed-filtered_R1.fastq.gz",
        R2 = f"{outdir}/samples/{{sample}}/markers-reads/{{marker}}/trimmed-filtered_R2.fastq.gz",
    output:
        merged = f"{outdir}/samples/{{sample}}/markers-reads/{{marker}}/merged.corrected.fastq.gz",
    params:
        fastp_merge_params = config.get("fastp_merge_params", "")
    log:
        fastp_html = f"{outdir}/samples/{{sample}}/markers-reads/{{marker}}/fastp.html",
        fastp_json = f"{outdir}/samples/{{sample}}/markers-reads/{{marker}}/fastp.json",
    shell:
        """
        fastp -A -Q -L -m --in1 {input.R1} --in2 {input.R2} --merged_out {output.merged} -j {log.fastp_json} -h {log.fastp_html} {params.fastp_merge_params}
        """

rule vsearch_filter_insert_length:
    group: "merge_filter_unique_{sample}"
    input:
        merged = f"{outdir}/samples/{{sample}}/markers-reads/{{marker}}/merged.corrected.fastq.gz",
    output:
        merged_fasta = f"{outdir}/samples/{{sample}}/markers-reads/{{marker}}/merged.fasta",
    params:
        min_max = lambda w: get_min_max_len(w.marker),
    log:
        vsearch_filter_log = f"{outdir}/samples/{{sample}}/markers-reads/{{marker}}/vsearch_filter.log",
    shell:
        """
        vsearch --fastx_filter {input.merged} --fastaout {output.merged_fasta} --fastq_minlen {params.min_max[0]} --fastq_maxlen {params.min_max[1]} > {log.vsearch_filter_log} 2>&1
        """

rule vsearch_filter_unique:
    group: "merge_filter_unique_{sample}"
    input:
        merged_fasta = f"{outdir}/samples/{{sample}}/markers-reads/{{marker}}/merged.fasta",
    output:
        uniq_fasta = f"{outdir}/samples/{{sample}}/markers-reads/{{marker}}/final_merged.fasta",
    log:
        vsearch_unique_log = f"{outdir}/samples/{{sample}}/markers-reads/{{marker}}/vsearch_unique.log"
    shell:
        """
        vsearch --fastx_uniques {input.merged_fasta} --fastaout {output.uniq_fasta} --fasta_width 0 --sizeout --sample {wildcards.sample} > {log.vsearch_unique_log} 2>&1
        """

rule unoise_denoise:
    group: "merge_filter_unique_{sample}"
    input:
        merged_fasta = f"{outdir}/samples/{{sample}}/markers-reads/{{marker}}/final_merged.fasta",
    output:
        denoised_fasta = f"{outdir}/samples/{{sample}}/markers-reads/{{marker}}/final_merged.denoised.fasta",
    log:
        unoise_log = f"{outdir}/samples/{{sample}}/markers-reads/{{marker}}/unoise.log"
    params:
        minsize = config.get("unoise_minsize", 4),  # minimum size of clusters to be kept
        minsize_ratio = config.get('unoise_minsize_ratio', 1), # final minsize if minsize_ratio < 1 = max(minsize_ratio * total size, minsize)
        unoise_alpha = config.get("unoise_alpha", 2),  # alpha parameter for unoise algorithm
    run:
        import math
        def get_size_from_fasta(fasta):
            with open(fasta, 'r') as f:
                for line in f:
                    if line.startswith(">"):
                        match = re.search(r"size=(\d+)", line)
                        if match:
                            yield int(match.group(1))
        if params.minsize_ratio < 1:
            total_read_in_fa = sum([read_size for read_size in get_size_from_fasta(input.merged_fasta)])
            unoise_minsize = max( math.floor(total_read_in_fa * params.minsize_ratio), params.minsize)
            print(f"Dynamic unoise minsize for {input.merged_fasta}: {unoise_minsize}")
        else:
            unoise_minsize = params.minsize
        
        shell(f"""
        vsearch --cluster_unoise {input.merged_fasta} --centroids {output.denoised_fasta} \
        --fasta_width 0  --minsize {unoise_minsize} --unoise_alpha {params.unoise_alpha} \
        --threads 1 > {log.unoise_log} 2>&1
        """)

rule uchime_denovo:
    group: "merge_filter_unique_{sample}"
    input:
        denoised_fasta = f"{outdir}/samples/{{sample}}/markers-reads/{{marker}}/final_merged.denoised.fasta",
    output:
        uchime_fasta = f"{outdir}/samples/{{sample}}/markers-reads/{{marker}}/final_merged.denoised_nochimera.fasta",
    log:
        uchime_log = f"{outdir}/samples/{{sample}}/markers-reads/{{marker}}/uchime.log"
    shell:
        """
        vsearch --uchime3_denovo {input.denoised_fasta} --nonchimeras {output.uchime_fasta} \
        --fasta_width 0 --qmask none --threads 1 > {log.uchime_log} 2>&1
        """

#   def aggregate_input(wildcards):
#       checkpoint_output = checkpoints.somestep.get(**wildcards).output[0]
#       return expand("my_directory/{i}.txt",
#                   i=glob_wildcards(os.path.join(checkpoint_output, "{i}.txt")).i)

# def aggregate_input(wildcards):
#     checkpoint_output = checkpoints.clustering.get(**wildcards).output[0]
#     return expand("post/{sample}/{i}.txt",
#            sample=wildcards.sample,
#            i=glob_wildcards(os.path.join(checkpoint_output, "{i}.txt")).i)

# def get_all_marker_input(wildcards):
#     output_dir = checkpoints.bam_to_marker_fastq.get(sample=wildcards.sample).output[0]
#     all_markers = glob_wildcards(os.path.join(output_dir, "{marker}", "target_R1.fastq.gz")).marker
#     return expand(f"{outdir}/samples/{{sample}}/markers-reads/{{marker}}/final_merged.denoised_nochimera.fasta",
#                   sample=IDs),
#     return expand(f"{outdir}/samples/{wildcards.sample}/markers-reads/{{marker}}/seqtab.tsv", marker=all_markers)

def get_same_marker_all_sample(wildcards):
    all_markers = []
    for sample in IDs:
        output_dir = checkpoints.bam_to_marker_fastq.get(sample=sample).output[0]
        all_markers.extend(glob_wildcards(os.path.join(output_dir, "{marker}", "target_R1.fastq.gz")).marker)
    all_markers = list(set(all_markers))
    return expand(f"{outdir}/samples/{{sample}}/markers-reads/{wildcards.marker}/final_merged.denoised_nochimera.fasta",
                  sample=IDs)


rule merge_per_marker:
    input:
        # expand(f"{outdir}/samples/{{sample}}/markers-reads/complete.flag", sample=IDs),
        expand(f"{outdir}/samples/{{sample}}/markers-reads/{{{{marker}}}}/final_merged.denoised_nochimera.fasta", sample=IDs),
    output:
        merged_denoise_nochim_nonuniq = temp(f"{outdir}/malamp/merged/final_merged_{{marker}}.denoised_nochimera.temp.fasta"),
        merged_denoise_nochim_unsorted = temp(f"{outdir}/malamp/merged/final_merged_{{marker}}.denoised_nochimera.temp1.fasta"),
        merged_denoise_nochim = f"{outdir}/malamp/merged/final_merged_{{marker}}.denoised_nochimera.fasta",
    shell:
        """
        cat {input} > {output.merged_denoise_nochim_nonuniq} &&
        vsearch --fastx_uniques {output.merged_denoise_nochim_nonuniq} --fastaout {output.merged_denoise_nochim_unsorted} --fasta_width 0 --sizein --sizeout --relabel_sha1 &&
        vsearch --sortbysize {output.merged_denoise_nochim_unsorted} --output {output.merged_denoise_nochim}
        """


rule indv_otutab_per_marker:
    input:
        fasta = f"{outdir}/malamp/merged/final_merged_{{marker}}.denoised_nochimera.fasta",
        sample_fa = f"{outdir}/samples/{{sample}}/markers-reads/{{marker}}/final_merged.fasta",
    output:
        otutab = f"{outdir}/samples/{{sample}}/markers-reads/{{marker}}/seqtab.tsv",
    params:
        sequence_to_denoise_criteria = config.get('sequence_to_denoise_criteria', "--id 0.98 --iddef 1 --target_cov 0.98 --query_cov 0.98")
    shell:
        """
        vsearch --usearch_global {input.sample_fa} --db {input.fasta} {params.sequence_to_denoise_criteria} \
        --otutabout {output.otutab} \
        --maxaccepts 0 --maxrejects 0 --maxhits 1 \
        --sizein --sizeout --fasta_width 0 --qmask none --dbmask none --threads {threads}
        """

# def get_per_sample_marker_input(wildcards):
#     all_markers = []
#     for sample in IDs:
#         output_dir = checkpoints.bam_to_marker_fastq.get(sample=sample).output[0]
#         all_markers.extend(glob_wildcards(os.path.join(output_dir, "{marker}", "target_R1.fastq.gz")).marker)
#     all_markers = list(set(all_markers))
#     return expand(f"{outdir}/samples/{wildcards.sample}/markers-reads/{{marker}}/seqtab.tsv", marker=all_markers)


rule merge_otutab_sample:
    input:
        expand(f"{outdir}/samples/{{{{sample}}}}/markers-reads/{{marker}}/seqtab.tsv", marker = Markers)
    output:
        otutab = f"{outdir}/samples/{{sample}}/markers-reads/seqtab.tsv",
    run:
        import pandas as pd
        all_dfs = [pd.read_table(f, sep="\t", header=0) for f in input]
        to_concat = [df for df in all_dfs if df.shape[0] > 0]
        if len(to_concat) == 0:
            # All empty
            otutab = pd.DataFrame(columns=["#OTU ID", wildcards.sample])
        else:
            otutab = pd.concat(to_concat, axis=0, ignore_index=True).fillna(0)
            otutab.columns = ["#OTU ID", wildcards.sample]
            otutab.iloc[:, 1] = otutab.iloc[:, 1].astype(int)
            otutab.to_csv(output.otutab, sep="\t", header=True, index=None)

rule merge_otutab_joint:
    input:
        otutabs = expand(f"{outdir}/samples/{{sample}}/markers-reads/seqtab.tsv", sample=IDs),
        denoise_fa = expand(f"{outdir}/malamp/merged/final_merged_{{marker}}.denoised_nochimera.fasta", marker = Markers)
    output:
        otutab = f"{outdir}/malamp/fastp/seqtab.tsv",
        filtered_fasta = f"{outdir}/malamp/fastp/dereplicated_counted_no_chimeras-unique.fa",
    shell:
        """
        cat {input.denoise_fa} > {output.filtered_fasta} &&
        python {microhaps_basedir}/scripts/merge_otutab.py \
         --denoised_fasta {output.filtered_fasta} \
         --out_seqtab {output.otutab} \
         --indv_seqtabs {input.otutabs}
        """
