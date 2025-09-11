dada2_bbmapmerge_bbmap = config.get('merge_map') # dada2, bbmap_merge or bbmerge

rule extract_sequences_from_seqtab:
    input:
        f"{outdir}/malamp/{dada2_bbmapmerge_bbmap}/seqtab.tsv"
    output:
        fasta = f"{outdir}/malamp/haplotypes.fasta",
    shell:
        "python {microhaps_basedir}/scripts/seqtab_to_fasta.py --table {input[0]} --output_fasta {output.fasta}"

rule msa_align_haplotype_to_reference:
    threads: 10
    input:
        fasta = f"{outdir}/malamp/haplotypes.fasta",
        paf = f"{outdir}/malamp/haplotypes.paf",
        reference = insertseq_file,
        seqtab = f"{outdir}/malamp/{dada2_bbmapmerge_bbmap}/seqtab.tsv",
        reference_idx = f"{insertseq_file.removesuffix('.fasta')}.1.bt2",
    output:
        outHaplotype_rm_ins = f"{outdir}/malamp/outputHaplotypes_rm_ins.tsv",
    params:
        script_dir = f"{microhaps_basedir}/scripts",
        cs_style = "long" if config.get("post_process", "cs_long") == "cs_long" else "short"
    run:
        from concurrent.futures import ThreadPoolExecutor
        import pandas as pd
        import sys

        sys.path.append(params.script_dir)
        import seq_utils

        paf_df = seq_utils.parse_paf_file(input.paf)
        haplotypes_seq = seq_utils.read_fasta(input.fasta)
        reference_seq = seq_utils.read_fasta(input.reference)
        seqtab_df = pd.read_table(input.seqtab, index_col=0)

        reference_fasta = []
        query_fasta = []
        output_fasta = []
        log_file = []
        all_markers_name = paf_df['tname'].unique()
        for marker in all_markers_name:
            haplotype_seq_names = paf_df.query("tname == @marker").qname.to_list()
            haplotype_seqs = {a: haplotypes_seq[a] for a in haplotype_seq_names}
            need_to_rc = paf_df.query("tname == @marker & orientation == '-'").qname.to_list()
            if len(need_to_rc) > 0:
                for n in need_to_rc:
                    haplotype_seqs[n] = seq_utils.reverse_complement(haplotype_seqs[n])
            seq_utils.write_fasta(haplotype_seqs, f"{outdir}/malamp/{marker}_haplotypes.fasta")
            seq_utils.write_fasta({marker: reference_seq[marker]}, f"{outdir}/malamp/{marker}_reference.fasta")
            reference_fasta.append(f"{outdir}/malamp/{marker}_reference.fasta")
            query_fasta.append(f"{outdir}/malamp/{marker}_haplotypes.fasta")
            output_fasta.append(f"{outdir}/malamp/{marker}_aligned.fasta")
            log_file.append(f"{outdir}/logs/mafft/{marker}_mafft.log")

        shell(f"mkdir -p {outdir}/logs/mafft")
        marker_per_mafft = 2
        # --threadit = 0 because The iterative refinement step can return different results for different runs, when using mutitple threads.  To obtain the same result for every run, add the --threadit 0 flag, which disables multithreading only for iterative refinement step. (added 2020/Dec)
        def run_mafft(ref_fa, query_fa, output_fa, log_file, marker_per_mafft = 2):
            shell(f"mafft-linsi --maxiterate 1000 --thread {marker_per_mafft} --threadit 0 --preservecase --keeplength --add {query_fa} {ref_fa} > {output_fa} 2> {log_file}")
            inf = seq_utils.read_fasta(output_fa, ordered = True)
            seq_utils.write_fasta(inf[1:], output_fa)
            return 1

        nworker = int(threads/marker_per_mafft) if int(threads/marker_per_mafft) > 0 else 1
        with ThreadPoolExecutor(max_workers=nworker) as executor:
            result = list(executor.map(run_mafft, reference_fasta, query_fasta, output_fasta, log_file,
                [marker_per_mafft]*len(reference_fasta)))
        
        if not all([res == 1 for res in result]):
            raise Exception("Error in aligning some markers, please check previous logs")

        joined_fasta = f"{outdir}/malamp/haplotypes_rm_ins_aligned.fasta"
        no_ins_paf = f"{outdir}/malamp/haplotypes_rm_ins.paf"
        shell(f"cat {' '.join(output_fasta)} > {joined_fasta}")

        output_sam = f"{outdir}/malamp/out.sam"
        shell(f"bowtie2 -f --end-to-end -x {input.reference.removesuffix('.fasta')} -U {joined_fasta} --mp 1,0 --rdg 2,1 --sam-opt-config 'md' -S {output_sam} 2> {outdir}/logs/align_haplotypes_to_insertseq.log")

        

        def alignment_to_cs_tags(sam, cs_style = "long"):
            import pysam
            import cstag
            
            long_cs_style = True if cs_style == "long" else False

            algs = pysam.AlignmentFile(sam, "r")
            qname_cstag = []
            for alg in algs:
                if not alg.has_tag("MD"):
                    raise ValueError("Missing MD tag")

                cs_tag = cstag.call(cigar = alg.cigarstring, md = alg.get_tag("MD"), seq = alg.query_sequence, long = long_cs_style)
                qname = int(alg.query_name.replace("col_", ""))
                qname_cstag.append({
                    "qname": qname, "cs_tag": f"{alg.reference_name},{cs_tag}"
                })
            qname_cstag = pd.DataFrame(qname_cstag)
            qname_cstag = qname_cstag.sort_values("qname")
            return qname_cstag
        
        qname_cstag = alignment_to_cs_tags(output_sam, params.cs_style)

        #check if any columns in seqtab_df that is not in qname_cstag

        def transform_seqtab(seqtab, rename_df):
            to_drop = [colid for colid in list(range(seqtab.shape[1])) if not colid in rename_df['qname']]
            to_drop_col = seqtab.columns[to_drop]

            if not to_drop_col.empty:
                print(f"Warning: The following sequences are not mapped to any reference sequence and will be dropped: {to_drop_col}", file=sys.stderr)
                seqtab = seqtab.drop(columns=to_drop_col)
        
            seqtab_df_renamed = seqtab.copy()
            seqtab_df_renamed.columns = rename_df["cs_tag"].values
            
            # check for same cstag without insertions
            for dup_column in seqtab_df_renamed.columns[seqtab_df_renamed.columns.duplicated(keep = 'first')]:
                series = seqtab_df_renamed[[dup_column]].sum(axis = 1)
                seqtab_df_renamed = seqtab_df_renamed.drop(columns=dup_column)
                seqtab_df_renamed[dup_column] = series

            column_order = ["sample"] + list(seqtab_df_renamed.columns)
            seqtab_df_renamed.index.name = 'sample'
            seqtab_df_renamed.reset_index(inplace=True)
            return seqtab_df_renamed[column_order]

        transformed_seqtab = transform_seqtab(seqtab_df, qname_cstag)
        transformed_seqtab.to_csv(output.outHaplotype_rm_ins, sep="\t", index=False, header=True)

        temp_files = reference_fasta + query_fasta + output_fasta

        shell(f"rm {' '.join(temp_files)}")

        # shell(f"minimap2 -x asm20 --secondary=no --end-bonus --cs={params.cs_style} {input.reference} {joined_fasta} --paf-no-hit -o {no_ins_paf}")

        # shell(f"""python {microhaps_basedir}/scripts/post_process_dada2.py \
        #     --paf {no_ins_paf} \
        #     --fasta {joined_fasta} \
        #     --seqtab {input.seqtab} \
        #     --output {output.outHaplotype_rm_ins} \
        #     --insert {input.reference} \
        #     --rename_columns_by_id""")

        # print(f"""python {microhaps_basedir}/scripts/post_process_dada2.py \
        #     --paf {no_ins_paf} \
        #     --fasta {joined_fasta} \
        #     --seqtab {input.seqtab} \
        #     --output {output.outHaplotype_rm_ins} \
        #     --insert {input.reference} \
        #     --rename_columns_by_id""")


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
        "minimap2 -x sr -t {threads} {params.minimap_params} --secondary=no --cs={params.cs_style} {params.reference} {input.fasta} --paf-no-hit -o {output.paf}"

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


