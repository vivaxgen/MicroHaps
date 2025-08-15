rule prepare_shorah_ref:
    input:
        ref = {refseq}
    output:
        ref_dir = expand(f"{outdir}/shorah/{{region}}", region = config.get(regions))
    params:
        chroms = config.get(regions, [])
    run:
        ref_chroms = []
        with open(f"{input.ref}.fai", "rb") as fai:
            for line in fai:
                seq_id, seq_len, offset = line.decode('utf-8').split("\t")[:3]
                if seq_id in params.chroms:
                    ref_chroms.append((seq_id, int(seq_len), int(offset)))

        with open(input.ref, "r") as ref:
            for chrom in ref_chroms:
                seq_id, seq_len, offset = chrom
                ref.seek(offset)
                seq = ref.read(seq_len)
                with open(f"{output.ref_dir}/{seq_id}.fa", "w+") as out_f:
                    out_f.write(f">{seq_id}\n{seq}\n")
        
        # Check that all regions available
        from glob import glob
        
        # index all the reference sequences
        


rule map_for_shorah:
    threads: 2
    input:
        read1 = f"{outdir}/samples/{{sample}}/mhaps-reads/trimmed-filtered_R1.fastq.gz",
        read2 = f"{outdir}/samples/{{sample}}/mhaps-reads/trimmed-filtered_R2.fastq.gz",
    output:
        bam = f"{outdir}/samples/{{sample}}/shorah/mapped.bam",
    params:
        rg = "-R '@RG\tID:{sample}\tSM:{sample}\tLB:LIB-{sample}\tPL:generic'"
    log:
        shorah_map =  f"{outdir}/samples/{{sample}}/logs/shorah_bwa.log",
    shell:
        "bwa-mem2 mem -M -t {threads} {params.rg} {refseq} {input.read1} {input.read2} 2> {log.shorah_map} | samtools sort -o {output.bam}"




SRR28327857,155 <- global cluster

SRR28327720,26 <- global custer