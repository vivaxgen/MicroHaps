
rule bam_to_fastq:
    threads: 4
    input:
        bam = "{pfx}/samples/{sample}/maps/final.bam"
    output:
        R1 = "{pfx}/samples/{sample}/mhaps-reads/target_R1.fastq.gz",
        R2 = "{pfx}/samples/{sample}/mhaps-reads/target_R2.fastq.gz"
    shell:
        "samtools collate -u -O {input.bam}"
        " | samtools fastq --thread 2 -1 {output.R1} -2 {output.R2} -0 /dev/null -s /dev/null -n"

# EOF
