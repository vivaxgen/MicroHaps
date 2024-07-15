import os
import pathlib
from ngs_pipeline import cerr, fileutils

microhaps_basedir = os.environ['MICROHAPS_BASEDIR']

reference = f'{microhaps_basedir}/{config.get('mit_reference', "refs/MIT/mit_synthetic_genome.fasta")}'
bedfile = f'{microhaps_basedir}/{config.get('mit_bedfile', "refs/MIT/mit_synthetic_genome.bed")}'
ref_spec = f'{microhaps_basedir}/{config.get('mit_ref_spec', "refs/MIT/mit_spec2.fasta")}'



out_dir = config['outdir']

in_dir = ''
singleton = config.get('singleton', None)
paired = config.get('paired', None)
if singleton:
    read_mode = fileutils.ReadMode.SINGLETON
elif paired:
    read_mode = fileutils.ReadMode.PAIRED_END
else:
    read_mode = None

read_files = fileutils.ReadFileDict(config['infiles'],
                                    underscore=config['underscore'],
                                    mode=read_mode)
IDs = read_files.keys()

def get_read_file(wildcards):
    return read_files._d[wildcards.sample][0]

rule all:
    input:
        f"{out_dir}/bedstats.txt",
        *[f"{out_dir}/analysis/{sample}/mapped/{sample}.bam" for sample in IDs],
        f"{out_dir}/combined_consensus.msa.fasta",
        f"{out_dir}/combined_pairwise_distance.tsv",
        f"{out_dir}/final_report.tsv",

rule concat_bedcov_spec:
    input:
        bedstats = f"{out_dir}/bedstats.txt",
        pairwise_out = f"{out_dir}/combined_pairwise_distance.tsv"
    output:
        final_result = f"{out_dir}/final_report.tsv"
    run:
        import pandas as pd
        bedstats = pd.read_csv(input.bedstats, sep="\t")
        pairwise = pd.read_csv(input.pairwise_out, sep="\t")
        result = pd.merge(bedstats, pairwise, on="sample", how = "left")
        result.to_csv(output.final_result, sep="\t", index=False)

checkpoint sample_pass_criteria:
    input:
        bedstat = f"{out_dir}/bedstats.txt",
    output:
        proceed_samples = temp(f"{out_dir}/proceed_samples.txt"),
        outdir2 = directory(f"{out_dir}/secondary")
    run:
        def sample_pass_criteria():
            import pandas as pd
            try:
                if os.path.exists(f"{out_dir}/bedstats.txt"):
                    data = pd.read_csv(f"{out_dir}/bedstats.txt", sep="\t")
                    result = data.query("depth >= 1 & cov > 0.9")['sample'].tolist()
                    return(result)
            except:
                return([])

        with open(output.proceed_samples, "w") as f:
            for sample in sample_pass_criteria():
                shell(f"mkdir -p {output.outdir2}")
                shell(f"ln -sr {out_dir}/analysis/{sample} {output.outdir2}/{sample}")
                f.write(f"{sample}\n")

def aggregate_input_fasta(wildcards):
    checkpoint_output = checkpoints.sample_pass_criteria.get(**wildcards).output[0]
    return expand(f"{out_dir}/secondary/{{sample}}/{{sample}}.consensus.updated.fasta", sample=open(checkpoint_output).read().strip().split("\n"))

def aggregate_input_pairwise_distance(wildcards):
    checkpoint_output = checkpoints.sample_pass_criteria.get(**wildcards).output[0]
    return expand(f"{out_dir}/secondary/{{sample}}/{{sample}}_pairwise_distance.tsv", sample=open(checkpoint_output).read().strip().split("\n"))

rule pairwise_most_similar:
    input:
        sample_consensus =  f"{out_dir}/secondary/{{sample}}/{{sample}}.consensus.updated.fasta",
        spec_fasta = ref_spec,
    output:
        result_file = f"{out_dir}/secondary/{{sample}}/{{sample}}_pairwise_distance.tsv"
    run:
        from Bio import Align
        from Bio import SeqIO
        import pandas as pd

        aligner = Align.PairwiseAligner()

        sample_sequence = next(SeqIO.parse(open(input.sample_consensus),'fasta'))
        spec_sequence = SeqIO.parse(open(input.spec_fasta),'fasta')

        score = dict()
        score['sample'] = sample_sequence.id
        highest_score = 0
        for species in spec_sequence:
            score[species.id] = round(aligner.score(sample_sequence.seq, species.seq)/len(species.seq),3)
            if score[species.id] > highest_score:
                highest_score = score[species.id]
        
        score["species"] = [k for k,v in score.items() if v == highest_score][0]        
        score = {k:("*"+str(v) if v == highest_score else str(v)) for k,v in score.items() }

        df = pd.DataFrame.from_dict([score])
        df.to_csv(output.result_file, sep="\t", index=False)

rule combine_pair_most_similar:
    input:
        aggregate_input_pairwise_distance,
    output:
        combined_most_similar = f"{out_dir}/combined_pairwise_distance.tsv"
    run:
        import pandas as pd
        df = pd.concat([pd.read_csv(f, sep="\t") for f in input])
        df.to_csv(output.combined_most_similar, sep="\t", index=False)
    
rule combine_consensus:
    input:
        aggregate_input_fasta,
        spec_fasta = ref_spec,
        ref = reference,
    output:
        combined_consensus = f"{out_dir}/combined_consensus.fasta"
    shell:
        "cat {input} >> {output.combined_consensus}"

rule msa_consensus:
    input:
        combined_consensus = f"{out_dir}/combined_consensus.fasta"
    output:
        msa = f"{out_dir}/combined_consensus.msa.fasta"
    shell:
        "muscle -in {input.combined_consensus} -out {output.msa}"

rule get_bedcov_mean:
    input:
        bam = f"{out_dir}/analysis/{{sample}}/mapped/{{sample}}.bam",
    output:
        bedstats = f"{out_dir}/analysis/{{sample}}/mapped/bedstats.txt",
    params:
        bed = bedfile,
        samplename = lambda w: w.sample,
    shell:
        "printf '{params.samplename}\t' > {output.bedstats} && "
        "bedtools coverage -a {params.bed} -b {input.bam} >> {output.bedstats}"

rule join_bedcov_mean:
    input:
        expand(f"{out_dir}/analysis/{{sample}}/mapped/bedstats.txt", sample=IDs)
    output:
        bedstats = f"{out_dir}/bedstats.txt"
    shell:
        "printf 'sample\tchrom\tstart\tend\tdepth\trlen\tqlen\tcov\n' > {output.bedstats} &&"
        "cat {input} >> {output.bedstats}"

rule update_consensus_sample_name:
    input:
        consensus = f"{out_dir}/secondary/{{sample}}/{{sample}}.consensus.fasta"
    output:
        updated_consensus = f"{out_dir}/secondary/{{sample}}/{{sample}}.consensus.updated.fasta"
    params:
        sample = lambda w: w.sample
    run:
        with open(input.consensus, "r") as f:
            lines = f.readlines()
        with open(output.updated_consensus, "w") as f:
            for line in lines:
                if line.startswith(">"):
                    f.write(f">{params.sample}\n")
                else:
                    f.write(line)

rule get_consensus_sequence:
    input:
        bam = f"{out_dir}/secondary/{{sample}}/mapped/{{sample}}.bam"
    output:
        consensus = f"{out_dir}/secondary/{{sample}}/{{sample}}.consensus.fasta"
    params:
        prefix = f"{out_dir}/secondary/{{sample}}/{{sample}}.consensus"
    shell:
        "samtools mpileup -aa -A -d 0 -Q 0 {input.bam} | ivar consensus -q 10 -m 1 -p {params.prefix} && "
        "mv {params.prefix}.fa {output.consensus}"



rule map_reads:
    input:
        lambda w: get_read_file(w),
    params:
        ref = reference,
        minimap2_extra_param = "--sam-hit-only -A2 -B4 -O4,24 -E2,1 --end-bonus 10" # -A2 -B4 -O4,24 -E2,1 --end-bonus 10 
    output:
        bam = f"{out_dir}/analysis/{{sample}}/mapped/{{sample}}.bam"
    shell:
        "minimap2 -ax sr {params.minimap2_extra_param} {params.ref} {input} | samtools view -bSh -f 3 | samtools sort -o {output.bam}"
