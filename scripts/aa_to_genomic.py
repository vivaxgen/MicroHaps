
import argparse

p = argparse.ArgumentParser()
p.add_argument('--aa_pos_file', type=str, required=True, help='File containing amino acid positions')
p.add_argument('--reference', type=str, required=True, help='Reference sequence file')
p.add_argument('--gff_file', type=str, required=True, help='GFF file')
p.add_argument('--output', type=str, required=True, help='Output file')
args = p.parse_args()

def convert_aa_pos_to_genomic(aa_pos_file, reference, gff_file, output):
    import gff_csq
    import pysam
    import pandas as pd

    interested_aa_pos = pd.read_table(aa_pos_file)
    fasta_file = pysam.FastaFile(reference)
    gff_ = gff_csq.GFF3(gff_file)

    ref_lookup = lambda seqid, pos: fasta_file.fetch(seqid, pos-1, pos).upper()

    results = []
    for _, row in interested_aa_pos.iterrows():
        gc = gff_csq.GeneCoords(gff_, row["gene_id"]).aa_to_genomic(row["aa_pos"])
        variant_id = f"{row['gene_id']}:{row['aa_pos']}"
        if len(gc) == 0 or len(gc) > 1:
            print(f"Warning: {variant_id} has {len(gc)} genomic coordinates, skipping")
            continue
        genomic_positions = gc.popitem()[1]
        for gpos in genomic_positions:
            chrom = gpos.seqid
            pos0 = gpos.pos - 1 # convert to 0-based
            pos = gpos.pos
            results.append({
                "chrom": chrom,
                "start": pos0,
                "end": pos,
                "variant_id": variant_id
            })
    results_df = pd.DataFrame(results).to_csv(output, sep="\t", index=False, header=False)

if __name__ == "__main__":
    convert_aa_pos_to_genomic(args.aa_pos_file, args.reference, args.gff_file, args.output)
# reference = "/home/data/malaria/work/ldwg/LATEST_NGSPL/vvg-MicroHaps/envs/MicroHaps/configs/refs/Pf/Pf3D7_v3/Pf3D7_v3.fasta"
# gff = "/home/data/malaria/work/ldwg/LATEST_NGSPL/vvg-MicroHaps/envs/MicroHaps/configs/refs/Pf/Pf3D7_v3/Pf3D7_v3.gff.bz2"
# aa_pos_file ="drug_markers_interested_aa_pos.tsv"
# interested_aa_pos = pd.read_table(aa_pos_file)
# fasta_file = pysam.FastaFile(reference)
# gff_ = gff_csq.GFF3(gff)