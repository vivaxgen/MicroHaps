
import argparse

p = argparse.ArgumentParser()
p.add_argument('--variant_list', type=str, required=True, help='File containing variant positions')
p.add_argument('--reference', type=str, required=True, help='Reference sequence file')
p.add_argument('--pseudohaplotypes', type=str, required=True, help='Pseudohaplotypes file')
p.add_argument('--gff_file', type=str, required=True, help='GFF file')
p.add_argument('--protein_predictions', type=str, required=True, help='Protein predictions file')
p.add_argument('--drug_resistance_report', type=str, required=True, help='Drug resistance report file')
p.add_argument('--min_support', type=int, default=5, help='Minimum support to consider a variant for drug resistance report')
p.add_argument('--sample', type=str, default="sample", help='Sample name to include in the drug resistance report')
args = p.parse_args()


def predict_aa_change(variant_list, pseudohaplotypes, reference, gff_file, protein_predictions, drug_resistance_report, sample, min_support = 5):
    import gff_csq
    import pysam
    import pandas as pd

    variants_ = pd.read_table(variant_list, header=None, names=["chrom", "pos0", "pos", "variant_id"])
    pseudohaplotypes_ = pd.read_table(pseudohaplotypes)
    fasta_file = pysam.FastaFile(reference)
    gff_ = gff_csq.GFF3(gff_file)

    ref_lookup = lambda seqid, pos: fasta_file.fetch(seqid, pos-1, pos).upper()

    fin_rows = []
    for _, row in pseudohaplotypes_.iterrows():
        gene_id = row["marker"].split(":")[0]
        aa_pos = int(row["marker"].split(":")[1])
        codons = {pos: val for pos, val in zip(variants_.query("variant_id == @row.marker")["pos"].values, list(row["haplotype"]))}
        aa = gff_csq.GeneCoords(gff_, gene_id).translate_variant(aa_pos, ref_lookup, codons)
        if len(aa) != 1:
            raise ValueError(f"Warning: {row['marker']} has {len(aa)} predicted amino acids, skipping")
        
        aa = aa.popitem()[1]
        updated_row = row.copy()

        updated_row["ref_codon"] = aa.ref_codon
        updated_row["ref_aa"] = aa.ref_aa
        updated_row["alt_codon"] = aa.alt_codon
        updated_row["alt_aa"] = aa.alt_aa
        updated_row["filled_from_ref"] = ",".join([str(a) for a in aa.ref_filled])
        fin_rows.append(updated_row.to_frame())
    if len(fin_rows) > 0:
        fin_df = pd.concat(fin_rows, axis=1).T
        fin_df.to_csv(protein_predictions, sep="\t", index=False)
    else:
        fin_df = pd.DataFrame(columns=["sample", "marker", "haplotype", "count", "ref_codon", "ref_aa", "alt_codon", "alt_aa", "filled_from_ref"])
        fin_df.to_csv(protein_predictions, sep="\t", index=False)
    
    final_report = []
    for id_ in variants_["variant_id"].unique():
        found_rows = fin_df.query("marker == @id_ and count >= @min_support")
        if len(found_rows) == 0:
            final_report.append({
                "sample": sample,
                "marker": id_,
                "codon": ".",
                "aa": ".",
                "count": 0,
                "filled": ""
            })
        else:
            for _, row in found_rows.iterrows():
                final_report.append({
                    "sample": sample,
                    "marker": id_,
                    "codon": row["alt_codon"],
                    "aa": row["alt_aa"],
                    "count": row["count"],
                    "filled": row["filled_from_ref"]
                })
    final_report_df = pd.DataFrame(final_report)
    final_report_df.to_csv(drug_resistance_report, sep="\t", index=False)

if __name__ == "__main__":
    predict_aa_change(args.variant_list, args.pseudohaplotypes, args.reference, args.gff_file, args.protein_predictions, args.drug_resistance_report, args.sample, args.min_support)

# variants_ = pd.read_table("haplotype_list.bed", header=None, names=["chrom", "pos0", "pos", "variant_id"])
# pseudohaplotypes_ = pd.read_table("results/ET-193/haplotype_pseudohaplotypes.tsv")
# fasta_file = pysam.FastaFile("/home/data/malaria/work/ldwg/LATEST_NGSPL/vvg-MicroHaps/envs/MicroHaps/configs/refs/Pf/Pf3D7_v3/Pf3D7_v3.fasta")
# gff_ = gff_csq.GFF3("/home/data/malaria/work/ldwg/LATEST_NGSPL/vvg-MicroHaps/envs/MicroHaps/configs/refs/Pf/Pf3D7_v3/Pf3D7_v3.gff.bz2")

# variants_["hap_pos"] = variants_.groupby("variant_id")["variant_id"].cumcount()