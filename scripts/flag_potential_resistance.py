import pandas as pd
import argparse as ap
from itertools import groupby
import gff_csq
import pysam

ap = ap.ArgumentParser(description="Flag potential drug resistance from samples aa call")
ap.add_argument("-i", "--input", required=True, help="input file with aa calls")
ap.add_argument("-o", "--output", required=True, help="output file with potential resistance flagged")
ap.add_argument("-d", "--database", required=True, help="database of known resistance mutations")
ap.add_argument("-m", "--min_depth", type=int, default=10, help="minimum depth to consider a call valid")
ap.add_argument('--reference', type=str, required=True, help='Reference sequence file')
ap.add_argument('--gff_file', type=str, required=True, help='GFF file')

if __name__ == "__main__":
    args = ap.parse_args()
    fasta_file = pysam.FastaFile(args.reference)
    gff_ = gff_csq.GFF3(args.gff_file)
    db = pd.read_table(args.database)
    df = pd.read_table(args.input).set_index("sample")
    
    gene_name_to_gene_id = {"Pfcrt": "PF3D7_0709000",
        "Pfdhfr": "PF3D7_0417200", "Pfdhps": "PF3D7_0810800",
        "Pfk13": "PF3D7_1343700", "Pfmdr1": "PF3D7_0523000"}

    def get_ref_aa(gene_id, aa_pos):
        ref_lookup = lambda seqid, pos: fasta_file.fetch(seqid, pos-1, pos).upper()
        aa = gff_csq.GeneCoords(gff_, gene_id).translate_variant(aa_pos, ref_lookup)
        if len(aa) != 1:
            raise ValueError(f"Warning: {row['marker']} has {len(aa)} predicted amino acids, skipping")
        aa = aa.popitem()[1]
        return aa.ref_aa

    parts = db["Alteration"].str.strip().str.extract(r"^([A-Za-z*])(\d+)([A-Za-z*]?)$")
    unparsed = parts[1].isna()
    if unparsed.any():
        raise ValueError(f"Cannot parse Alteration(s): {db.loc[unparsed, 'Alteration'].tolist()}")
    no_sensitive = parts[2] == ""
    db["interested_aa"] = db["Gene"] + ":" + parts[1]
    db["interested_var"] = parts[2].where(~no_sensitive, parts[0])
    db["gene_id"] = db["Gene"].map(gene_name_to_gene_id)
    db["ref_aa"] = db.apply(lambda row: get_ref_aa(row["gene_id"], int(row["interested_aa"].split(":")[1])), axis=1)

    results = db.copy()
    for i, row in df.iterrows():
        sample = i
        results.loc[:, sample] = "?"
        for db_i, db_row in results.iterrows():
            if db_row["interested_aa"] in row.index:
                if db_row["interested_aa"] not in row.index:
                    continue
                gts, depths = row[db_row["interested_aa"]].split("|")
                gts = gts.split(",")
                depths = list(map(int, depths.split(",")))
                if db_row["interested_var"] in gts:
                    call = "+"
                    depth = depths[gts.index(db_row["interested_var"])]
                else:
                    call = "-"
                    depth = sum(depths)
                if depth < args.min_depth:
                    call += " (low_depth)"
                results.loc[db_i, sample] = call

    results = results.rename(columns={"interested_var": "resistant_aa"})
    results.drop(columns=["interested_aa", "gene_id"]).to_csv(f"{args.output.replace('.tsv', '.details.tsv')}", sep="\t", index=False)

    sample_cols = list(df.index)
    meta_cols = ["Classification", "Drug"]

    def join_by_gene(group, col):
        return "&".join(
            group.groupby("Gene", sort=False)[col]
                .apply(lambda s: "+".join(s.astype(str)))
        )

    def aggregate_group(group):
        row = group[meta_cols].iloc[0].to_dict()
        group = group.sort_values(["Gene", "interested_aa"])
        row["Gene(s)"]         = "&".join(group["Gene"].unique())
        row["Alterations"]     = join_by_gene(group, "Alteration")
        row["resistant_aa(s)"] = join_by_gene(group, "resistant_aa")
        row["ref_aa(s)"]       = join_by_gene(group, "ref_aa")
        for sample in sample_cols:
            values = group[sample]
            observed = sorted({v.split(" ")[0] for v in values})
            low_depth = values.str.contains("low depth").any()
            match observed:
                case _ if "?" in observed:
                    call = "?"
                case ["+"]:
                    call = "+"
                case ["-"]:
                    call = "-"
                case ["+", "-"]:
                    call = "-/+"
                case _:
                    call = "?"
            row[sample] = call + "*" if low_depth else call
        return pd.Series(row)

    results["interested_aa"] = results["interested_aa"].str.split(":").str[1].astype(int)
    grouped = results.groupby("Index").apply(aggregate_group, include_groups=False).reset_index()
    grouped.to_csv(args.output, sep="\t", index=False)