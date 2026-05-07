import pandas as pd
import argparse as ap

ap = ap.ArgumentParser(description="Flag potential drug resistance from samples aa call")
ap.add_argument("-i", "--input", required=True, help="input file with aa calls")
ap.add_argument("-o", "--output", required=True, help="output file with potential resistance flagged")
ap.add_argument("-d", "--database", required=True, help="database of known resistance mutations")
ap.add_argument("-m", "--min_depth", type=int, default=10, help="minimum depth to consider a call valid")

if __name__ == "__main__":
    args = ap.parse_args()
    db = pd.read_table(args.database)
    df = pd.read_table(args.input).set_index("sample")
    
    translate = {'CRT': "Pfcrt", 'DHFR': "Pfdhfr", 'DHPS': "Pfdhps", 'K13': "Pfk13", 'MDR1': "Pfmdr1"}
    translate_rev = {v: k for k, v in translate.items()}
    
    db["interested_aa"] = db["Gene"].map(translate_rev) + ":" + db["Alteration"].str[1:-1]
    db["interested_var"] = db["Alteration"].str[-1]

    results = db.copy()
    for i, row in df.iterrows():
        sample = i
        results.loc[:, sample] = "?"
        for db_i, db_row in results.iterrows():
            if db_row["interested_aa"] in row.index:
                gts, depths = row[db_row["interested_aa"]].split("|")
                gts = gts.split(",")
                depths = list(map(int, depths.split(",")))
                if db_row["interested_var"] in gts:
                    var_index = gts.index(db_row["interested_var"])
                    if depths[var_index] >= args.min_depth:
                         results.loc[db_i, sample] = "+"
                    else:
                        results.loc[db_i, sample] = "+ (low depth)"
                else:
                    results.loc[db_i, sample] = "-"
            else:
                continue

    results.drop(columns=["interested_aa", "interested_var"]).to_csv(f"{args.output.replace('.tsv', '.details.tsv')}", sep="\t", index=False)

    sample_cols = list(df.index)
    meta_cols = ["Classification", "Drug", "Gene"]

    def aggregate_group(group):
        row = group[meta_cols].iloc[0].to_dict()
        row["Alterations"] = ", ".join(group["Alteration"].tolist())
        for sample in sample_cols:
            vals = group[sample].tolist()
            row[sample] = "+" if all(v == "+" for v in vals) else \
                ("-" if all(v == "-" for v in vals) else "?")
        return pd.Series(row)

    grouped = results.groupby("Index").apply(aggregate_group, include_groups=False).reset_index()
    grouped.to_csv(args.output, sep="\t", index=False)