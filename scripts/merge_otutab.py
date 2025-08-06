
if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Convert merged read-pairs from different samples to a denoised sequence table.")
    parser.add_argument("--denoised_fasta", help="Fasta file with denoised sequences", required=True)
    parser.add_argument("--out_seqtab", help="Output sequence table file", required=True)
    parser.add_argument("--indv_seqtabs", nargs="+", help="OTU table in tsv files", required=True)
    args = parser.parse_args()
    
    denoised_fasta = args.denoised_fasta
    out_seqtab = args.out_seqtab
    indv_seqtabs = args.indv_seqtabs

    from seq_utils import read_fasta
    import re
    import pandas as pd
        
    dfs = [pd.read_csv(f, sep="\t", index_col=0) for f in indv_seqtabs]
    merged = pd.concat(dfs, join = "outer", axis = 1).fillna(0).astype(int).reset_index()

    denoised_seqs = read_fasta(denoised_fasta)

    def parse_header(header):
        matches = re.match(r"^([a-zA-Z0-9]{40});size=(\d+)$", header)
        if not matches:
            raise ValueError(f"Header {header} does not match expected format")
        return matches.group(1), int(matches.group(2))

    def denoised_seq_to_tab(dseq):
        ks, v = zip(*dseq.items())
        oids, expected_counts = zip(*[parse_header(k) for k in ks])
        df = pd.DataFrame(
            {
                "#OTU ID": oids,
                "Expected count": expected_counts,
                "Sequence": v
            }
        )
        return df
    
    denoised_df = denoised_seq_to_tab(denoised_seqs)

    new_merged = pd.merge(merged, denoised_df, on="#OTU ID", how="outer", suffixes=("", "_denoised"))

    # sanity check
    assert new_merged["#OTU ID"].is_unique, "Merged OTU IDs are not unique"
    assert denoised_df["#OTU ID"].isin(merged["#OTU ID"]).all(), "Not all denoised OTU IDs are in the merged table"

    new_merged = new_merged.drop(columns=["#OTU ID"])
    new_merged = new_merged.set_index("Sequence")
    new_merged.index = new_merged.index.str.upper()

    samples = [col for col in new_merged.columns if col != "Expected count"]

    assert (new_merged.loc[:, samples].sum(axis = 1) >= new_merged.loc[:, "Expected count"]).all(), "Expected counts do not at least equate to the sum of counts across samples"

    to_write = new_merged.loc[:, samples].astype(int).T
    to_write.to_csv(out_seqtab, sep="\t", index=True, header=True)