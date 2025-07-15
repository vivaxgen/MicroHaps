if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Extract sequence from DADA2 output.")
    parser.add_argument("-t", "--table", help="Path to the seqtab.tsv file from DADA2.", required=True)
    parser.add_argument("-o", "--output_fasta", help="Path to save the output FASTA file.", required=True)
    parser.add_argument("-s", "--sample_id", type=str, default="index",
                        help="Column name for sample IDs in the seqtab.tsv file. [default: 'index'] - using the actual index")

    args = parser.parse_args()

    import pandas as pd
    seqtab_df = pd.read_table(args.table)

    if args.sample_id == "index":
        sample_id = seqtab_df.index
        sequences = seqtab_df.columns
    else:
        assert args.sample_id in seqtab_df.columns, f"Column '{args.sample_id}' not found in the seqtab.tsv file."
        sample_id = seqtab_df[args.sample_id]
        sequences = seqtab_df.drop(columns=[args.sample_id]).columns

    def write_to_fasta(sequences, output_file):
        with open(output_file, "w") as fasta_file:
            for i, seq in enumerate(sequences):
                fasta_file.write(f">col_{i}\n{seq}\n")

    write_to_fasta(sequences, output_file=args.output_fasta)