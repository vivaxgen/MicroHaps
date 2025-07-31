if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Extract sequence from DADA2 output.")
    parser.add_argument("-t", "--table", help="Path to the seqtab.tsv file from DADA2.", required=True)
    parser.add_argument("-o", "--output_fasta", help="Path to save the output FASTA file.", required=True)

    args = parser.parse_args()

    import pandas as pd
    seqtab_df = pd.read_table(args.table, index_col=0)  # seqtab files constraint: first column is index

    sample_id = seqtab_df.index
    sequences = seqtab_df.columns

    def write_to_fasta(sequences, output_file):
        with open(output_file, "w") as fasta_file:
            for i, seq in enumerate(sequences):
                fasta_file.write(f">col_{i}\n{seq}\n")

    write_to_fasta(sequences, output_file=args.output_fasta)