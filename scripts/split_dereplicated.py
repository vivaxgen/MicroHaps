import pandas as pd
import argparse
from seq_utils import read_fasta, read_uc_table, write_fasta


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Split Fasta file into multiple Fasta based on UC table")
    parser.add_argument("-f", help="Fasta to split", required=True)
    parser.add_argument("-u", help="UC table to use for splitting", required=True)
    parser.add_argument("-o", help="Output directory", required=True)
    parser.add_argument("-i", help="Insert sequence file",  required=True)
    args = parser.parse_args()
    
    fasta = args.f
    uc = args.u
    output_dir = args.o
    insert_seq = args.i

    fasta_dict = read_fasta(fasta)
    all_targets = read_fasta(insert_seq).keys()
    uc_table = read_uc_table(uc)

    filtered_table = uc_table.query("label_target != '*'")[['label_query', 'label_target']]
    for target in all_targets:
        target_table = filtered_table[filtered_table['label_target'] == target]
        write_fasta(
            {k: fasta_dict[k] for k in target_table['label_query'] if k in fasta_dict},
            f"{output_dir}/{target}-split.fa"
        )
