
from seq_utils import * 

if __name__ == "__main__":
    import argparse
    from ngs_pipeline import cerr
    from typing import Literal

    parser = argparse.ArgumentParser(description="Post-process DADA2 output.")
    parser.add_argument("-p", "--paf", help="Path to the PAF mapping the haplotypes to reference.", required=True)
    parser.add_argument("-s", "--seqtab", help="Path to the seqtab.tsv file from DADA2.", required=True)
    parser.add_argument("-f", "--fasta", help="Path to the haplotype FASTA file.", required=True)
    parser.add_argument("-i", "--insert", help="Path to the insert sequence FASTA file.", required=True)
    parser.add_argument("-o", "--output", help="Path to save the output TSV file.", required=True)
    parser.add_argument("--rename_columns_by_id", action='store_true', default=False, help="Rename columns by ID, do not check for identical sequence.")
    args = parser.parse_args()

    import pandas as pd
    seqtab_df = pd.read_table(args.seqtab, index_col=0) # seqtab files constraint: first column is index
    paf_df = parse_paf_file(args.paf)
    fasta_seqs = read_fasta(args.fasta)
    insert_seqs = read_fasta(args.insert)

    assert 'cs' in paf_df.columns, "PAF file must contain 'cs' column"
    
    # sanity check
    # log sequence mapped to multiple reference sequences & retain only the one with the highest match count
    paf_df = paf_df.sort_values(by=['qname', 'nmatch'], ascending=[True, False])
    dup_qname = (paf_df.value_counts('qname') > 1)
    dup_qname_mask = paf_df['qname'].isin(dup_qname[dup_qname].index)
    if dup_qname_mask.any():
        cerr("Warning: Some sequences are mapped to multiple reference sequences. Retaining the one with the highest match count.")
        cerr(f"Duplicated sequences:\n{paf_df[dup_qname_mask]}")
    paf_df = paf_df.drop_duplicates(subset='qname', keep='first')

    # each sequence is only mapped to a single reference sequence
    assert paf_df['qname'].is_unique, "Each sequence must be mapped to only one reference sequence"
    
    # non sequence not mapped to any reference sequence
    # warn and discard these sequences
    if paf_df.query("tname == '*'").shape[0] > 0:
        cerr("Warning: Some sequences are not mapped to any reference sequence. These will be discarded.")
        cerr(f"Sequences not mapped to any reference:\n{paf_df.query('tname == \'*\'')}")
        paf_df = paf_df.query("tname != '*'")

    assert paf_df.query("tname == '*'").shape[0] == 0, "Some sequences are not mapped to any reference sequence"

    def add_unmapped_qseq_to_cs(row, get: Literal["start", "end"]):
        if row['orientation'] == '-':
            rel_t_seq = reverse_complement(insert_seqs[row['tname']])
            rel_q_seq = reverse_complement(fasta_seqs[row['qname']])
        else:
            rel_t_seq = insert_seqs[row['tname']]
            rel_q_seq = fasta_seqs[row['qname']]
        if get == "start":
            start_s = ""
            if row['qstart'] > 0: # there are bases before the mapped region
                start_s += f"+{rel_q_seq[:row['qstart']]}"
            if row['tstart'] > 0: # there are bases before the mapped region in the inserts but not in the query
                start_s += f"-{rel_t_seq[:row['tstart']]}"
            return start_s
        if get == "end":
            end_s = ""
            if row['qend'] < row['qlen']: # there are bases after the mapped region
                end_s += f"+{rel_q_seq[row['qend']:]}"
            # At this stage, all the query sequences are used up
            if row['tend'] < row['tlen']: # there are bases after the mapped region in the inserts but not in the query
                end_s += f"-{rel_t_seq[row['tend']:]}"
            return end_s


    # add start to final representation
    paf_df["non_mapped_start"] = paf_df.apply(
        lambda row: add_unmapped_qseq_to_cs(row, get="start"), axis=1
    )

    # add end to final representation
    paf_df["non_mapped_end"] = paf_df.apply(
        lambda row: add_unmapped_qseq_to_cs(row, get="end"), axis=1
    )

    paf_df["cs"] = paf_df["non_mapped_start"] + paf_df["cs"] + paf_df["non_mapped_end"]

    # 0 == tstart because the query is mapped to the start of the target sequence
    # if already mapped to the start of target, the deletion/insertion is actual
    # if not, deletion/insertion is the "artifact"
    paf_df.loc[:, "start_s"] = paf_df.apply(lambda row: "s" if 0 == row["tstart"] else f'{row["qstart"]}:{row["tstart"]}', axis=1)
    
    # tlen == tend because if the whole query sequence is mapped to the end
    # any extra sequence will be actual insertion
    # if less, then deletion/insertion is the "artifact"
    paf_df.loc[:, "end_s"] = paf_df.apply(lambda row: "e" if (row["tlen"] == row["tend"]) else f'{row["qend"]}:{row["tend"]}', axis=1)


    paf_df.loc[:, "final_representation"] = paf_df.loc[:, "tname"] + "," + paf_df.loc[:, "cs"] + "_" + paf_df.loc[:, "start_s"].astype(str) + "-" + paf_df.loc[:, "end_s"].astype(str)

    assert paf_df.query("orientation == '+'")['final_representation'].is_unique, "Final representation for the same orientation must be unique"
    assert paf_df.query("orientation == '-'")['final_representation'].is_unique, "Final representation for the same orientation must be unique"

    seq_mapping = pd.DataFrame.from_dict(fasta_seqs, orient = 'index', columns = ['sequence'])
    seq_mapping = seq_mapping.merge(paf_df[['qname', 'final_representation']], left_index=True, right_on='qname', how='left')
    # drop columns that cannot be converted
    to_drop = seq_mapping[seq_mapping['final_representation'].isna()]['sequence'].tolist()
    if to_drop:
        cerr(f"Warning: The following sequences are not mapped to any reference sequence and will be dropped: {to_drop}")
        seqtab_df = seqtab_df.drop(columns=to_drop)
        seq_mapping = seq_mapping.dropna(subset=['final_representation'])

    if args.rename_columns_by_id:
        seq_mapping["col"] = seq_mapping["qname"].str.replace("col_", "").astype(int)
        seq_mapping = seq_mapping.sort_values("col")
        assert all(seq_mapping["col"] == list(range(seqtab_df.shape[1])))
        seqtab_df_renamed = seqtab_df.copy()
        seqtab_df_renamed.columns = seq_mapping["final_representation"].values
    else:
        rename_mapping = seq_mapping.set_index('sequence')['final_representation'].to_dict()
        seqtab_df_renamed = seqtab_df.rename(columns=rename_mapping)
    column_order = ["sample"] + list(seqtab_df_renamed.columns)
    seqtab_df_renamed.index.name = 'sample'
    seqtab_df_renamed.reset_index(inplace=True)
    seqtab_df_renamed[column_order].to_csv(args.output, sep="\t", index=False, header=True)

