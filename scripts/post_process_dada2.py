
def parse_paf_line(line):
    fields = line.strip().split("\t")
    if len(fields) < 12:
        raise ValueError("PAF line does not contain enough fields")
    
    initial_fields = {
        "qname": fields[0],
        "qlen": int(fields[1]),
        "qstart": int(fields[2]),
        "qend": int(fields[3]),
        "orientation": fields[4],
        "tname": fields[5],
        "tlen": int(fields[6]),
        "tstart": int(fields[7]),
        "tend": int(fields[8]),
        "nmatch": int(fields[9]),
        "nbases": int(fields[10]),
        "mapq": int(fields[11])
    }
    
    if len(fields) > 12:
        for i in range(12, len(fields)):
            tag_name = fields[i][:2]
            tag_dtype = fields[i][3]
            tag_value = fields[i][5:]
            if tag_dtype == "i":
                tag_value = int(tag_value)
            elif tag_dtype == "f":
                tag_value = float(tag_value)
            initial_fields[tag_name] = tag_value

    return initial_fields

def parse_paf_file(paf_file):
    import pandas as pd

    lines = open(paf_file, 'r').read().strip().split("\n")

    paf_data = [parse_paf_line(line) for line in lines]
    return pd.DataFrame(paf_data)

def read_fasta(fasta_file, ordered = False):
    import re
    reference = open(fasta_file, 'r').read().strip()
    ref_groups = re.findall(r">(.*)\n([^>]*)", reference)
    if ordered:
        seqs =[(group[0], group[1].replace("\n", "")) for group in ref_groups]
    else:
        seqs = {group[0]: group[1].replace("\n", "") for group in ref_groups }
    return(seqs)

def reverse_complement(seq):
    complement = str.maketrans("ACGT", "TGCA")
    return seq.translate(complement)[::-1]

if __name__ == "__main__":
    import argparse
    from ngs_pipeline import cerr
    from typing import Literal

    parser = argparse.ArgumentParser(description="Post-process DADA2 output.")
    parser.add_argument("-p", "--paf", help="Path to the PAF mapping the haplotypes to reference.", required=True)
    parser.add_argument("-s", "--seqtab", help="Path to the seqtab.tsv file from DADA2.", required=True)
    parser.add_argument("-f", "--fasta", help="Path to the haplotype FASTA file.", required=True)
    parser.add_argument("-o", "--output", help="Path to save the output TSV file.", required=True)

    args = parser.parse_args()

    import pandas as pd
    seqtab_df = pd.read_table(args.seqtab, index_col=0) # seqtab files constraint: first column is index
    paf_df = parse_paf_file(args.paf)
    fasta_seqs = read_fasta(args.fasta)

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
            rel_q_seq = reverse_complement(fasta_seqs[row['qname']])
        else:
            rel_q_seq = fasta_seqs[row['qname']]
        if get == "start":
            if row['qstart'] > 0:
                return f"+{rel_q_seq[:row['qstart']]}"
            else:
                return ""
        if get == "end":
            if row['qend'] < row['qlen']:
                return f"+{rel_q_seq[row['qend']:]}"
            else:
                return ""

    # add start to final representation
    paf_df["non_mapped_start"] = paf_df.apply(
        lambda row: add_unmapped_qseq_to_cs(row, get="start"), axis=1
    )

    # add end to final representation
    paf_df["non_mapped_end"] = paf_df.apply(
        lambda row: add_unmapped_qseq_to_cs(row, get="end"), axis=1
    )

    paf_df["cs"] = paf_df["non_mapped_start"] + paf_df["cs"] + paf_df["non_mapped_end"]

    paf_df.loc[:, "final_representation"] = paf_df.loc[:, "tname"] + "," + paf_df.loc[:, "cs"]

    assert paf_df['final_representation'].is_unique, "Final representation must be unique"

    seq_mapping = pd.DataFrame.from_dict(fasta_seqs, orient = 'index', columns = ['sequence'])
    seq_mapping = seq_mapping.merge(paf_df[['qname', 'final_representation']], left_index=True, right_on='qname', how='left')
    # drop columns that cannot be converted
    to_drop = seq_mapping[seq_mapping['final_representation'].isna()]['sequence'].tolist()
    if to_drop:
        cerr(f"Warning: The following sequences are not mapped to any reference sequence and will be dropped: {to_drop}")
        seqtab_df = seqtab_df.drop(columns=to_drop)
        seq_mapping = seq_mapping.dropna(subset=['final_representation'])

    rename_mapping = seq_mapping.set_index('sequence')['final_representation'].to_dict()
    seqtab_df_renamed = seqtab_df.rename(columns=rename_mapping)
    column_order = ["sample"] + list(seqtab_df_renamed.columns)
    seqtab_df_renamed.index.name = 'sample'
    seqtab_df_renamed.reset_index(inplace=True)
    seqtab_df_renamed[column_order].to_csv(args.output, sep="\t", index=False, header=True)

