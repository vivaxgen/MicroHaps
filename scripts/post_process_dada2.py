
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


if __name__ == "__main__":
    import argparse
    from ngs_pipeline import cerr

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
    # each sequence is only mapped to a single reference sequence
    assert paf_df['qname'].is_unique, "Each sequence must be mapped to a single reference sequence"
    
    # non sequence not mapped to any reference sequence
    assert paf_df.query("tname == '*'").shape[0] == 0, "Some sequences are not mapped to any reference sequence"

    paf_df.loc[:, "final_representation"] = paf_df.loc[:, "tname"] + "," + paf_df.loc[:, "cs"]

    assert paf_df['final_representation'].is_unique, "Final representation must be unique"

    seq_mapping = pd.DataFrame.from_dict(fasta_seqs, orient = 'index', columns = ['sequence'])
    seq_mapping = seq_mapping.merge(paf_df[['qname', 'final_representation']], left_index=True, right_on='qname', how='left')
    rename_mapping = seq_mapping.set_index('sequence')['final_representation'].to_dict()
    seqtab_df_renamed = seqtab_df.rename(columns=rename_mapping)
    column_order = ["sample"] + list(seqtab_df_renamed.columns)
    seqtab_df_renamed.index.name = 'sample'
    seqtab_df_renamed.reset_index(inplace=True)
    seqtab_df_renamed[column_order].to_csv(args.output, sep="\t", index=False, header=True)

