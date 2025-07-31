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