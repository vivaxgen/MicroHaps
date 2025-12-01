import pysam 
import cstag
import argparse 
import pandas as pd

def new_cigar_md_seq(alg):
    alg_cigar = []
    current_op = None
    current_length = 0
    seqs = []
    for alg_pair in alg.get_aligned_pairs(with_seq=True, with_cigar=True):
        query_pos, ref_pos, base, cigar_op = alg_pair
        if ref_pos is None:
            continue

        seqs.append(base)
        if cigar_op != current_op:
            if current_op is not None:
                alg_cigar.append((current_op, current_length))
            current_op = cigar_op
            current_length = 1
        else:
            current_length += 1

    if current_op is not None and current_length > 0:
        alg_cigar.append((current_op, current_length))
    
    CIGAR_OPERATIONS = "MIDNSHP=X"
    cigar_string = ""
    for operation_code, length in alg_cigar:
        operation_char = CIGAR_OPERATIONS[operation_code]
        cigar_string += str(length) + operation_char
    return cigar_string, alg.get_tag("MD"), "".join(seqs)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract overlapping regions from bam to seqtab")
    parser.add_argument("--input", help="Input BAM file")
    parser.add_argument("--output", help="Out seqtab")
    args = parser.parse_args()

    bamfile = pysam.AlignmentFile(args.input, "rb")

    results = []
    for alg in bamfile:
        if alg.is_unmapped:
            continue
        else:
            cigar, md, seq = new_cigar_md_seq(alg)
            cs_str = cstag.call(cigar, md, alg.query, long=False)

            result = pd.DataFrame({
                "qseq": [alg.query_sequence],
                "qname": [alg.query_name],
                "locus": [alg.reference_name],
                "allele": [cs_str],
            })
            results.append(result)
    
    final_df = pd.concat(results, ignore_index=True)
    final_df.loc[:, "allele_str"] = final_df.loc[: , "locus"] + "," + final_df.loc[:, "allele"]
    final_df[["qname", "qseq", "allele_str"]].to_csv(args.output, sep="\t", index=False)

    bamfile.close()