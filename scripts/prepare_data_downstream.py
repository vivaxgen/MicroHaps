import pandas as pd
import argparse 

p = argparse.ArgumentParser()
p.add_argument('--input', type=str, required=True, help='outputHaplotype.tsv')
p.add_argument('--output', type=str, required=True, help='filtered.tsv')
p.add_argument('--min_allele_count', type=int, default=5, help='minimum allele count to keep')
p.add_argument('--min_locus_count', type=int, default=5, help='minimum locus count to keep')
p.add_argument('--min_marker', type=int, default=70, help='minimum markers per sample to keep')
args = p.parse_args()

fin_cols = ["sample_id", "locus", "allele", "count"]

# sample_id, locus, allele, count
def to_long_format(df):
    df = df.melt(id_vars=["sample"], var_name="haplotype", value_name="count")
    df.rename({"sample": "sample_id"}, axis=1, inplace=True)
    df['locus'] = df['haplotype'].str.split(',').str[0]
    df['allele'] = df['haplotype'].str.split(',').str[1]
    return df[fin_cols]

def filter_haplotype_by_count(df, min_allele_count = 5, min_locus_count = 10, min_marker = 111):
    df_long = df.copy()
    
    # 1 allele count filter
    df_long = df_long.query("count >= @min_allele_count") 
    
    # 2 locus count filter
    sample_markers_total = df_long.groupby(["sample_id", "locus"]).sum("count").reset_index().rename({"count":"locus_count"}, axis = 1)
    df_long = df_long.merge(sample_markers_total[["sample_id", "locus", "locus_count"]], on = ["sample_id", "locus"])
    df_long.loc[:, "kept_locus"] = df_long.loc[:, "locus_count"] >= min_locus_count
    df_long = df_long.loc[df_long['kept_locus'], :]
    
    # 3 min marker filter
    passed_sample = df_long.query("count >= @min_allele_count")[fin_cols].groupby(["sample_id", "locus"]).sum("count").reset_index().groupby("sample_id")[['locus']].nunique().reset_index().query('locus >= @min_marker')['sample_id']
    df_long = df_long.query("sample_id.isin(@passed_sample)")
    return df_long[fin_cols]


def filter_haplotype_by_pc_count(df, filter_pc = 0.05, population_aware = True):
    df_long = df.copy()
    sample_markers_total = df_long.groupby(["sample_id", "locus"]).sum("count").reset_index()
    sample_markers_total.loc[:, "min"] = sample_markers_total.loc[:, "count"] * filter_pc
    
    df_long = df_long.merge(sample_markers_total[["sample_id", "locus", "min"]], on = ["sample_id", "locus"])
    df_long.loc[:, "kept_allele"] = df_long.loc[:, "count"] >= df_long.loc[:, "min"]

    if not population_aware:
        df_long.loc[df_long['kept_allele'] == False, "count"] = 0
        return df_long[['sample_id', 'locus', 'allele', 'count']]

    fin_allele = df_long.query('kept_allele == True')[["locus", "allele"]].drop_duplicates()
    fin_allele['final_allele'] = True
    
    df_long = df_long.merge(fin_allele, on=["locus", "allele"])
    df_long.loc[df_long['final_allele'] == False, "count"] = 0
    return df_long[['sample_id', 'locus', 'allele', 'count']]

if __name__ == "__main__":
    df = pd.read_csv(args.input, sep="\t")
    df_long = to_long_format(df)
    df_filt = filter_haplotype_by_count(df_long, min_allele_count=args.min_allele_count, min_locus_count=args.min_locus_count, min_marker=args.min_marker)
    df_filt.to_csv(args.output, sep="\t", index=False)