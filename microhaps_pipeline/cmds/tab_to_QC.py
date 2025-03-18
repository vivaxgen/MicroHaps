# tab_to_QC.py - Microhaps command
# [https://github.com/vivaxgen/Microhaps]

__author__ = "Hidayat Trimarsanto"
__copyright__ = "(C) 2025, Hidayat Trimarsanto"
__email__ = "trimarsanto@gmail.com,hidayat.trimarsanto@menzies.edu.au"
__license__ = "MIT"

# to improve the responsiveness during bash autocomplete, do not import heavy
# modules (such as numpy, pandas, etc) here, but instead import them within the
# functions that require the respective heavy modules

import sys
import os
import pathlib
from ngs_pipeline import cerr, cexit, arg_parser


def init_argparser():
    p = arg_parser("run joint microhaplotype caller")
    p.add_argument("-d", "--mindepth", type=int, default=10)
    p.add_argument("--outdir", default="")
    p.add_argument("infile")
    return p


def plot_missingness(values, ylabel, xlabel, outfile):

    import seaborn as sns
    import matplotlib.pyplot as plt

    ax = sns.scatterplot(
        x=range(len(values)), y=values.sort_values(), s=3, linewidth=0, alpha=1.0
    )
    ax.set_ylabel(ylabel)
    ax.set_xlabel(xlabel)
    plt.savefig(outfile, dpi=600)
    plt.close()


def tab_to_QC(args):

    import pathlib
    import pandas as pd

    outdir = pathlib.Path(args.outdir)
    df = pd.read_table(args.infile)

    tab_df = pd.DataFrame(df["sample"])

    for col in df.columns[1:]:

        marker = col.split(",")[0]

        if not marker in tab_df.columns:
            tab_df[marker] = 0

        tab_df[marker] += df[col]

    missingness = tab_df.iloc[:, 1:] < args.mindepth

    sample_missingness = pd.DataFrame(df["sample"])
    sample_missingness["MISS"] = missingness.sum(axis=1)
    sample_missingness["F_MISS"] = sample_missingness["MISS"] / (
        len(tab_df.columns) - 1
    )

    marker_missingness = (
        missingness.sum(axis=0)
        .reset_index()
        .rename(columns={"index": "MARKER", 0: "MISS"})
    )
    marker_missingness["F_MISS"] = marker_missingness["MISS"] / len(df["sample"])

    sample_missingness.to_csv(outdir / "sample_missingness.tsv", sep="\t", index=False)
    marker_missingness.to_csv(outdir / "marker_missingness.tsv", sep="\t", index=False)

    plot_missingness(
        sample_missingness["F_MISS"],
        "Microhaps marker missingness",
        "Sample Index",
        outdir / "sample_missingness.png",
    )

    plot_missingness(
        marker_missingness["F_MISS"],
        "Microhaps sample missingness",
        "Marker Index",
        outdir / "marker_missingness.png",
    )
    markers = tab_df.columns[1:]
    chrom = [m.split(":")[0] for m in markers]
    start = [m.split(":")[1].split("-")[0] for m in markers]
    end = [m.split(":")[1].split("-")[1] for m in markers]
    amplicon = [m.split(":")[2] for m in markers]
    new_df = pd.DataFrame(
        {
            "Chr": chrom,
            "Start": start,
            "End": end,
            "Amplicon_name": amplicon,
        }
    )
    new_df = pd.concat([new_df, pd.DataFrame(tab_df.T.values[1:,:], columns=tab_df.T.values[0,:])], axis = 1)
    new_df["temp"] = new_df["Chr"].str.split("_").str[1]
    new_df = new_df.sort_values(by=["temp", "Start"])
    new_df = new_df.drop(columns=["temp"])
    new_df.to_csv(outdir / "depths.tsv", sep="\t", index=False)


def main(args):
    tab_to_QC(args)


# EOF
