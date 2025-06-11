#!/usr/bin/env ngs-pl

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

    p = arg_parser("Generate primer balancer table")
    p.add_argument("--targetcolumn", default="Amplicon_name")
    p.add_argument("--filtercolumn", default="")
    p.add_argument("--outfile", default="balanced-parameter.txt")
    p.add_argument("--reportfile", default="balanced-report.txt")
    p.add_argument("infile")

    return p


def balance_primers(df):
    """df should have been indexed by targets"""

    from scipy import stats
    import pandas as pd

    df_normed = df.div(df.sum(axis=0), axis=1)
    median_fraction = df_normed.median(axis=1)
    scale_100 = median_fraction / median_fraction.sum()
    pool_weighting = scale_100**-0.561
    scale_to_min = pool_weighting / pool_weighting.min()
    clip_to_max = scale_to_min.clip(upper=10)

    Sum = clip_to_max.sum()
    # cut 0.25 from head & tail, ie use 50% data
    IQM = stats.trim_mean(clip_to_max, 0.25)
    centre = IQM / Sum * 250000
    dilution_factor = centre / 40

    table = pd.DataFrame(
        {
            "Median_Frac": median_fraction,
            "Scale100": scale_100,
            "PoolWeighting": pool_weighting,
            "ScaleToMin": scale_to_min,
            "ClipToMax10": clip_to_max,
        },
        index=df.index,
    )

    results = {"Sum": Sum, "IQM": IQM, "[centre]": centre, "DF": dilution_factor}

    return (table, results)


def generate_balanced_primer_table(args):

    import pandas as pd

    in_df = pd.read_table(args.infile, sep="\t")

    # filtering
    # regex = "^%s$" % args.targetcolumn
    # if args.filtercolumn:
    #    regex = regex + "|" + args.filtercolumn

    # filtered_df = in_df.filter(regex=regex).set_index(args.targetcolumn)
    filtered_df = in_df.iloc[:, 3:].set_index(args.targetcolumn)
    print("Columns:", filtered_df.columns)

    table, results = balance_primers(filtered_df)

    table.to_csv(args.outfile, sep="\t", index=True, header=True)
    with open(args.reportfile, "w") as fout:
        fout.write(
            "Sum      : %d\nIQM     : %f\n[centre] : %f\nDilution factor: %f\n"
            % (results["Sum"], results["IQM"], results["[centre]"], results["DF"])
        )


def main(args):
    generate_balanced_primer_table(args)


# EOF
