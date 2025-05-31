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
    p = arg_parser("generate plots from tabular QC depth table")
    p.add_argument(
        "--additional-title", default="", help="additional title for the plot"
    )
    p.add_argument("--index-column", default="Amplicon_name", help="index column name")
    p.add_argument(
        "--nolog",
        action="store_true",
        default=False,
        help="do not use log scale for heatmap",
    )
    p.add_argument(
        "--value-is-ratio", action="store_true", default=False, help="values are ratios"
    )
    p.add_argument("--outheatmap", help="output heatmap plot for depth")
    p.add_argument("--outmarkers")
    p.add_argument("--outsamples")
    p.add_argument("infile")
    return p


def tab_to_plots(args):

    import pandas as pd
    import numpy as np
    import seaborn as sns
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm

    depth = pd.read_table(args.infile)

    try:
        column_idx = depth.columns.get_loc(args.index_column)
    except KeyError:
        cexit(f"Error: Column '{args.index_column}' not found in the input file.")

    depth_map = depth.iloc[:, (column_idx + 1) :].set_index(depth[args.index_column])

    if args.value_is_ratio:
        args.nolog = True
        vmin, vmax = 0, 1
    else:
        vmin, vmax = depth_map.min().min(), depth_map.max().max()

    if args.outheatmap:

        # plot heatmap for depth
        plt.figure(figsize=(len(depth_map) * 0.2 + 5, len(depth_map.columns) * 0.2 + 5))
        sns.heatmap(
            depth_map,
            cmap="YlGnBu",
            xticklabels=1,
            yticklabels=1,
            vmin=vmin,
            vmax=vmax,
            cbar_kws={"label": "Depth"},
            norm=LogNorm() if not args.nolog else None,
        )
        plt.title(f"{args.additional_title} Depth Heatmap")
        plt.tight_layout()
        plt.savefig(args.outheatmap, dpi=300)
        plt.close()

    if args.outmarkers:

        if not args.value_is_ratio:
            total_per_sample = depth_map.sum(axis=0)
            marker_ratio = depth_map.div(total_per_sample, axis=1)
            marker_ratio = marker_ratio.fillna(0)
        else:
            marker_ratio = depth_map.copy()
        marker_ratio["__median__"] = marker_ratio.median(axis=1)
        marker_ratio = marker_ratio.sort_values(by="__median__")
        del marker_ratio["__median__"]
        marker_ratio_df = marker_ratio.stack().reset_index()
        marker_ratio_df.rename(columns={0: "ratio"}, inplace=True)
        marker_ratio_df["log_ratio"] = np.log2(marker_ratio_df["ratio"]).fillna(0)

        plt.figure(figsize=(20, 80))

        sns.catplot(
            data=marker_ratio_df,
            y="Amplicon_name",
            x="log_ratio" if not args.nolog else "ratio",
            hue="level_1",
            palette="YlGnBu",
            legend=False,
            s=3,
        )
        plt.yticks(size=3)
        plt.ylabel("Amplicon Name")
        plt.xlabel("Log2 Ratio" if not args.nolog else "Ratio")
        plt.title(f"{args.additional_title} Marker Ratio")
        plt.tight_layout()
        plt.savefig(args.outmarkers, dpi=600)
        plt.close()


def main(args):
    tab_to_plots(args)


# EOF
