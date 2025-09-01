# version.py - ngs-pipeline command line
# [https://github.com/vivaxgen/ngs-pipeline]

__copyright__ = "(c) 2025, Hidayat Trimarsanto <trimarsanto@gmail.com>"
__license__ = "MIT"

# to improve the responsiveness during bash autocomplete, do not import heavy
# modules (such as numpy, pandas, etc) here, but instead import them within the
# functions that require the respective heavy modules

import os
from ngs_pipeline import cout
from ngs_pipeline.cmds import version as ngs_pl_version


init_argparser = ngs_pl_version.init_argparser


def version(args):

    ngs_pl_version.version(args)
    git_hash = ngs_pl_version.get_git_hash(
        os.getenv("MICROHAPS_BASEDIR"), "MicroHaps", verbose=args.verbose
    )
    cout(git_hash)


def main(args):
    version(args)


# EOF
