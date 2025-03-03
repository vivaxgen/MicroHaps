# run_discovery_variant_caller.py - ngs-pipeline command
# [https://github.com/vivaxgen/ngs-pipeline]

__copyright__ = "(C) 2024, Hidayat Trimarsanto <trimarsanto@gmail.com>"
__license__ = "MIT"

# to improve the responsiveness during bash autocomplete, do not import heavy
# modules (such as numpy, pandas, etc) here, but instead import them within the
# functions that require the respective heavy modules


from ngs_pipeline import cexit
from ngs_pipeline.cmds import run_multistep_variant_caller


def init_argparser():
    return run_multistep_variant_caller.init_argparser()


def run_discovery_variant_caller(args):
    if not args.panel:
        cexit("Please provide a panel to use using --panel argument")
    args.paired = True if not args.single else False
    run_multistep_variant_caller.run_multistep_variant_caller(args)


def main(args):
    run_discovery_variant_caller(args)


# EOF
