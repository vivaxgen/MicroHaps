__author__ = 'Hidayat Trimarsanto'

# set variables and configurations
import os
import pathlib
from ngs_pipeline import fileutils


microhaps_basedir = os.environ['MICROHAPS_BASEDIR']

# resetting parameters from ngs-pipeline
deduplicate = False
keep_paired_bam = True
keep_final_bam = True


def get_abspath(p, prefix=microhaps_basedir):
    if p is None:
        return p

    filepath = p.as_posix() if isinstance(p, pathlib.Path) else p
    if (
        filepath.startswith("/")
        or filepath.startswith("./")
        or filepath.startswith("../")
    ):
        return p

    prefix = pathlib.Path(prefix) if isinstance(prefix, str) else prefix
    return (prefix / p).absolute().as_posix()


def is_nextseq_or_novaseq():
    return config['instrument'].lower().startswith('nextseq') or config['instrument'].lower().startswith('novaseq')


targetregion_file = get_abspath(fn if (fn := config.get('targetregion_file')) else None, microhaps_basedir)
insertseq_file = get_abspath(config['insertseq'], microhaps_basedir)
primer_fw_file = get_abspath(config['primer_fw'], microhaps_basedir)
primer_rev_file = get_abspath(config['primer_rev'], microhaps_basedir)

# define all output files 

outdir = config['outdir']
#in_dir = config['indir']
indir=''

# set read mode
singleton = config.get('singleton', None)
paired = config.get('paired', None)
if singleton:
    read_mode = fileutils.ReadMode.SINGLETON
elif paired:
    read_mode = fileutils.ReadMode.PAIRED_END
else:
    read_mode = None


read_files = fileutils.ReadFileDict(config['infiles'],
                                    underscore=config['underscore'],
                                    mode=read_mode,
                                    skip_list=config.get('skip_list', []),
)
IDs = read_files.keys()


# EOF
