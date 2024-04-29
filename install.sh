#!/usr/bin/bash

# installation script for vivaxgen MicroHaps pipeline [https://github.com/vivaxgen/MicroHaps]

# optional variable:
# - BASEDIR
# - OMIT

set -eu

# run the base.sh
# Detect the shell from which the script was called
parent=$(ps -o comm $PPID |tail -1)
parent=${parent#-}  # remove the leading dash that login shells have
case "$parent" in
  # shells supported by `micromamba shell init`
  bash|fish|xonsh|zsh)
    shell=$parent
    ;;
  *)
    # use the login shell (basename of $SHELL) as a fallback
    shell=${SHELL##*/}
    ;;
esac

# Parsing arguments
if [ -t 0 ] && [ -z "${BASEDIR:-}" ]; then
  printf "Pipeline base directory? [./vvg-MicroHaps] "
  read BASEDIR
fi

# default value
BASEDIR="${BASEDIR:-./vvg-MicroHaps}"

OMIT="${OMIT:-}"
uMAMBA_ENVNAME='muhaps'
source <(curl -L https://raw.githubusercontent.com/vivaxgen/ngs-pipeline/main/install.sh)

echo "Installing latest bedtools"
micromamba -y install bedtools -c conda-forge -c bioconda -c defaults

echo "Installing datamash"
micromamba -y install datamash -c conda-forge -c bioconda -c defaults

echo "Installing trimmomatic"
micromamba -y install trimmomatic -c conda-forge -c bioconda -c defaults

echo "Installing mosdepth"
micromamba -y install mosdepth -c conda-forge -c bioconda -c defaults

echo "Installing samclip"
micromamba -y install samclip -c conda-forge -c bioconda -c defaults

echo "Installing muscle version v3.8.1551"
micromamba -y install "muscle=3.8.1551" -c conda-forge -c bioconda -c defaults

echo "Installing required R packages"
micromamba -y install r-ggplot2 r-BiocManager r-RCurl r-argparse r-data.table r-seqinr r-doMC -c conda-forge -c bioconda -c defaults

# Install additional R packages
R --no-save << EOF
BiocManager::install("GenomeInfoDb")
BiocManager::install("GenomicRanges")
BiocManager::install("Biostrings")
BiocManager::install("Rsamtools")
BiocManager::install("SummarizedExperiment")
BiocManager::install("GenomicAlignments")
BiocManager::install("ShortRead")
BiocManager::install("dada2")
BiocManager::install("limma")
EOF

echo "installing required Python modules"

# to use latest of all python-related stuff, uncomment below and remove the conda parts
pip3 install biopython
pip3 install cutadapt
pip3 install tqdm

# prepare MicroHaps pipeline environment

echo Cloning vivaxGEN MicroHaps pipeline
git clone https://github.com/vivaxgen/MicroHaps.git ${ENVS_DIR}/MicroHaps
ln -sr ${ENVS_DIR}/MicroHaps/etc/bashrc.d/50-microhaps ${BASHRC_DIR}/

echo "Reloading profiles"
reload_vvg_profiles

echo "Initialize enviroment"
ngs-pl initialize --target wgs

echo ""
echo "vivaxGEN MicroHaps pipeline has been successfully installed."
echo "Please run the activation file with the following command:"
echo ""
echo "    `realpath ${BINDIR}/activate`"
echo ""
echo "or source the activation file (eg. inside a script):"
echo ""
echo "    source `realpath ${BINDIR}/activate`"
echo ""

# EOF
