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
source <(curl -L https://raw.githubusercontent.com/vivaxgen/ngs-pipeline/refs/heads/dev/install.sh)

# prepare MicroHaps pipeline environment

echo Cloning vivaxGEN MicroHaps pipeline
git clone --depth 1  https://github.com/vivaxgen/MicroHaps.git ${ENVS_DIR}/MicroHaps
ln -sr ${ENVS_DIR}/MicroHaps/etc/bashrc.d/50-microhaps ${BASHRC_DIR}/

source ${ENVS_DIR}/MicroHaps/etc/inst-scripts/inst-deps.sh

echo "MicroHaps" >> ${ETC_DIR}/installed-repo.txt

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
