# MicroHaps-specific activation script

## find our own path
_script="$(readlink -f ${BASH_SOURCE[0]})"

## Delete last component from $_script ##
_mydir="$(dirname $_script)"

export MICROHAPS_BASEDIR=${_mydir}
export NGSENV_BASEDIR=${_mydir}

PATH=${NGSENV_BASEDIR}/bin:${PATH}
PYTHONPATH=${NGSENV_BASEDIR}:${PYTHONPATH}

## variable for vivaxGEN NGS-Pipeline
export NGS_PIPELINE_CMD_MODS=muhaps_pipeline.cmds:microhaps_pipeline.cmds:${NGS_PIPELINE_CMD_MODS}
export NGS_PIPELINE_FORCE=1
export NGS_PIPELINE_NO_CONFIG_CASCADE=1

## uncomment below if you want to  avoid running "python setup.py install"
PATH=${MICROHAPS_BASEDIR}/scripts:${PATH}

## set prompt
PS1=$'(\xc2\xb5haps) [\\u@\\h \\W]$ '

export HISTFILE=${HOME}/.bash_history.muhaps
history -c; history -r

