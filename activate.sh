# MicroHaps-specific activation script

## find our own path
_script="$(readlink -f ${BASH_SOURCE[0]})"

## Delete last component from $_script ##
_mydir="$(dirname $_script)"

export MICROHAPS_BASEDIR=${_mydir}

## set everything else here
NGS_PIPELINE_BASH_SOURCE=${VVG_BASEDIR}/env/ngs-pipeline/bin/activate.sh

export NGSENV_BASEDIR=${_mydir}
. ${NGS_PIPELINE_BASH_SOURCE}
PATH=${NGSENV_BASEDIR}/bin:${PATH}
PYTHONPATH=${NGSENV_BASEDIR}:${PYTHONPATH}
export NGS_PIPELINE_CMD_MODS=muhaps_pipeline.cmds:${NGS_PIPELINE_CMD_MODS}

## uncomment below if you want to  avoid running "python setup.py install"
PATH=${MICROHAPS_BASEDIR}/scripts:${PATH}


