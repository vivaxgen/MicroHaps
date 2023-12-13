# MicroHaps-specific activation script

## find our own path
_script="$(readlink -f ${BASH_SOURCE[0]})"

## Delete last component from $_script ##
_mydir="$(dirname $_script)"

export MICROHAPS_BASEDIR=${_mydir}

## set everything else here

## uncomment below if you want to  avoid running "python setup.py install"
PATH=${MICROHAPS_BASEDIR}/scripts:${PATH}


