#!/bin/bash
##********************************************************************************
# REQUIRED INPUTS!!
###############################################################################

export OFCASE=${PWD##*/}
# path of openfoam case (normally the "current directory"
export OFPATH=${PWD}
# path of extra libraries for python scripts
export PYTHON_COMMON_LIBRARIES=${OFPATH}/pyFiles/lib/

###############################################################################
# SCRIPT
###############################################################################
echo "postprocessing in pyhton case: "${OFCASE}
python ./pyFiles/OFsets2nc.py > log.of2nc

