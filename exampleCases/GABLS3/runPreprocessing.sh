#!/bin/bash
##********************************************************************************
## Script to perform the preprocessing steps for running NEWAFoam
##
##  Author:  Roberto Chavez (CENER)  rchavez@cener.com
##  Created  13/Dic/2017
###############################################################################
# REQUIRED INPUTS!!
###############################################################################

# path of openfoam case (normally the "current directory"
export OFPATH=${PWD}   
#fileName of the mesoscale tendencies (modify if necessary to where your input is)
export tendenciesFile=${OFPATH}/inputData/GABLS3_tendencies_d02_YSU_w60_L9000.nc  
# Initial date for the simulation
export datefrom="2006-07-01 08:00"    # the dates format HAS TO BE: "yyyy-mm-dd HH:MM"
# End date for the simulation
export dateto="2006-07-02 12:00"      # the dates format HAS TO BE: "yyyy-mm-dd HH:MM"
# Plot the WRF tendencies and CFD inputs?
export plotncData=0

###############################################################################
# SCRIPT
###############################################################################
export PYTHON_COMMON_LIBRARIES=${OFPATH}/pyFiles/lib/

rm -rf 0
cp -r 0.orig 0

echo "Initialize GABLS3 settings"
python ./pyFiles/createForcing.py > log.forcing
if [ $? -ne 0 ]; then  exit 10; fi

echo "Creating the mesh (blockMesh)"
blockMesh > log.blockMesh

echo "set the initial fields"
setRANSfieldsABL > log.setRANSfieldsABL


