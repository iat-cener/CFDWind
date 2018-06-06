#!/bin/bash
##********************************************************************************
##  Author:  Roberto Chavez (CENER)  rchavez@cener.com
##  Created  13/Jun/2017
###############################################################################

nCores=4   # example of 4 cores
solver=ABLRANS_Solver

echo "decompose"
decomposePar -force > log.decompose

echo "running the model"
mpirun -np ${nCores} ${solver} -parallel > log.modelRun

