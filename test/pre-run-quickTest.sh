#!/bin/bash
##********************************************************************************
##  Author:  Paweł Gancarski (CENER)  pgancarski@cener.com
##  Created  11/Jun/2018
###############################################################################


###
### To run the script from the host OS use:
### docker run -it test /bin/bash -c "source ~/OpenFOAM/OpenFOAM-2.4.0/etc/bashrc && ./test/pre-run-quickTest.sh"
###


cd exampleCases/GABLS3/

# we only want to run few iterations, so we change the "theEnd" number of iterations to 600
sed -i "s 100770 600 g"  inputParameters

./runPreprocessing.sh 


# run the case and test the time

# Expected results for Intel® Core™ i5-7300U CPU @ 2.60GHz × 4  (2 phisical cores)
# real	1m17.789s
# user	4m14.593s
# sys	0m42.955s

time ./runCase.sh


# check some prelimenary results

echo "SourceTHistory, expected: -1.13099866569e-05"
tail -c 20 ./postProcessing/SourceHistory/0/SourceTHistory

echo "SourceUXHistory, expected: -0.000138757210702"
tail -c 20 ./postProcessing/SourceHistory/0/SourceUXHistory

echo "SourceUYHistory, expected: -0.000186824395247"
tail -c 20 ./postProcessing/SourceHistory/0/SourceUYHistory

echo "SourceUZHistory, expected: 0 0 0 0 0 0 0 0 0 0"
tail -c 20 ./postProcessing/SourceHistory/0/SourceUZHistory 


