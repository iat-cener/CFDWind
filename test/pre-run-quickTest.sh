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

## If exist, remove the previous results
rm processor0 -R
rm processor1 -R
rm processor2 -R
rm processor3 -R

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

echo "epsilon" 
echo "   0.00190490069538"
tail -c 20 ./postProcessing/probes/0/epsilon
echo ""

echo "k" 
echo " 0.724389754221"
tail -c 16 ./postProcessing/probes/0/k
echo ""

echo "nut" 
echo "  8.2640609659"
tail -c 15 ./postProcessing/probes/0/nut
echo ""

echo "p" 
echo "  -1501.17742639"
tail -c 17 ./postProcessing/probes/0/p
echo ""

echo "U" 
echo "   (-5.96517598304 -0.303756960366 -1.59169069063e-07)"
tail -c 55 ./postProcessing/probes/0/U
echo ""

echo "T " 
echo " 294.353402946"
tail -c 15 ./postProcessing/probes/0/T
echo ""

echo "SourceU" 
echo "   (-0.00036463370292 -0.000882508253524 0)"
tail -c 44 ./postProcessing/probes/0/SourceU
echo ""

echo "SourceT"
echo "   -7.8929842225e-05"
tail -c 21 ./postProcessing/probes/0/SourceT




