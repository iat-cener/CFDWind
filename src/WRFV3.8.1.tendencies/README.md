# WRFV3.8.1.tendencies 
Javier Sanz Rodrigo, June 2017

## Introduction

This directory provides WRF source files to add tendencies to the standard output (wrfout). Original files come from Lehner (2012), upgraded to WRFV3.8.1 and added potential temperature tendencies. 

See Registry.EM_COMMON.tendencies for definition of output variables. For instance, ru_tend_adv is the u-advection tendency, directly computed by WRF using a mass-coupled system of equations. To decouple the tendencies you need to divide each term with the corresponding dry air mass (mut at mass-points, muu at u-points or muv at v-points). This is done by post-processing in the WRFtend.py script. This script will produce de-staggered time-height fields at a given site with interpolated or spatial-averaged quantities (Sanz Rodrigo et al., 2017). 

Additionaly, wrf_timeseries.F.addedWW, adds the W-component to the time-series output when using tslist.

## References

Lehner M (2012) Extracting  terms  of  the  horizontal  momentum  and  thermodynamic equations in the WRF model code. http://home.chpc.utah.edu/~u0653546/WRF_docs/Output_tendency_terms.pdf

Lehner M (2012) Observations and large-eddy simulations of the thermally driven cross-basin circulation in a small, closed basin, Ph.D. thesis, University of Utah, 2012

Sanz Rodrigo J, Churchfield M, Kosović B (2017) A methodology for the design and testing of atmospheric boundary layer models for wind energy applications. Wind Energ. Sci. 2: 1-20, doi:10.5194/wes-2-1-2017
