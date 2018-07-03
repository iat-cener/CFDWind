# CFDWind3 

## Overview
![CFDWindLogo](https://windbench.net/system/files/cfdwind3.png)

**CFDWind3** is a tool developed at CENER for the simulation of atmospheric flows whose latests implementations are carried out under the framework of the [New European Wind Atlas project (NEWA)](http://www.neweuropeanwindatlas.eu/).
The current system is built under the open-source Computational Fluid Dynamics (CFD) tool-kit [OpenFOAM version 2.4.0](https://openfoam.org/download/2-4-0-ubuntu/).

The model is designed to solve the unsteady Reynolds Navier-Stokes equations (URANS) for incompressible flows in which turbulence closure is achieved using eddy-viscosity theory and the modified (two-equation) k-ε closure scheme as proposed by Sogachev et al. (2012).

The solver is based on the Boussinesq approximation by including a buoyancy term in the momentum equations which, together with the solution of the energy-transport equation and additional source/sink terms in the turbulence closure, allows simulating the evolution of the diurnal cycle. On the other hand, only dry air is considered so neither humidity transport equation nor heat transfer by radiation or phase changes is included.

The Boussinesq approximation for incompressible flow is approached in OpenFOAM through the *buoyantBoussinesqPimpleFoam* solver which introduces the PISO-SIMPLE algorithm to solve the pressure-velocity-temperature coupling. This algorithm does not solve the continuity equation; instead, it solves a pressure Poisson equation that enforces continuity.  Details about the solver algorithm for OpenFOAM can be found in [OpenFOAM-PISO wiki description](https://openfoamwiki.net/index.php/BuoyantBoussinesqPisoFoam).

The original OpenFOAM solver is modified following the method proposed by Sanz-Rodrigo et al. 2017b. That is, Coriolis apparent force and realistic large scale (synoptic and mesoscale) fields are used as model forcing. 
The large scale forcing comes from the terms in WRF momentum & energy budget associated to the pressure gradient and the advection of momentum and temperature. These tendencies are obtained by a previously run WRF simulation in which a small modification to WRF's source code is made in order to include these fields in the standard output of WRF solution as described in Lehner (2018a,b).

So far, the tendencies are only height and time-dependent which are valid for flat-terrain sites. Prior to be introduced in the microscale model, these terms are stored and averaged horizontally in a 45km area as described in Chavez-Arroyo et al. 2018. 
The implementations for incorporating the tendencies in OpenFOAM are built on top of the [SOWFA project](https://github.com/NREL/SOWFA) (Churchfield et al. 2014) developed at the U.S. National Renewable Energy Laboratory (NREL).

The surface conditions comply with the Monin-Obukhov Similarity Theory (MOST) for neutral-stability condition as proposed by Richards & Hoxey (1993) and Parente et al. (2011). 
These conditions are applied as wall functions to the fields of eddy-viscosity ηt, kinematic thermal conductivity αt,, turbulence dissipation ε, and Turbulent Kinetic Energy (TKE).
So far, it has been found only small variations in the results obtained when including stability functions in the turbulence fields, namely αt  and ηt. Therefore, up to this release, only mesoscale Dirichlet conditions are prescribed for the variations of temperature at the ground which are introduced through a wall temperature flux field that uses the algorithm proposed by Basu et al. (2008) to account for the atmospheric stability in the surface layer. To this end, following the code released as part of the SOWFA system, the dynamic values of velocity and temperature scales are computed based on the local flow, following MOST and Etling (1996) stability functions.

## Code structure
CFDWind3 consist of the following list of OpenFOAM libraries and applications that can be easily compiled via the Allwmake file, once the proper environmental variables of OpenFOAM2.4.1 are loaded.
### **Applications**:
1. **Solvers**
   - *ABLRANS_Solver*: Main solver for atmospheric boundary layer turbulent flow. Built on top of ABLSolver from SOWFA project, which in turn is based on the OpenFOAM buoyantBoussinesqPimpleFoam solver. It includes the Coriolis, tendencies forcing with the functionality of including surface roughness, stability and wind speed and direction. As in SOWFA, so far it can only handle structured hexahedral meshes and flat terrain sites.
2. **Utilities**:
   - *setRANSFieldABL*: Departs from the setFieldsABL from SOWFA modified to initialize the flow field in the RANS simulations. It specifies the initial profile of velocity and temperature based either on height-dependent tables, logarithm profile or geostrophic values. Turbulence fields, namely TKE and epsilon are prescribed as in the Ekmann layer.
   - *writeCellCoordinates*: A utility to write in txt format the x,y,z coordinates of the internalField and all patches. Useful for creating the inputs for tools or programs outside the openfoam classes environment.
### **Libraries**
1. **Incompressible/turbulenceModels**: Contains custom boundary conditions based on SOWFA system:
   - *specifiedSurfaceTemperatureq2*: surface temperature flux / heating - Specify a surface cooling/heating rate and the appropriate temperature flux is calculated.
   - *timeVaryingMappedFixed* value with inlet outlet. Useful when WRF results are also applied as inlet boundary conditions.
2. **Incompressible/RAS**:
   - *kEpsilonSogachevPrecursor*. a modified version of OpenFOAM standard kEpsilon model which now include Pr_t and nu_t sensitive to atmospheric stability as proposed by Sogachev et al. 2012 and Koblitz et al. 2013.
   - derivedFvPatchFields/wallFunctions: wall functions to comply with MOST for epsilon, nu_t and epsilon as in Parente et al (2011). Although their current implementation only considers neutral stability, few tests have shown very small differences in the results when extended stability-based MOST functions are included.
### **External code**
1. **WRF**: Contain the modifications required in the WRF source code to obtain the tendencies in the standard output of WRF following the method described in Lehner (2018a,b).

## Installation/Compiling
Once OpenFOAM is installed and CFDWind3 downloaded, you should follow these steps:
1. Make sure you have loaded the OpenFOAM v2.4.1 environment.
2. Change directory to the CFDWind3 local directory
3. Run `./Allwclean`
4. Run `./Allwmake`
5. Make sure that no error messages appeared and that all libraries and applications are listed as "up to date."

Even though the example case is preconfigured with the forcing files required to run the tutorial, external libraries of netCFD and Python programming language are needed if the user wants to modify the pre-processing and, in particular for post-processing the OpenFOAM outputs to replicate/analyze/share the results in the format required in the [GABLS3 benchmarking exercise](http://windbench.net/gabls-3).
netCFD4 system libraries, as well as its Python bindings are available in most Linux-based OS repositories. In the specific case of Debian-type systems you can install them as follows:

1. Install netCDF system libs:
> `sudo apt install libnetcdf-dev netcdf-bin`
2. Install (or make sure you have) systems libs:
> `sudo apt install build-essential`
3. Install (or make sure you have) pip library management system:
> `sudo apt install python-pip`
4. Make sure pip is updated
> `pip install --upgrade pip`
5. Install netCDF4 python bindings
> `pip install netCDF4`


The scripts are known to work with Python v2.7 which in addition requires the following specific Python libraries (and version):
- matplotlib (1.5.1)
- netCDF4 (1.2.9)
- numexpr (2.6.3)
- numpy (1.14.0)
- pandas (0.22.0)
- pip (9.0.1)
- scipy (1.0.0)

## Tutorials
The tutorials are included in the folder “example cases” folder. So far only one case is included (GABLS3), but the list will be growing as new benchmarks from the NEWA experiments and other cases from external contributions are carried out.

### GABLS3: Diurnal cycle driven by mesoscale data

This site is part of the Cabauw Experimental Site for Atmospheric Research (CESAR) in the Netherlands. The example is taken from the benchmark initiative of the third GEWEX Atmospheric Boundary Layer Studies (GABLS3) case revisited for wind energy applications (Sanz-Rodrigo et al. 2017a).
The inter-comparison study aimed to benchmark different meso-microscale coupling strategies with different turbulence-closure fidelity of the underlying microscale models.

The example case aims to contribute to the “open-science” initiative already presented in Sanz-Rodrigo et al. (2017a) so the case is prepared to reproduce the [GABLS3 benchmark](http://windbench.net/gabls-3) described in that exersice whose results were publised in: http://iopscience.iop.org/article/10.1088/1742-6596/854/1/012037
The terrain is flat with laterally periodic boundaries. Due to these conditions and given the RANS approach, the case can be seen as an horizontally-homogeneous case. Thus only a 4x4x200 cells are needed to solve the 1D-type of problem, though mesh configuration can be selected in the *./inputParameters file*. 

### To run the tutorial
In order to run the GABLS3 example case you can follow the next steps:
1. Change directory to the GABLS3 folder
2. As described in the [benchmark](http://windbench.net/gabls-3), download the tendencies input file [*GABLS3_tendencies_d02_YSU_w60_L9000.nc*](https://b2share.eudat.eu/records/22e419b663cb4ffca8107391b6716c1b). So far the preprocessing scripts consider this file to be located in *./inputData/*
If you want to change the directory then edit the new path in the `./runPreprocessing.sh` scrip (tendenciesFile variable)
3. Edit the number of processors to use 
4. Run the `./runPreprocessing.sh` script to set up the mesh and initial fields.
5. Run the `./runCase.sh` script
6. Run the `./runPostprocessing.sh` script to extract the results. 


## License
The files are based on the OpenFOAM software and are either new files or modifications of files in the OpenFOAM source code distribution. 
Therefore, the published under the GNU General Public License version 3.0. This requires anyone distributing the code or a derivative work to make the source available under the same terms, and also provides an express grant of patent rights from contributors to users (GNU Operating System, 2017). The source code shall be downloaded from the OpenFOAM  Foundation (2018) website.
Please see the included OpenFOAM readme file ("README.OpenFOAM") and the GPL licence information ("COPYING"). 

WRF is developed by the National Center for Atmospheric Research (NCAR) and released with a public domain license, which means that the code can be used by any person or entity for any purpose without any fee or charge. The source code shall be downloaded from NCAR’s WRF Model User’s Page (2018) website. 

## References 
- Basu S., et al. (2008). An inconvenient “truth” about using sensible heat flux as a surface boundary condition in models under stably stratified regimes. *Acta Geophysica* (56)1 pp. 88-99. doi: 10.2478/s11600-007-0038-y.
- Chávez-Arroyo, R.A. et al. (2018). Analysis and validation of Weather Research and Forecasting model tendencies for meso-to-microscale modelling of the atmospheric boundary layer. *Submitted for J. Phys.: Conf. Ser.*
- Churchfield M. et al. (2014) Overview of the Simulator fOr Wind Farm Application (SOWFA). Available at: https://nwtc.nrel.gov/system/files/SOWFA_tutorial_05-20-2014.pdf   Source code available at: https://github.com/NREL/SOWFA
- Etling D. (1996). Modelling of Atmospheric Flow Field. Chapter 2 - Modelling the Vertical ABL Structure pp. 56—57. *World Scientific. ISBN: 978-981-02-2509-4*
- GNU Operating System (2018) GNU General Public License v3.0, https://www.gnu.org/licenses/gpl-3.0.en.html, last accessed February 2018
- Hahmann A.N, et al. (2017) Description of the probabilistic wind atlas methodology, NEWA deliverable D3.1, July 2017
- Lehner M (2012a) Extracting terms of the horizontal momentum and thermodynamic equations in the WRF model code. http://home.chpc.utah.edu/~u0653546/WRF_docs/Output_tendency_terms.pdf, last accessed February 2018
- Lehner M (2012b) Observations and large-eddy simulations of the thermally driven cross-basin circulation in a small, closed basin, *Ph.D. thesis, University of Utah*, 2012
- Mann J, et al. (2017) Complex terrain experiments in the New European Wind Atlas. *Phil. Trans. R. Soc. A*, 20160101, doi: 10.1098/rsta.2016.0101
- The OpenFOAM Foundation (2018) OpenFOAM 2.4.0 Released, https://openfoam.org/version/2-4-0/, last accessed February 2018
- Parente A et al. (2011). Improved k–ε model and wall function formulation for the RANS simulation of ABL flows.*J. Wind Eng. Ind. Aerodyn.* 99 267–278 ISSN 01676105
- Petersen EL, Troen I (2012) Wind conditions and resource assessment.* WIREs Energy Environ*, 1: 206–217 doi: 10.1002/wene.4
- Richards P. and Hoxey R (1993). Appropriate boundary conditions for computational wind engineering models using the k–eps turbulence model. *J. Wind Eng. Ind. Aerodyn.* 46-47, 145–153.
- Sanz Rodrigo J, et al. (2016a) Mesoscale to Microscale Wind Farm Modelling and Evaluation. *WIREs Energy Environ*, 6: e214, doi: 10.1002/wene.214 
- Sanz Rodrigo J, et al. (2017a) Results of the GABLS3 diurnal cycle benchmark for wind energy applications. *J. Phys.: Conf. Ser*, 854: 012037, doi :10.1088/1742-6596/854/1/012037
- Sanz Rodrigo, J., Churchfield, M., & Kosovic, B. (2017b). A methodology for the design and testing of atmospheric boundary layer models for wind energy applications. *Wind Energy Science*, 2(1), 35–54. https://doi.org/10.5194/wes-2-35-2017
- Sanz Rodrigo J (2018) Windbench/GABLS3 benchmark. http://windbench.net/gabls-3, last accessed February 2018
- Sogachev A. et al. (2012). Consistent Two-Equation Closure Modelling for Atmospheric Research: Buoyancy and Vegetation Implementations. *Boundary-Layer Meteorol.* 145: 307–327, doi:10.1007/s10546-012-9726-5
- Skamarock WC, et al (2008) A Description of the Advanced Research WRF Version 3. NCAR Tech. Note NCAR/TN-475+STR, 113 pp., doi:10.5065/D68S4MVH
- WRF Model User’s Page (2018) http://www2.mmm.ucar.edu/wrf/users/, last accessed February 2018 


## Acknowledgements 
The NEWA project is supported by a European Commission’s ERA-Net Plus project, number 618122, joining national projects from 9 funding agencies from 8 member states:
- Public Service of Wallonia, Department of Energy and Sustainable Building (Belgium)
- Department of Economy, Science and Innovation Flemish Government (Belgium)
- Danish Energy Authority (Denmark)
- Federal Ministry for the Economic Affairs and Energy, on the basis of the decision by the German Bundestag (Germany)
- Latvijas Zinatnu Akademija (Latvia)
- Fundação para a Ciência e a Tecnologia (Portugal)
- Ministerio de Economía, Industria y Competitividad (Spain)
- The Swedish Energy Agency (Sweden)
- The Scientific and Technological Research Council of Turkey (Turkey)
