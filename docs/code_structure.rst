

Code structure
--------------

CFDWind3 consist of the following list of OpenFOAM libraries and
applications that can be easily compiled via the Allwmake file, once the
proper environmental variables of OpenFOAM2.4.1 are loaded.

Applications:
^^^^^^^^^^^^^

    1. **Solvers** 
        - *ABLRANS\_Solver*: Main solver for atmospheric boundary layer turbulent flow. Built on top of ABLSolver from SOWFA project, which in turn is based on the OpenFOAM buoyantBoussinesqPimpleFoam solver. It includes the Coriolis, tendencies forcing with the functionality of including surface roughness, stability and wind speed and direction. As in SOWFA, so far it can only handle structured hexahedral meshes and flat terrain sites.         

    2. **Utilities**: 
        - *setRANSFieldABL*: Departs from the setFieldsABL from SOWFA modified to initialize the flow field in the RANS simulations. It specifies the initial profile of velocity and temperature based either on height-dependent tables, logarithm profile or geostrophic values. Turbulence fields, namely TKE and epsilon are prescribed as in the Ekmann layer. 
        - *writeCellCoordinates*: A utility to write in txt format the x,y,z coordinates of the internalField and all patches. Useful for creating the inputs for tools or programs outside the openfoam classes environment.

Libraries
^^^^^^^^^

    1. **Incompressible/turbulenceModels** - contains custom boundary conditions based on SOWFA system: 
        - *specifiedSurfaceTemperatureq2*: surface temperature flux / heating - Specify a surface cooling/heating rate and the appropriate temperature flux is calculated. 
        - *timeVaryingMappedFixed* value with inlet outlet.Useful when WRF results are also applied as inlet boundary conditions.
    2. **Incompressible/RAS**: 
        - *kEpsilonSogachevPrecursor*. a modified version of OpenFOAM standard kEpsilon model which now include Pr\_t andnu\_t sensitive to atmospheric stability as proposed by Sogachev et al. 2012 and Koblitz et al. 2013. 
        - derivedFvPatchFields/wallFunctions: wall functions to comply with MOST for epsilon, nu\_t and epsilon as in Parente et al (2011). Although their current implementation only considers neutral stability, few tests have shown very smalldifferences in the results when extended stability-based MOST functions are included. 

External code
^^^^^^^^^^^^^

    1. **WRF**: Contain the modifications required in the WRF source code to obtain the tendencies in the standard output of WRF following the method described in Lehner (2018a,b).



