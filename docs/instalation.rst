

Installation/Compiling
----------------------

Once OpenFOAM is installed and CFDWind3 downloaded, you should follow
these steps: 

1. Make sure you have loaded the OpenFOAM v2.4.1 environment. 
2. Change directory to the CFDWind3 local directory
3. Run ``./Allwclean`` 
4. Run ``./Allwmake`` 
5. Make sure that no error messages appeared and that all libraries and applications are listed as "up to date."

Even though the example case is preconfigured with the forcing files
required to run the tutorial, external libraries of netCFD and Python
programming language are needed if the user wants to modify the
pre-processing and, in particular for post-processing the OpenFOAM
outputs to replicate/analyze/share the results in the format required in
the `GABLS3 benchmarking exercise <http://windbench.net/gabls-3>`__.
netCFD4 system libraries, as well as its Python bindings are available
in most Linux-based OS repositories. In the specific case of Debian-type
systems you can install them as follows:

1. Install netCDF system libs: >
   ``sudo apt install ibnetcdf-dev netcdf-bin``
2. Install (or make sure you have) systems libs: >
   ``sudo apt install build-essential``
3. Install (or make sure you have) pip library management system: >
   ``sudo apt install python-pip``
4. Make sure pip is updated > ``pip install --upgrade pip``
5. Install netCDF4 python bindings > ``pip install netCDF4``

The scripts are known to work with Python v2.7 which in addition
requires the following specific Python libraries (and version): 

- matplotlib (1.5.1) 
- netCDF4 (1.2.9) 
- numexpr (2.6.3) 
- numpy (1.14.0)
- pandas (0.22.0) 
- pip (9.0.1) 
- scipy (1.0.0)

