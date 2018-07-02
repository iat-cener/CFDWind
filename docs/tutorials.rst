
Tutorials
---------


The tutorials are included in the folder “example cases” folder. So far
only one case is included (GABLS3), but the list will be growing as new
benchmarks from the NEWA experiments and other cases from external
contributions are carried out.

GABLS3: Diurnal cycle driven by mesoscale data
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This site is part of the Cabauw Experimental Site for Atmospheric
Research (CESAR) in the Netherlands. The example is taken from the
benchmark initiative of the third GEWEX Atmospheric Boundary Layer
Studies (GABLS3) case revisited for wind energy applications
(Sanz-Rodrigo et al. 2017a). The inter-comparison study aimed to
benchmark different meso-microscale coupling strategies with different
turbulence-closure fidelity of the underlying microscale models.

The example case aims to contribute to the “open-science” initiative
already presented in Sanz-Rodrigo et al. (2017a) so the case is prepared
to reproduce the `GABLS3 benchmark <http://windbench.net/gabls-3>`__
described in that exersice whose results were publised in:
http://iopscience.iop.org/article/10.1088/1742-6596/854/1/012037 The
terrain is flat with laterally periodic boundaries. Due to these
conditions and given the RANS approach, the case can be seen as an
horizontally-homogeneous case. Thus only a 4x4x200 cells are needed to
solve the 1D-type of problem, though mesh configuration can be selected
in the *./inputParameters file*.

To run the tutorial
"""""""""""""""""""

In order to run the GABLS3 example case you can follow the next steps:

1. Change directory to the GABLS3 folder.
2. As described in the `benchmark <http://windbench.net/gabls-3>`__, download the tendencies input file `*GABLS3\_tendencies\_d02\_YSU\_w60\_L9000.nc* <https://b2share.eudat.eu/records/22e419b663cb4ffca8107391b6716c1b>`__. So far the preprocessing scripts consider this file to be located in *./inputData/* If you want to change the directory then edit the new path in the ``./runPreprocessing.sh`` scrip (tendenciesFile variable) 
3. Edit the number of processors to use 
4. Run the ``./runPreprocessing.sh`` script to set up the mesh and initial fields.
5. Run the ``./runCase.sh`` script 
6. Run the ``./runPosprocessing.sh`` script to extract the results.



