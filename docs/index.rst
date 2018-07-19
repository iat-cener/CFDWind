

CFDWind3
========

.. figure:: https://windbench.net/system/files/cfdwind3.png
   :alt: CFDWindLogo


Overview
--------

**CFDWind3** is a tool developed at CENER for the simulation of
atmospheric flows whose latests implementations are carried out under
the framework of the `New European Wind Atlas project
(NEWA) <http://www.neweuropeanwindatlas.eu/>`__. The current system is
built under the open-source Computational Fluid Dynamics (CFD) tool-kit
`OpenFOAM version
2.4.0 <https://openfoam.org/download/2-4-0-ubuntu/>`__.

The model is designed to solve the unsteady Reynolds Navier-Stokes
equations (URANS) for incompressible flows in which turbulence closure
is achieved using eddy-viscosity theory and the modified (two-equation)
k-ε closure scheme as proposed by Sogachev et al. (2012).

The solver is based on the Boussinesq approximation by including a
buoyancy term in the momentum equations which, together with the
solution of the energy-transport equation and additional source/sink
terms in the turbulence closure, allows simulating the evolution of the
diurnal cycle. On the other hand, only dry air is considered so neither
humidity transport equation nor heat transfer by radiation or phase
changes is included.

The Boussinesq approximation for incompressible flow is approached in
OpenFOAM through the *buoyantBoussinesqPimpleFoam* solver which
introduces the PISO-SIMPLE algorithm to solve the
pressure-velocity-temperature coupling. This algorithm does not solve
the continuity equation; instead, it solves a pressure Poisson equation
that enforces continuity. Details about the solver algorithm for
OpenFOAM can be found in `OpenFOAM-PISO wiki
description <https://openfoamwiki.net/index.php/BuoyantBoussinesqPisoFoam>`__.

The original OpenFOAM solver is modified following the method proposed
by Sanz-Rodrigo et al. 2017b. That is, Coriolis apparent force and real
large scale forcing are used as model forcing. Forcing comes from the
terms in WRF momentum& energy budget associated to the pressure gradient
and the advection of momentum and temperature. These tendencies are
obtained in the standard output of WRF following the method described in
Lehner (2018a,b).

So far, the tendencies are only height and time-dependent which are
valid for flat-terrain sites. Prior to be introduced in the microscale
model, these terms are stored and averaged horizontally in a 45km area
as described in Chavez-Arroyo et al. 2018. The implementations for
incorporating the tendencies in OpenFOAM are built on top of the `SOWFA
project <https://github.com/NREL/SOWFA>`__ (Churchfield et al. 2014)
developed at the U.S. National Renewable Energy Laboratory (NREL).

The surface conditions comply with the Monin-Obukhov Similarity Theory
(MOST) for neutral-stability condition as proposed by Richards & Hoxey
(1993) and Parente et al. (2011). These conditions are applied as wall
functions to the fields of eddy-viscosity ηt, kinematic thermal
conductivity αt,, turbulence dissipation ε, and Turbulent Kinetic Energy
(TKE). So far, it has been found only small variations in the results
obtained when including stability functions in the turbulence fields,
namely αt and ηt. Therefore, up to this release, only mesoscale
Dirichlet conditions are prescribed for the variations of temperature at
the ground which are introduced through a wall temperature flux field
that uses the algorithm proposed by Basu et al. (2008) to account for
the atmospheric stability in the surface layer. To this end, following
the code released as part of the SOWFA system, the dynamic values of
velocity and temperature scales are computed based on the local flow,
following MOST and Etling (1996) stability functions.



.. toctree::
   :hidden:
   :maxdepth: 2

   code_structure
   instalation
   tutorials
   get_involved
   license
   references
   acnowladgements


