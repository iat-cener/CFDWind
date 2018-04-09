/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    NEWAFoamv0

Description
    Transient solver for buoyant, turbulent flow of incompressible 
    atmospheric flows. It is strongly based on the SOWFA ABLSolver:
    https://github.com/NREL/SOWFA
    (which in turn is based on the buoyantBoussinesqPimpleFoam solver) 
    sligthly modified for the specific RAS turbulence model of Sogachev et al. 2012 

    Uses the Boussinesq approximation:
    \f[
        rho_{k} = 1 - beta(T - T_{ref})
    \f]

    where:
        \f$ rho_{k} \f$ = the effective (driving) kinematic density
        beta = thermal expansion coefficient [1/K]
        T = temperature [K]
        \f$ T_{ref} \f$ = reference temperature [K]

    Valid when:
    \f[
        \frac{beta(T - T_{ref})}{rho_{ref}} << 1
    \f]

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
//#include "RASModel.H"
#include "turbulenceModel.H"

#include "fvIOoptionList.H"
#include "pimpleControl.H"
#include "fixedFluxPressureFvPatchScalarField.H"

//#include "IFstream.H"
//#include "OFstream.H"
//#include "wallDist.H"
//#include "interpolateXY.H"
//#include "interpolateSplineXY.H"
#include "interpolate2D.H"
#include "windRoseToCartesian.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
	#include "createPostProcessingDir.H"
	#include "findVerticalCellLevels.H"
    #include "readGravitationalAcceleration.H"
    #include "createFields.H"
	#include "createSourceTerms.H"

    #include "createFvOptions.H"
    #include "initContinuityErrs.H"
    #include "readTimeControls.H"
    #include "CourantNo.H"
    #include "setInitialDeltaT.H"

    pimpleControl pimple(mesh);

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;
        Info << "Time Step = " << runTime.timeIndex() << endl;

        #include "readTimeControls.H"
        #include "CourantNo.H"
        #include "setDeltaT.H"

        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
            Info << "   Predictor..." << endl;
            #include "UEqn.H"
            #include "turbulenceCorrect.H"
            #include "TEqn.H"

            // --- Pressure corrector loop
            int corr = 0;
            while (pimple.correct())
            {
            	Info << "   Corrector Step " << corr << "..." << endl;
                #include "pEqn.H"
                // In SOWFA code, Lapointe-Theriault and Churchfield et al. 2014
                // found that the temperature eq. should be called also in the 
                // corrector loop to better couple the system. However in 
                // RANS application this step was found to introduce numerical 
                // instabilities and since the temp-momentum couple was weaker 
                // it was decided to leave it  as in the original buoyantBoussinesqPimpleFoam solver
                //#include "turbulenceCorrect.H"  // SOWFA ABLSolver 
                //#include "TEqn.H"  // SOWFA ABLSolver 
                corr++;
            }


            // --- Update the source terms
            #include "correctSourceTerms.H"   


            //if (pimple.turbCorr())
            //{
            //    turbulence->correct();
            //}

            // --- Update temperature flux conditions
            //qwall.correctBoundaryConditions();

        }

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
