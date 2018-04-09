/*---------------------------------------------------------------------------*\
This code was developed at the National Renewable Energy Centre of Spain (CENER) 
on Nov, 2017 within the framework of the NEWAFoam package
Roberto Chavez (rchavez@cener.com)
\*---------------------------------------------------------------------------*/

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
    writeInletPatchCoordinates

Description
    Application to write coordinates of inlet patch of an openfoam case
    Created originally by Samuel Chang (Franhoufer IWES)
    Then eddited and updated by  Roberto Chavez (CENER) rchavez@cener.com

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "fvc.H"
#include "IFstream.H"
#include "OFstream.H"

#include "Time.H"
#include "argList.H"
#include "nearWallDist.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


int main(int argc, char *argv[])
{
#   include "setRootCase.H"
#   include "createTime.H"
#   include "createMesh.H"

	const surfaceVectorField& cellFaces = mesh.Cf();

	// Reading ABL dictionary (R.Chavez)
    IOdictionary ABLProperties_
    (
        IOobject
        (
            "ABLProperties",
            runTime.time().constant(),
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    // GENERATING THE INTERNAL FIELDS
    Info<< "Dump face centres of the internal Field\n" << endl;

    OFstream outPutFileInternal
    (
        runTime.path()/"coordinatesInternalField.dat"
    );

    forAll(mesh.cells(), cellI)
    {
        outPutFileInternal	<<
        mesh.cellCentres()[cellI].x() << tab <<
        mesh.cellCentres()[cellI].y() << tab <<
        mesh.cellCentres()[cellI].z() << '\n';
    }

    // GENERATING THE INLET BOUNDARY
    word inletPatchName_(ABLProperties_.lookup("inletPatchName"));
    Info<< "Dump face centres of the " << inletPatchName_<< " patch\n" << endl;
    label inletPatchID = mesh.boundaryMesh().findPatchID(inletPatchName_);
    const fvsPatchVectorField& facesInlet = cellFaces.boundaryField()[inletPatchID];

    OFstream outPutFileInlet
    (
        runTime.path()/"coordinatesInletPatch.dat"
    );

	forAll(facesInlet, faceI)
	{
		outPutFileInlet	<<
		facesInlet[faceI].x() << tab <<
		facesInlet[faceI].y() << tab <<
		facesInlet[faceI].z() << '\n';
	}

    // GENERATING THE GROUND BOUNDARY
	word groundPatchName_(ABLProperties_.lookup("groundPatchName"));
    Info<< "Dump face centres of the " << groundPatchName_<< " patch\n" << endl;
    label groundPatchID = mesh.boundaryMesh().findPatchID(groundPatchName_);
    const fvsPatchVectorField& facesGround = cellFaces.boundaryField()[groundPatchID];

    OFstream outPutFileGround
    (
        runTime.path()/"coordinatesGroundPatch.dat"
    );

	forAll(facesGround, faceI)
	{
		outPutFileGround	<<
		facesGround[faceI].x() << tab <<
		facesGround[faceI].y() << tab <<
		facesGround[faceI].z() << '\n';
	}



    // GENERATING THE TOP BOUNDARY

    Info<< "Dump face centres of the top patch\n" << endl;

    label topPatchID = mesh.boundaryMesh().findPatchID("top");
    const fvsPatchVectorField& facesTop = cellFaces.boundaryField()[topPatchID];

    OFstream outPutFileTop
    (
        runTime.path()/"coordinatesTopPatch.dat"
    );

	forAll(facesTop, faceI)
	{
		outPutFileTop	<<
		facesTop[faceI].x() << tab <<
		facesTop[faceI].y() << tab <<
		facesTop[faceI].z() << '\n';
	}
	//nearWallDist wDistance(mesh);
	//const volScalarField::GeometricBoundaryField& y=wDistance.y();
    //Info<< "y wall dist:" << y[groundPatchID] <<endl;

}

// ************************************************************************* //
