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

\*---------------------------------------------------------------------------*/

#include "nutAtmWallFunctionFvPatchScalarField.H"
#include "incompressible/turbulenceModel/turbulenceModel.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"

#define z0min 0.0002  //rchavez (2016)

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

tmp<scalarField> nutAtmWallFunctionFvPatchScalarField::calcNut() const
{
    const label patchI = patch().index();

    const turbulenceModel& turbulence =
        db().lookupObject<turbulenceModel>("turbulenceModel");
    const scalarField& y = turbulence.y()[patchI];
    const tmp<volScalarField> tk = turbulence.k();
    const volScalarField& k = tk();
    const tmp<volScalarField> tnu = turbulence.nu();
    const volScalarField& nu = tnu();
    const scalarField& nuw = nu.boundaryField()[patchI];
    
    // Velocity field on Patch (boundary field)
    const fvPatchVectorField& Uw = turbulence.U().boundaryField()[patchI]; //rca
    // Magnitude of internal field  
    //const scalarField Up = mag(Uw.patchInternalField());  //rca

    // The flow velocity at the adjacent cell centre
    const scalarField magUp(mag(Uw.patchInternalField() - Uw));  //rca
    
    //const scalarField uStarRHfield(Up*kappa_ /log( (y+z0_)/z0_));

    const scalar Cmu25 = pow025(Cmu_);

    tmp<scalarField> tnutw(new scalarField(*this));
    scalarField& nutw = tnutw();

    forAll(nutw, faceI)
    {
        label faceCellI = patch().faceCells()[faceI];
        
        scalar z0value = max(z0min,z0_[faceI]);

        //scalar EdashT = (y[faceI] + z0_[faceI]) / z0_[faceI]; //rca
		scalar EdashT = (y[faceI] + z0value) / z0value; //rca
        
        // Calculate ustar with logarithmic wall law for rough walls (Richards and Hoxley method)
        scalar uStarRH = magUp[faceI]*kappa_ / log( max(EdashT, 1 + 1e-4) );  //rca
        
        // Calculate ustar assuming Equilibrium
        scalar uStarEq = Cmu25*sqrt(k[faceCellI]);  //rca
        scalar yPlus = uStarEq*y[faceI]/nuw[faceI];
        
        //scalar tauw = uStarRH*uStarRH;  // rchavez
        scalar tauw = uStarRH*uStarEq;   //Sorensen and Sumner and RH 
        
        //nutw[faceI] = tauw*kappa_*z0_[faceI]*Edash /uStarRH - nuw[faceI]; //rca
        //nutw[faceI] = tauw*y[faceI]/(max(magUp[faceI],1e-4 ) ) - nuw[faceI];  //rca
        nutw[faceI] = max( tauw*y[faceI]/( max(magUp[faceI],1e-4) ) - nuw[faceI] , 0.0);  //rca

        
        /*  //ORIGINAL OF
        scalar uStar = Cmu25*sqrt(k[faceCellI]);
        scalar yPlus = uStar*y[faceI]/nuw[faceI];
        scalar Edash = (y[faceI] + z0_[faceI])/(z0_[faceI] + 1e-4);
        nutw[faceI] =
            nuw[faceI]*(yPlus*kappa_/log(max(Edash, 1 + 1e-4)) - 1);
        */
        
//        Info << ";  nutw[" <<faceI<<"]="<<nutw<< endl;
        Info << "Cmu=" <<Cmu_<< endl;
        Info << "kappa="<<kappa_<< endl;
    
    }

    return tnutw;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

nutAtmWallFunctionFvPatchScalarField::
nutAtmWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    nutkWallFunctionFvPatchScalarField(p, iF),
    zZeroField_   //Build z0_ Field  (R.Chavez)
    (
        new volScalarField
        (
            IOobject
            (
	            "z0",
                db().time().timeName(),
                patch().boundaryMesh().mesh(),
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            patch().boundaryMesh().mesh()
        )
    ),
    z0_(zZeroField_->boundaryField()[patch().index()])

{}


nutAtmWallFunctionFvPatchScalarField::
nutAtmWallFunctionFvPatchScalarField
(
    const nutAtmWallFunctionFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    nutkWallFunctionFvPatchScalarField(ptf, p, iF, mapper),
    z0_(ptf.z0_, mapper)
{}


nutAtmWallFunctionFvPatchScalarField::
nutAtmWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    nutkWallFunctionFvPatchScalarField(p, iF, dict),
    zZeroField_   //Build z0_ Field  (R.Chavez)
    (
        new volScalarField
        (
            IOobject
            (
                "z0",
                db().time().timeName(),
                patch().boundaryMesh().mesh(),
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            patch().boundaryMesh().mesh()
        )
    ),
    z0_(zZeroField_->boundaryField()[patch().index()])

{}


nutAtmWallFunctionFvPatchScalarField::
nutAtmWallFunctionFvPatchScalarField
(
    const nutAtmWallFunctionFvPatchScalarField& rwfpsf
)
:
    nutkWallFunctionFvPatchScalarField(rwfpsf),
    z0_(rwfpsf.z0_)
{}


nutAtmWallFunctionFvPatchScalarField::
nutAtmWallFunctionFvPatchScalarField
(
    const nutAtmWallFunctionFvPatchScalarField& rwfpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    nutkWallFunctionFvPatchScalarField(rwfpsf, iF),
    z0_(rwfpsf.z0_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void nutAtmWallFunctionFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    nutkWallFunctionFvPatchScalarField::autoMap(m);
    z0_.autoMap(m);
}


void nutAtmWallFunctionFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    nutkWallFunctionFvPatchScalarField::rmap(ptf, addr);

    const nutAtmWallFunctionFvPatchScalarField& nrwfpsf =
        refCast<const nutAtmWallFunctionFvPatchScalarField>(ptf);

    z0_.rmap(nrwfpsf.z0_, addr);
}


void nutAtmWallFunctionFvPatchScalarField::write(Ostream& os) const
{
    fvPatchField<scalar>::write(os);
    writeLocalEntries(os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    nutAtmWallFunctionFvPatchScalarField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
