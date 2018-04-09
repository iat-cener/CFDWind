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

#include "alphatAlletoWallFunctionFvPatchScalarField.H"
#include "incompressible/turbulenceModel/turbulenceModel.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "wallFvPatch.H"
#include "addToRunTimeSelectionTable.H"

#define z0min 0.0002  //rchavez (2016)

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

//scalar alphatAlletoWallFunctionFvPatchScalarField::tolerance_ = 0.01;
//label alphatAlletoWallFunctionFvPatchScalarField::maxIters_ = 10;

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void alphatAlletoWallFunctionFvPatchScalarField::checkType()
{
    if (!isA<wallFvPatch>(patch()))
    {
        FatalErrorIn
        (
            "alphatAlletoWallFunctionFvPatchScalarField::checkType()"
        )   << "Invalid wall function specification" << nl
            << "    Patch type for patch " << patch().name()
            << " must be wall" << nl
            << "    Current patch type is " << patch().type() << nl << endl
            << abort(FatalError);
    }
}


/*
scalar alphatAlletoWallFunctionFvPatchScalarField::yPlusTherm
(
    const scalar P,
    const scalar Prat
) const
{
    scalar ypt = 11.0;

    for (int i=0; i<maxIters_; i++)
    {
        scalar f = ypt - (log(E_*ypt)/kappa_ + P)/Prat;
        scalar df = 1.0 - 1.0/(ypt*kappa_*Prat);
        scalar yptNew = ypt - f/df;

        if (yptNew < VSMALL)
        {
            return 0;
        }
        else if (mag(yptNew - ypt) < tolerance_)
        {
            return yptNew;
        }
        else
        {
            ypt = yptNew;
        }
     }

    return ypt;
}
*/

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

alphatAlletoWallFunctionFvPatchScalarField::
alphatAlletoWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    Prts_(0),
    Cmu_(0),
    kappa_(0),
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

{
    checkType();
}


alphatAlletoWallFunctionFvPatchScalarField::
alphatAlletoWallFunctionFvPatchScalarField
(
    const alphatAlletoWallFunctionFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    Prts_(ptf.Prts_),
    Cmu_(ptf.Cmu_),
    kappa_(ptf.kappa_),
	z0_(ptf.z0_) //R. Chavez
{
    checkType();
}


alphatAlletoWallFunctionFvPatchScalarField::
alphatAlletoWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF, dict),
    Prts_(readScalar(dict.lookup("Prts"))), // force read to avoid ambiguity
    //Cmu_(dict.lookupOrDefault<scalar>("Cmu", 0.0256)),
    //kappa_(dict.lookupOrDefault<scalar>("kappa", 0.4))
	Cmu_(readScalar(dict.lookup("Cmu"))),
	kappa_(readScalar(dict.lookup("kappa"))),
    zZeroField_    //Build z0_ Field  (R.Chavez)
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

{
    checkType();
}


alphatAlletoWallFunctionFvPatchScalarField::
alphatAlletoWallFunctionFvPatchScalarField
(
    const alphatAlletoWallFunctionFvPatchScalarField& wfpsf
)
:
    fixedValueFvPatchScalarField(wfpsf),
    Prts_(wfpsf.Prts_),
    Cmu_(wfpsf.Cmu_),
    kappa_(wfpsf.kappa_),
    z0_(wfpsf.z0_) // R. Chavez
{
    checkType();
}


alphatAlletoWallFunctionFvPatchScalarField::
alphatAlletoWallFunctionFvPatchScalarField
(
    const alphatAlletoWallFunctionFvPatchScalarField& wfpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(wfpsf, iF),
    Prts_(wfpsf.Prts_),
    Cmu_(wfpsf.Cmu_),
    kappa_(wfpsf.kappa_),
	z0_(wfpsf.z0_) // R. Chavez
{
    checkType();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void alphatAlletoWallFunctionFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const label patchi = patch().index();

    // Retrieve turbulence properties from model
    const turbulenceModel& turbModel =
        db().lookupObject<turbulenceModel>("turbulenceModel");
    //const scalar Cmu25 = pow(Cmu_, 0.25);
    const scalarField& y = turbModel.y()[patchi];
    const tmp<volScalarField> tnu = turbModel.nu();
    const volScalarField& nu = tnu();
    const scalarField& nuw = nu.boundaryField()[patchi];
    const tmp<volScalarField> tk = turbModel.k();
    const volScalarField& k = tk();

    const IOdictionary& transportProperties =
        db().lookupObject<IOdictionary>("transportProperties");

	//my shit
	const fvPatchVectorField& Uw = turbModel.U().boundaryField()[patchi];
    //const scalarField magGradUw(mag(Uw.snGrad()));
	// The flow velocity at the adjacent cell centre
    const scalarField magUp(mag(Uw.patchInternalField() - Uw));  //rca
	//const fvPatchScalarField& Tw = turbModel.T().boundaryField()[patchi];
	//finish my shit

    // Molecular Prandtl number
    /*const scalar Pr
    (
        dimensionedScalar(transportProperties.lookup("Pr")).value()
    );*/


    // Populate boundary values
    scalarField& alphatw = *this;
    forAll(alphatw, faceI)
    {
        label faceCellI = patch().faceCells()[faceI];
        
        //scalar z0value = 0.1;
        scalar z0value = max(z0min,z0_[faceI]);

        scalar EdashT = 1.0 + max( y[faceI]/z0value , 1e-5); //rca
        scalar uStarRH = magUp[faceI]*kappa_ / log(EdashT);  //rca
        //scalar uStarEq = Cmu25*sqrt(k[faceCellI]);  //rca
        alphatw[faceI] = uStarRH * y[faceI] * kappa_/( log(EdashT) * Prts_ );

		//Info << "La puta uStarEq="<< uStarEq<<endl;
		//Info << "La puta uStarRH="<< uStarRH<<endl;
		//Info << "La puta alphatw["<<faceI<<"]"<< alphatw[faceI]<<endl;
		//Info << "La puta magUp["<<faceI<<"]"<< magUp[faceI]<<endl;
		//Info << "La puta magGradUw["<<faceI<<"]"<< magGradUw[faceI]<<endl;
        Info << "La puta Prts"<< Prts_ <<endl;
		
		
		/*   //original Crap
        // y+
        scalar yPlus = Cmu25*sqrt(k[faceCellI])*y[faceI]/nuw[faceI];

        // Molecular-to-turbulent Prandtl number ratio
        scalar Prat = Pr/Prt_;

        // Thermal sublayer thickness
        scalar P = Psmooth(Prat);
        scalar yPlusTherm = this->yPlusTherm(P, Prat);

        // Update turbulent thermal conductivity
        if (yPlus > yPlusTherm)
        {
            scalar nu = nuw[faceI];
            scalar kt = nu*(yPlus/(Prt_*(log(E_*yPlus)/kappa_ + P)) - 1/Pr);
            alphatw[faceI] = max(0.0, kt);
        }
        else
        {
            alphatw[faceI] = 0.0;
        }
        */
    }

    fixedValueFvPatchField<scalar>::updateCoeffs();
}


void alphatAlletoWallFunctionFvPatchScalarField::write(Ostream& os) const
{
    fvPatchField<scalar>::write(os);
    os.writeKeyword("Prts") << Prts_ << token::END_STATEMENT << nl;
    os.writeKeyword("Cmu") << Cmu_ << token::END_STATEMENT << nl;
    os.writeKeyword("kappa") << kappa_ << token::END_STATEMENT << nl;
    //os.writeKeyword("E") << E_ << token::END_STATEMENT << nl;
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    alphatAlletoWallFunctionFvPatchScalarField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
