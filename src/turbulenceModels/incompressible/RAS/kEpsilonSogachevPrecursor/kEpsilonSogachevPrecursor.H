/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2014 OpenFOAM Foundation
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

Class
    Foam::incompressible::RASModels::kEpsilonSogachevPrecursor

Group
    grpIcoRASTurbulence

Description
    Standard k-epsilon turbulence model for incompressible flows.

    Reference:
    \verbatim
        Launder, B.E., and Spalding, D.B.,
        "Mathematical Models of Turbulence",
        Academic Press, 1972.

        Launder, B.E., and Spalding, D.B.,
        "Computational Methods for Turbulent Flows",
        Comp. Meth in Appl Mech & Eng'g, Vol 3, 1974, pp 269-289.
    \endverbatim

    The default model coefficients correspond to the Dettering and Ettling set as follows:
    \verbatim
        kEpsilonSogachevPrecursorCoeffs
        {
            Cmu         0.0256;
            C1          1.13;
            C2          1.9;
            sigmaEps    1.2987;
            sigmaK      0.7407;
        }
    \endverbatim

SourceFiles
    kEpsilonSogachevPrecursor.C

\*---------------------------------------------------------------------------*/

#ifndef kEpsilonSogachevPrecursor_H
#define kEpsilonSogachevPrecursor_H

#include "RASModel.H"
#include "uniformDimensionedFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
                           Class kEpsilonSogachevPrecursor Declaration
\*---------------------------------------------------------------------------*/

class kEpsilonSogachevPrecursor
:
    public RASModel
{

protected:

    // Protected data

        // Model coefficients
            dimensionedScalar Cmu_;
            dimensionedScalar C1_;
            dimensionedScalar C2_;
            dimensionedScalar sigmaEps_;
            dimensionedScalar sigmaK_;             // Not considered in v2.4 version (R. Chavez)
            dimensionedScalar kappa_;             // Not considered in v2.4 version (R. Chavez)

        // ABL dictionaries data (R. Chavez)
            const IOdictionary ABLProperties_;              // ABL dictionary of input data
            const IOdictionary transportProperties_;        // transport properties dictionary

        // Fields
            volScalarField k_;
            volScalarField epsilon_;
            volScalarField nut_;
            volScalarField lMY_;
            volScalarField lm_;
            volScalarField C1ast_;
            volScalarField alphaB_;
            const dimensionedScalar Prts_;
            volScalarField Prt_;
            volScalarField C3_;
            volScalarField Rig_;
            volScalarField RiG_;
            volScalarField gradTz_;
            volScalarField wTz_;
            volScalarField uStar_;
            volScalarField L_;

            const uniformDimensionedVectorField& g_; // gravitational aceleration
            Switch strongLimitForEpsilon_;                //
            Switch limitTopAmbientValues_;
            const scalar minHeightForAmbientValues_;
            const dimensionedScalar lmAmb_;
            const dimensionedScalar kAmb_;         // clip for k
            const dimensionedScalar epsilonAmb_;    // maximum lm
            Switch calculatefc_;                      // compute fc or it is provided??
            scalar latitude_;                         // latiude

            const dimensionedScalar omegaEarth_;     // Earth's rotation rate(7.2921e−5 rad/s)
            dimensionedScalar fc_;                    // Coriolis factor
            dimensionedVector Ug_;                    // Geostrophic Wind
            dimensionedScalar lmax_;
            const dimensionedScalar TRef_;                  // Reference potentital temperature
            const vector lMYreferencePoint_;                 // Base point for k integration for mellor-yamada lm
            const label iCellMY;
            const vector zDir_;
            const dimensionedScalar alphaBcoeff_;

            const volVectorField cellCentres_;
            const volScalarField z_;

        // Private member functions (R. Chavez)
           // tmp<volScalarField> calclm(); //compute mixing length

public:

    //- Runtime type information
    TypeName("kEpsilonSogachevPrecursor");

    // Constructors

        //- Construct from components
        kEpsilonSogachevPrecursor
        (
            const volVectorField& U,
            const surfaceScalarField& phi,
            transportModel& transport,
            const word& turbulenceModelName = turbulenceModel::typeName,
            const word& modelName = typeName
        );


    //- Destructor
    virtual ~kEpsilonSogachevPrecursor()
    {}


    // Member Functions

        //- Return the turbulence viscosity
        virtual tmp<volScalarField> nut() const
        {
            return nut_;
        }


        //- Return the turbulent Prandtl number
        virtual tmp<volScalarField> Prt() const
        {
            return Prt_;
        }

        //- Return the effective diffusivity for k
        tmp<volScalarField> DkEff() const
        {
            return tmp<volScalarField>
            (
                //new volScalarField("DkEff", nut_ + nu())
                new volScalarField("DkEff", nut_/sigmaK_ + nu())   //R. Chavez
            );
        }

        //- Return the effective diffusivity for epsilon
        tmp<volScalarField> DepsilonEff() const
        {
            return tmp<volScalarField>
            (
                new volScalarField("DepsilonEff", nut_/sigmaEps_ + nu())
            );
        }

        //- Return the turbulence kinetic energy
        virtual tmp<volScalarField> k() const
        {
            return k_;
        }

        //- Return the turbulence kinetic energy dissipation rate
        virtual tmp<volScalarField> epsilon() const
        {
            return epsilon_;
        }

        //- Return the Reynolds stress tensor
        virtual tmp<volSymmTensorField> R() const;

        //- Return the effective stress tensor including the laminar stress
        virtual tmp<volSymmTensorField> devReff() const;

        //- Return the source term for the momentum equation
        virtual tmp<fvVectorMatrix> divDevReff(volVectorField& U) const;

        //- Return the source term for the momentum equation
        virtual tmp<fvVectorMatrix> divDevRhoReff
        (
            const volScalarField& rho,
            volVectorField& U
        ) const;

        //- Solve the turbulence equations and correct the turbulence viscosity
        virtual void correct();

        //- Read RASProperties dictionary
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace incompressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
