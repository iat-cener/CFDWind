/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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
    Foam::incompressible::nutAtmWallFunctionFvPatchScalarField

Group
    grpIcoWallFunctions

Description
    This boundary condition provides a turbulent kinematic viscosity for
    atmospheric velocity profiles.  It is desinged to be used in conjunction
    with the atmBoundaryLayerInletVelocity boundary condition.  The values
    are calculated using:

        \f[
            U = frac{U_f}{K} ln(\frac{z + z_0}{z_0})
        \f]

    where

    \vartable
        U_f | frictional velocity
        K   | Von Karman's constant
        z_0 | surface roughness length
        z   | vertical co-ordinate
    \endvartable

    \heading Patch usage

    \table
        Property     | Description             | Required    | Default value
        z0           | surface roughness length| yes         |
    \endtable

    Example of the boundary condition specification:
    \verbatim
    myPatch
    {
        type            nutAtmWallFunction;
        z0              uniform 0;
    }
    \endverbatim

SeeAlso
    Foam::nutkWallFunctionFvPatchField

SourceFiles
    nutAtmWallFunctionFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef nutAtmWallFunctionFvPatchScalarField_H
#define nutAtmWallFunctionFvPatchScalarField_H

#include "nutkWallFunctionFvPatchScalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{

/*---------------------------------------------------------------------------*\
          Class nutAtmWallFunctionFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class nutAtmWallFunctionFvPatchScalarField
:
    public nutkWallFunctionFvPatchScalarField
{
protected:

    // Protected data

        //- z0 volScalarField
        // A new field is introduced in order to easily  consider varying 
        // roughness.  rchavez
        volScalarField* zZeroField_;   

        //- Surface roughness length
        scalarField z0_;  //R. Chavez
	    
	    
    // Protected Member Functions


        //- Calculate the turbulence viscosity
        virtual tmp<scalarField> calcNut() const;


public:

    //- Runtime type information
    TypeName("nutAtmWallFunction");


    // Constructors

        //- Construct from patch and internal field
        nutAtmWallFunctionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        nutAtmWallFunctionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  nutAtmWallFunctionFvPatchScalarField
        //  onto a new patch
        nutAtmWallFunctionFvPatchScalarField
        (
            const nutAtmWallFunctionFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        nutAtmWallFunctionFvPatchScalarField
        (
            const nutAtmWallFunctionFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new nutAtmWallFunctionFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        nutAtmWallFunctionFvPatchScalarField
        (
            const nutAtmWallFunctionFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new nutAtmWallFunctionFvPatchScalarField(*this, iF)
            );
        }


    // Member functions

        // Access functions

            // Return z0
            scalarField& z0()
            {
                return z0_;
            }

        

        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap(const fvPatchFieldMapper&);

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap
            (
                const fvPatchScalarField&,
                const labelList&
            );


        // I-O

            //- Write
            virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace incompressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
