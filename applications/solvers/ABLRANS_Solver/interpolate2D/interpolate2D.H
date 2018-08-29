/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

InNamespace
    Foam

Description
    Interpolates function values from a surface to a given point on the 
    surface using bilinear interpolation.

SourceFiles
    interpolate2D.C

\*---------------------------------------------------------------------------*/

#ifndef interpolate2D_H
#define interpolate2D_H

#include "scalar.H"
#include "primitiveFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
List<List<Type> > interpolate2D
(
    const List<scalar>& x,
    const List<scalar>& y,
    const List<scalar>& xs,
    const List<scalar>& ys,
    const List<List<Type> >& f
);

template<class Type>
List<Type> interpolate2D
(
    const scalar& x,
    const List<scalar>& y,
    const List<scalar>& xs,
    const List<scalar>& ys,
    const List<List<Type> >& f
);

template<class Type>
List<Type> interpolate2D
(
    const scalar& x,
    const List<scalar>& y,
    const label& xLow,
    const List<label>& yLow,
    const List<scalar>& xs,
    const List<scalar>& ys,
    const List<List<Type> >& f
);

template<class Type>
List<Type> interpolate2D
(
    const List<scalar>& x,
    const scalar& y,
    const List<scalar>& xs,
    const List<scalar>& ys,
    const List<List<Type> >& f
);

template<class Type>
Type interpolate2D
(
    const scalar x,
    const scalar y,
    const List<scalar>& xs,
    const List<scalar>& ys,
    const List<List<Type> >& f
);


// interpolation with precomuted indexes for xs and ys lists
template<class Type>
Type interpolate2D
(
    const scalar x,
    const scalar y,
    const label xLow,
    const label yLow,
    const List<scalar>& xs,
    const List<scalar>& ys,
    const List<List<Type> >& f
);




// helper function O(log n)
scalar findIndex
(
    const scalar x,
    const label nx,
    const List<scalar>& xs
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "interpolate2D.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
