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

\*---------------------------------------------------------------------------*/

#include "interpolate2D.H"
#include "primitiveFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
List<List<Type> > interpolate2D
(
    const List<scalar>& xi,
    const List<scalar>& yi,
    const List<scalar>& x,
    const List<scalar>& y,
    const List<List<Type> >& f
)
{
    // Get size of interpolation point lists.
    label nxi = xi.size();
    label nyi = yi.size();

    // Interpolate element by element.
    List<List<Type> > fi(nxi,List<Type>(nyi));
    forAll(xi, i)
    {
        forAll(yi, j)
        {
            fi[i][j] = interpolate2D(xi[i],yi[j],x,y,f);
        }
    }
    return fi;
}


template<class Type>
List<Type> interpolate2D
(
    const scalar& xi,
    const List<scalar>& yi,
    const List<scalar>& x,
    const List<scalar>& y,
    const List<List<Type> >& f
)
{
    // Get size of interpolation point lists.
    label nyi = yi.size();

    // Interpolate element by element.
    List<Type> fi(nyi);
    forAll(yi, j)
    {
        fi[j] = interpolate2D(xi,yi[j],x,y,f);
    }
    return fi;
}


template<class Type>
List<Type> interpolate2D
(
    const List<scalar>& xi,
    const scalar& yi,
    const List<scalar>& x,
    const List<scalar>& y,
    const List<List<Type> >& f
)
{
    // Get size of interpolation point lists.
    label nxi = xi.size();

    // Interpolate element by element.
    List<Type> fi(nxi);
    forAll(xi, i)
    {
        fi[i] = interpolate2D(xi[i],yi,x,y,f);
    }
    return fi;
}

scalar findIndex
(
    const scalar x,
    const label nx,
    const List<scalar>& xs
)
{
	label iLow = 0;
	label iHigh = nx-1;
	label i;
	
	// take care of out of range cases
	//if(xs[iLow] >= x) iHigh = iLow;
	//if(xs[iHigh] <= x) iLow = iHigh;
	
	while ((iHigh-iLow) > 1)
    {
        i = label((iHigh+iLow)/2);
              
        if(xs[i] < x)
        {
			iLow = i;
		}
		else if(xs[i] > x)
		{
			iHigh = i;
		}
		else 
		{
			iLow = i;
			iHigh = i;
		}

    }
    return iLow;
}

template<class Type>
Type interpolate2D
(
    const scalar xi,
    const scalar yi,
    const List<scalar>& x,
    const List<scalar>& y,
    const List<List<Type> >& f
)
{
    // Get interpolation data size.
    label nx = x.size();
    label ny = y.size();
    
 	// TODO - this step should be done when initializing, not here
/*   
    label nxf = f.size();
    label nyf = f[0].size();

    // Check to make sure data sizes all match up.  Does the size of x
    // match the x-index size of f, and same with y.  Give error message
    // and exit if not.
    if ((nx != nxf) || (ny != nyf))
    {
        FatalErrorIn
        (
            "interpolate2d"
        ) << "Sizes of input x and y vectors do not match size of f array:" << endl <<
             "x size: " << nx << endl <<
             "y size: " << ny << endl <<
             "f(x,y) size: " << nxf << ", " << nyf <<
             abort(FatalError);
    }
*/

    // Find bounding indices in x.
    label xLow = findIndex(xi,nx,x);
    label xHigh = min(xLow + 1,nx - 1);

    // Find bounding indices in y.
	label yLow = findIndex(yi,ny,y);
    label yHigh = min(yLow + 1,ny - 1);

	// When out of range we dont want linear extrapolation, instead it is safer to use the nearest value
    /*// If, the data point lies outside the x or y given data ranges
    // set up the high and low indices so that linear extrapolation
    // will be done.
    if (xHigh == xLow)
    {
        if (xHigh == 0)
        {
            xHigh++;
        }
        else if (xHigh == nx - 1)
        {
            xLow--;
        }
    }
   
    if (yHigh == yLow)
    {
        if (yHigh == 0)
        {
            yHigh++;
        }
        else if (yHigh == ny - 1)
        {
            yLow--;
        }
    }*/


    Type m1;
    Type m2;
    Type f1;
    Type f2;
    Type fi;

    if (ny > 1)
    {
        // First, interpolate in x.
        if (nx > 1)
        {
            m1 = (f[xHigh][yLow]  - f[xLow][yLow])  / (x[xHigh] - x[xLow]);
            m2 = (f[xHigh][yHigh] - f[xLow][yHigh]) / (x[xHigh] - x[xLow]);

            f1 = f[xLow][yLow]  + (m1 * (xi - x[xLow]));
            f2 = f[xLow][yHigh] + (m2 * (xi - x[xLow]));
        }
        else
        {
            f1 = f[xLow][yLow];
            f2 = f[xLow][yHigh];
        }

        // Then, interpolate in y.
        Type n = (f2 - f1) / (y[yHigh] - y[yLow]);

        fi = f1 + (n * (yi - y[yLow]));
    }
    else
    {
        // Interpolate in x only.
        if (nx > 1)
        {
            m1 = (f[xHigh][yLow] - f[xLow][yLow]) / (x[xHigh] - x[xLow]);
            fi = f[xLow][yLow] + (m1 * (xi - x[xLow]));
        }
        else
        {
            fi = f[xLow][yLow];
        }
    }

    return fi;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
