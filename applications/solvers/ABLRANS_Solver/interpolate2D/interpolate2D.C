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
    const List<scalar>& x,
    const List<scalar>& y,
    const List<scalar>& xs,
    const List<scalar>& ys,
    const List<List<Type> >& f
)
{
    // Get size of interpolation point lists.
    label nxi = x.size();
    label nyi = y.size();

    // Interpolate element by element.
    List<List<Type> > fi(nxi,List<Type>(nyi));
    forAll(x, i)
    {
        forAll(y, j)
        {
            fi[i][j] = interpolate2D(x[i],y[j],xs,ys,f);
        }
    }
    return fi;
}


template<class Type>
List<Type> interpolate2D
(
    const scalar& x,
    const List<scalar>& y,
    const List<scalar>& xs,
    const List<scalar>& ys,
    const List<List<Type> >& f
)
{
    // Get size of interpolation point lists.
    label nyi = y.size();

    // Interpolate element by element.
    List<Type> fi(nyi);
    forAll(y, j)
    {
        fi[j] = interpolate2D(x,y[j],xs,ys,f);
    }
    return fi;
}

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
)
{
    // Get size of interpolation point lists.
    label nyi = y.size();

    // Interpolate element by element.
    List<Type> fi(nyi);
    forAll(y, j)
    {
        fi[j] = interpolate2D(x,y[j],xLow,yLow[j],xs,ys,f);
    }
    return fi;
}

template<class Type>
List<Type> interpolate2D
(
    const List<scalar>& x,
    const scalar& y,
    const List<scalar>& xs,
    const List<scalar>& ys,
    const List<List<Type> >& f
)
{
    // Get size of interpolation point lists.
    label nxi = x.size();

    // Interpolate element by element.
    List<Type> fi(nxi);
    forAll(x, i)
    {
        fi[i] = interpolate2D(x[i],y,xs,ys,f);
    }
    return fi;
}


scalar findIndex
(
    const scalar x,
    const List<scalar>& xs
)
{
	//findIndex helper function O(log n)


    label nx = xs.size();
	label iLow = 0; 
	label iHigh = nx-1;

	// to speed up convergance, assume linear distribution and try to guess the correct index
	label i = label( (iHigh-iLow) * (x-xs[iLow]) / (xs[iHigh]-xs[iLow]) );
	
	// now lets check if we guessed it right
	label gHigh = min(i + 1,nx - 1);
	label gLow = max(i - 1, 0);
	if((x >= xs[gLow]) && (x <= xs[gHigh])){
		//we got the solution, now just refine it if necessarry
		if(xs[i] <= x) return i;
		return gLow;
	}

	while ((iHigh-iLow) > 1)
    {     
    	if(xs[i] < x)
				iLow = i;
		else if(xs[i] > x)
				iHigh = i;
		else 
		{
			iLow = i;
			iHigh = i;
		}
		i = label((iHigh+iLow)/2);
    }

	return iLow;
}


template<class Type>
Type interpolate2D
(
    const scalar x,
    const scalar y,
    const List<scalar>& xs,
    const List<scalar>& ys,
    const List<List<Type> >& f
)
{
    // Find low bounding indices in x.
	label xIndex = findIndex(x,xs);

    // Find low bounding indices in y.
	label yIndex = findIndex(y,ys);


    return interpolate2D(x,y,xIndex,yIndex,xs,ys,f);
}


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
)
{
    // Get interpolation data size.
    label nx = xs.size();
    label ny = ys.size();
    
    // Find high bounding indices in x.
    label xHigh = min(xLow + 1,nx - 1);

    // Find high bounding indices in y.
    label yHigh = min(yLow + 1,ny - 1);

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
            m1 = (f[xHigh][yLow]  - f[xLow][yLow])  / (xs[xHigh] - xs[xLow]);
            m2 = (f[xHigh][yHigh] - f[xLow][yHigh]) / (xs[xHigh] - xs[xLow]);

            f1 = f[xLow][yLow]  + (m1 * (x - xs[xLow]));
            f2 = f[xLow][yHigh] + (m2 * (x - xs[xLow]));
        }
        else
        {
            f1 = f[xLow][yLow];
            f2 = f[xLow][yHigh];
        }

        // Then, interpolate in y.
        Type n = (f2 - f1) / (ys[yHigh] - ys[yLow]);

        fi = f1 + (n * (y - ys[yLow]));
    }
    else
    {
        // Interpolate in x only.
        if (nx > 1)
        {
            m1 = (f[xHigh][yLow] - f[xLow][yLow]) / (xs[xHigh] - xs[xLow]);
            fi = f[xLow][yLow] + (m1 * (x - xs[xLow]));
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
