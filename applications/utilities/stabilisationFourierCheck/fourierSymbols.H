/*---------------------------------------------------------------------------*\
License
    This file is part of solids4foam.

    solids4foam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    solids4foam is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with solids4foam.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#ifndef fourierSymbols_H
#define fourierSymbols_H

#include "dictionary.H"
#include "Switch.H"
#include "vector.H"
#include <cmath>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

//- Independent mathematical symbols of the Cartesian operators.
//  No production discretisation or field evaluation belongs in this function.
inline scalar expectedEigenvalue
(
    const dictionary& modelDict,
    const vector& theta,
    const vector& h,
    const label (&N)[3]
)
{
    const word type(modelDict.lookup("type"));
    const scalar s = type == "laplacian"
      ? modelDict.lookupOrDefault<scalar>("scaleFactor", 1.0)
      : readScalar(modelDict.lookup("scaleFactor"));

    vector q;
    for (direction d = 0; d < 3; ++d)
    {
        q[d] = sqr(std::sin(theta[d]/2));
    }

    const Switch normalise
    (
        modelDict.lookupOrDefault<Switch>("normalise", false)
    );
    if (normalise)
    {
        const label rStar
        (
            readLabel(modelDict.lookup("referenceNyquistDirections"))
        );
        if (rStar <= 0)
        {
            FatalErrorInFunction
                << "referenceNyquistDirections must be positive"
                << exit(FatalError);
        }

        scalar hRef = -1;
        scalar A = 0;
        scalar split = 0;
        for (direction d = 0; d < 3; ++d)
        {
            if (N[d] > 1)
            {
                if (hRef < 0)
                {
                    hRef = h[d];
                }
                else if (mag(h[d] - hRef) > SMALL*max(h[d], hRef))
                {
                    FatalErrorInFunction
                        << "Normalised symbols require isotropic active spacing"
                        << exit(FatalError);
                }
            }
            A += q[d];
            split += sqr(q[d]);
        }
        if (hRef <= 0)
        {
            FatalErrorInFunction << "No active mesh direction"
                << exit(FatalError);
        }

        if (type == "diffStencilLaplacian" || type == "RhieChow")
        {
            // Operator-specific split-m=2 spectral shape. This is not the
            // coupled paper-m=1 formula, despite sharing the factor 1/(4r*).
            return -s*split/(rStar*sqr(hRef));
        }

        label paperOrder = 1;
        if (type == "JamesonSchmidtTurkel")
        {
            paperOrder = 2;
        }
        else if (type == "generalisedEvenOrderLaplacian")
        {
            const label laplacianPower
            (
                readLabel(modelDict.lookup("laplacianPower"))
            );
            if (laplacianPower < 0)
            {
                FatalErrorInFunction << "Negative laplacianPower"
                    << exit(FatalError);
            }
            paperOrder = laplacianPower + 1;
        }
        else if (type != "laplacian")
        {
            FatalErrorInFunction << "Unsupported normalised symbol: " << type
                << exit(FatalError);
        }

        // The normalised m=1 path is the ideal Laplacian stencil, not the
        // legacy raw h^2-scaled branch.
        return -s*std::pow(A/rStar, paperOrder)/sqr(hRef);
    }

    // The raw formulas below are the trusted legacy verification path.
    if (type == "diffStencilLaplacian" || type == "RhieChow")
    {
        scalar split = 0;
        for (direction d = 0; d < 3; ++d)
        {
            split += sqr(q[d])/sqr(h[d]);
        }
        return -4*s*split;
    }

    label m = 0;
    if (type == "JamesonSchmidtTurkel")
    {
        m = 1;
    }
    else if (type == "generalisedEvenOrderLaplacian")
    {
        m = readLabel(modelDict.lookup("laplacianPower"));
        if (m < 0)
        {
            FatalErrorInFunction << "Negative laplacianPower"
                << exit(FatalError);
        }
    }
    else if (type != "laplacian")
    {
        FatalErrorInFunction << "Unsupported raw symbol: " << type
            << exit(FatalError);
    }

    // m=0 is a special production branch, with an extra numerical h^2.
    if (m == 0)
    {
        return -4*s*(q.x() + q.y() + q.z());
    }

    // A is |lambda_Delta|; B is the direction-dependent outer face scaling.
    // This includes the anisotropic form without copying a stencil.
    scalar A = 0;
    scalar B = 0;
    for (direction d = 0; d < 3; ++d)
    {
        A += 4*q[d]/sqr(h[d]);
        B += 4*q[d]*std::pow(h[d], 2*m - 2);
    }
    return -s*std::pow(A, m)*B;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#endif

// ************************************************************************* //
