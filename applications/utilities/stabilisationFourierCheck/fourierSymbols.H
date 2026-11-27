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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

//- Analytical eigenvalue (Fourier symbol) of a stabilisation model for the
//  mode theta on a uniform Cartesian mesh of spacings h and cell counts N.
//  It is written from the mathematical operators, independently of the
//  model implementations, so that it can verify them.
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
        q[d] = sqr(Foam::sin(theta[d]/2));
    }

    const Switch normalise
    (
        modelDict.lookupOrDefault<Switch>("spectralNormalisation", false)
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
            FatalErrorInFunction
                << "No active mesh direction"
                << exit(FatalError);
        }

        if (type == "diffStencilLaplacian" || type == "RhieChow")
        {
            // Directionally split m = 1 (bi-Laplacian) shape, normalised
            // by its raw reference response 4*rStar
            return -s*split/(rStar*sqr(hRef));
        }

        label nLaplacianFactors = 1;
        if (type == "JamesonSchmidtTurkel")
        {
            nLaplacianFactors = 2;
        }
        else if (type == "generalisedEvenOrderLaplacian")
        {
            const label laplacianPower
            (
                readLabel(modelDict.lookup("laplacianPower"))
            );
            if (laplacianPower < 0)
            {
                FatalErrorInFunction
                    << "Negative laplacianPower"
                    << exit(FatalError);
            }
            nLaplacianFactors = laplacianPower + 1;
        }
        else if (type != "laplacian")
        {
            FatalErrorInFunction
                << "Unsupported normalised symbol: " << type
                << exit(FatalError);
        }

        // Coupled repeated Laplacian with m + 1 factors; for m = 0 the
        // normalised model uses the plain snGrad stencil, without the |d|^2
        // of the raw branch. At the reference mode, A = rStar, and the
        // eigenvalue is -s/hRef^2
        return -s*Foam::pow(A/rStar, scalar(nLaplacianFactors))/sqr(hRef);
    }

    // Raw (unnormalised) symbols
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
            FatalErrorInFunction
                << "Negative laplacianPower"
                << exit(FatalError);
        }
    }
    else if (type != "laplacian")
    {
        FatalErrorInFunction
            << "Unsupported raw symbol: " << type
            << exit(FatalError);
    }

    // m = 0 carries the legacy face |d|^2 = h^2 factor, which cancels the
    // 1/h^2 of the Laplacian symbol
    if (m == 0)
    {
        return -4*s*(q.x() + q.y() + q.z());
    }

    // A is the symbol magnitude of the Laplacian; B that of the outer
    // Laplacian with its face scaling h_d^(2m)
    scalar A = 0;
    scalar B = 0;
    for (direction d = 0; d < 3; ++d)
    {
        A += 4*q[d]/sqr(h[d]);
        B += 4*q[d]*Foam::pow(h[d], scalar(2*m - 2));
    }
    return -s*Foam::pow(A, scalar(m))*B;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#endif

// ************************************************************************* //
