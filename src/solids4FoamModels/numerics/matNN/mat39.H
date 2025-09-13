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

Class
    Foam::mat39

  Description
      Lightweight fixed-size 3x9 matrix type (row-major storage) for material
      tangent operators or similar data in Voigt notation.

      - Provides (i,j) element accessors.
      - Provides m() and n() methods for compatibility with existing code.
      - Includes helpers to zero, fill with a constant, or set to identity.
      - Supports initialisation from an initializer_list of up to 27 values.

      Designed for use in List<mat39> containers as a portable replacement for
      scalarSquareMatrix when a fixed 3x9 matrix is sufficient. Avoids the
      need for pTraits and tmp<>, and works across foam-extend, OpenFOAM-
      foundation, and OpenFOAM-esi forks.

Author
    Philip Cardiff, UCD.

\*---------------------------------------------------------------------------*/

#ifndef mat39_H
#define mat39_H

#include "scalar.H"
#include "label.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

struct mat39
{
    // Storage in row-major order: a[i*9 + j]
    scalar a_[27];

    static const label nRows_ = 3;
    static const label nCols_ = 9;

    // Element access
    inline scalar& operator()(label i, label j)
    {
#ifdef FULLDEBUG
        if (i >= nRows_ || j >= nCols_)
        {
            FatalErrorInFunction
                << "Out of bounds" << exit(FatalError);
        }
#endif

        return a_[i*nCols_ + j];
    }
    inline const scalar& operator()(label i, label j) const
    {
#ifdef FULLDEBUG
        if (i >= nRows_ || j >= nCols_)
        {
            FatalErrorInFunction
                << "Out of bounds" << exit(FatalError);
        }
#endif

        return a_[i*nCols_ + j];
    }

    // Dimensions
    inline label m() const
    {
        return nRows_;
    }
    inline label n() const
    {
        return nCols_;
    }

    // Set all entries to a value
    inline void fill(scalar val)
    {
        for (label k=0; k<27; ++k)
        {
            a_[k] = val;
        }
    }

    // Zero matrix
    inline void clear()
    {
        fill(0);
    }

    // Construct from initializer_list, e.g. {1,2,...,27}
    mat39(std::initializer_list<scalar> vals)
    {
        label n = std::min<label>(27, vals.size());
        auto it = vals.begin();
        for (label k=0; k<n; ++k, ++it)
        {
            a_[k] = *it;
        }
        for (label k=n; k<27; ++k)
        {
            a_[k] = 0; // pad with zero if fewer values
        }
    }

    // Default constructor leaves values uninitialised
    mat39() = default;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
