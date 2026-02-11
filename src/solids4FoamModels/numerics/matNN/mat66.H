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
    Foam::mat66

  Description
      Lightweight fixed-size 6x6 matrix type (row-major storage) for material
      tangent operators or similar data in Voigt notation.

      - Provides (i,j) element accessors.
      - Provides m() and n() methods for compatibility with existing code.
      - Includes helpers to zero, fill with a constant, or set to identity.
      - Supports initialisation from an initializer_list of up to 36 values.

      Designed for use in List<mat66> containers as a portable replacement for
      scalarSquareMatrix when a fixed 6x6 matrix is sufficient. Avoids the
      need for pTraits and tmp<>, and works across foam-extend, OpenFOAM-
      foundation, and OpenFOAM-esi forks.

Author
    Philip Cardiff, UCD.

\*---------------------------------------------------------------------------*/

#ifndef mat66_H
#define mat66_H

#include "scalar.H"
#include "label.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

struct mat66
{
    // Storage in row-major order: a[i*6 + j]
    scalar a_[36];

    static const label dim_ = 6;

    // Element access
    inline scalar& operator()(label i, label j)
    {
#ifdef FULLDEBUG
        if (i >= dim_ || j >= dim_)
        {
            FatalErrorInFunction
                << "Out of bounds" << exit(FatalError);
        }
#endif

        return a_[i*dim_ + j];
    }
    inline const scalar& operator()(label i, label j) const
    {
#ifdef FULLDEBUG
        if (i >= dim_ || j >= dim_)
        {
            FatalErrorInFunction
                << "Out of bounds" << exit(FatalError);
        }
#endif

        return a_[i*dim_ + j];
    }

    // Dimensions
    inline label m() const
    {
        return dim_;
    }
    inline label n() const
    {
        return dim_;
    }

    // Set all entries to a value
    inline void fill(scalar val)
    {
        for (label k=0; k<36; ++k)
        {
            a_[k] = val;
        }
    }

    // Zero matrix
    inline void clear()
    {
        fill(0);
    }

    // Identity matrix
    inline void identity()
    {
        clear();
        for (label i=0; i<dim_; ++i)
        {
            (*this)(i,i) = 1;
        }
    }

    // Construct from initializer_list, e.g. {1,2,...,36}
    mat66(std::initializer_list<scalar> vals)
    {
        label n = std::min<label>(36, vals.size());
        auto it = vals.begin();
        for (label k=0; k<n; ++k, ++it)
        {
            a_[k] = *it;
        }
        for (label k=n; k<36; ++k)
        {
            a_[k] = 0; // pad with zero if fewer values
        }
    }

    // Default constructor leaves values uninitialised
    mat66() = default;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
