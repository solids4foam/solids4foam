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
    weightFunction

Description
    Helper class for evaluating weight functions used in moving least squares.

    The class reads the selected weighting kernel and its coefficients from a
    dictionary and provides a fast runtime evaluation of the weight value.

    Supported weight functions:

        one                            Constant weight
        linear                         Linear weighting
        inverseDistance                Inverse-distance weighting
        radiallySymmetricExponential   Radially symmetric exponential kernel

SourceFile
    weightFunction.C

Author
    Ivan Batistic, UCD. All rights reserved.
    Philip Cardiff, UCD. All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef weightFunction_H
#define weightFunction_H

#include "dictionary.H"
#include "Enum.H"
#include "scalar.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class weightFunction
{

public:
        //- Weight functions types
        enum weightFunctionType
        {
            ONE,
            LINEAR,
            INV_DIST,
            RAD_SYMM_EXP
        };

        //- Weight funtions names
        static const Enum<weightFunctionType> weightFunctionNames_;

private:

    // Private data

        //- Selected weight function type
        weightFunctionType type_;

        //- Kernel parameter for radial symmetric exponential kernel
        const scalar k_;

        //- Storing exp(-k^2) to avoid computing each time weight func is used
        const scalar expSqrK_;

        //- Storing 1/(1-expSqrK_) to avoid computing each time
        const scalar invDenom_;

        //- First parameter for inverse-distance weighting
        const scalar s_;

        //- Second parameter for inverse-distance weighting
        const scalar b_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        weightFunction(const weightFunction&) = delete;

        //- Disallow default bitwise assignment
        void operator=(const weightFunction&) = delete;


public:

    // Constructors
    explicit weightFunction (const dictionary& dict);

    // Destructor
    ~weightFunction() = default;

    // Member Functions

        //- Return weight for distance d and characteristic scale
        scalar weight(const scalar d, const scalar scale) const;

        //- Return selected type
        inline weightFunctionType type() const
        {
            return type_;
        }

        //- Return kernel parameter
        inline scalar k() const
        {
            return k_;
        }

        //- Return first inverse-distance parameter
        inline scalar s() const
        {
            return s_;
        }

        //- Return second inverse-distance parameter
        inline scalar b() const
        {
            return b_;
        }
};

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
