/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     |
    \\  /    A nd           | For copyright notice see file Copyright
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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
    Foam::fv::pointGaussLeastSquaresGrad

Description
    Point Gauss gradient calculation where the cell values are interpolated to
    the points using a least squares method.

SourceFiles
    pointGaussLeastSquaresGrad.C

Author
    Zeljko Tukovic Zagreb/UCD
    Re-organised by Philip Cardiff UCD/UT

\*---------------------------------------------------------------------------*/

#ifndef pointGaussLeastSquaresGrad_H
#define pointGaussLeastSquaresGrad_H

#include "gradScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

/*---------------------------------------------------------------------------*\
                      Class pointGaussLeastSquaresGrad Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class pointGaussLeastSquaresGrad
:
    public fv::gradScheme<Type>
{
    // Private Member Functions

        //- Disallow default bitwise copy construct
        pointGaussLeastSquaresGrad(const pointGaussLeastSquaresGrad&);

        //- Disallow default bitwise assignment
        void operator=(const pointGaussLeastSquaresGrad&);


public:

    //- Runtime type information
    TypeName("pointGaussLeastSquares");


    // Constructors

        //- Construct from mesh
        pointGaussLeastSquaresGrad(const fvMesh& mesh)
        :
            gradScheme<Type>(mesh)
        {}

        //- Construct from Istream
        pointGaussLeastSquaresGrad(const fvMesh& mesh, Istream&)
        :
            gradScheme<Type>(mesh)
        {}


    // Member Functions

        //- Calculate and return the grad of the given field.
        //  Used by grad either to recalculate the cached gradient when it is
        //  out of date with respect to the field or when it is not cached.
        //- Volume gradient
        virtual tmp
        <
            GeometricField
            <typename outerProduct<vector, Type>::type, fvPatchField, volMesh>
        > calcGrad
        (
            const GeometricField<Type, fvPatchField, volMesh>&,
            const word& name
        ) const;

        // Surface gradient
        // tmp
        // <
        //     GeometricField
        //     <
        //         typename outerProduct<vector, Type>::type,
        //         fvsPatchField,
        //         surfaceMesh
        //     >
        // > fGrad
        // (
        //     const GeometricField<Type, fvPatchField, volMesh>&
        // ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "pointGaussLeastSquaresGrad.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
