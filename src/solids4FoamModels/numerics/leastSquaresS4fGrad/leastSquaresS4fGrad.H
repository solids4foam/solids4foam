/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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
    Foam::fv::leastSquaresS4fGrad

Group
    grpFvGradSchemes

Description
    Second-order gradient scheme using least-squares.

SourceFiles
    leastSquaresS4fGrad.C

\*---------------------------------------------------------------------------*/

#ifndef leastSquaresS4fGrad_H
#define leastSquaresS4fGrad_H

#include "gradScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

/*---------------------------------------------------------------------------*\
                      Class leastSquaresS4fGrad Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class leastSquaresS4fGrad
:
    public fv::gradScheme<Type>
{
    // Private Member Functions

        //- No copy construct
        leastSquaresS4fGrad(const leastSquaresS4fGrad&) = delete;

        //- No copy assignment
        void operator=(const leastSquaresS4fGrad&) = delete;


public:

    //- Runtime type information
    TypeName("leastSquaresS4f");


    // Constructors

        //- Construct from mesh
        leastSquaresS4fGrad(const fvMesh& mesh)
        :
            gradScheme<Type>(mesh)
        {}

        //- Construct from Istream
        leastSquaresS4fGrad(const fvMesh& mesh, Istream&)
        :
            gradScheme<Type>(mesh)
        {}


    // Member Functions

        //- Return the gradient of the given field to the gradScheme::grad
        //  for optional caching
        virtual tmp
        <
            GeometricField
            <typename outerProduct<vector, Type>::type, fvPatchField, volMesh>
        > calcGrad
        (
            const GeometricField<Type, fvPatchField, volMesh>& vsf,
            const word& name
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "leastSquaresS4fGrad.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
