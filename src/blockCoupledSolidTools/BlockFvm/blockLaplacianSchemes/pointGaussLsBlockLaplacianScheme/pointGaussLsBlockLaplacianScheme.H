/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     |
    \\  /    A nd           | For copyright notice see file Copyright
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::fv::pointGaussLsDivSigmaScheme

Description
    Second order discretisation of laplacian(U) where non-orthogonal correction
    is implicitly discretised within the block matrix

    The normal derivatives are discretised using central differencing, and the
    tangential derivatives are discretised using the face Gauss theorem where
    the vertex values are interpolated from the surrounding cell centres using
    a linear least squares fit.

SourceFiles
    pointGaussLsBlockLaplacianScheme.C
    pointGaussLsBlockLaplacianSchemeTang.C

Author
    Philip Cardiff UCD

\*---------------------------------------------------------------------------*/

#ifndef pointGaussLsBlockLaplacianScheme_H
#define pointGaussLsBlockLaplacianScheme_H

#include "blockLaplacianScheme.H"
#include "newLeastSquaresVolPointInterpolation.H"
#include "solidPolyMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

/*---------------------------------------------------------------------------*\
               Class pointGaussLsBlockLaplacianScheme Declaration
\*---------------------------------------------------------------------------*/

class pointGaussLsBlockLaplacianScheme
:
    public fv::blockLaplacian
{
    // Private Member Functions

        // Return const reference to the volume-to-point interpolator
        const newLeastSquaresVolPointInterpolation& volToPointInterp() const
        {
            // return volToPointInterp_;
            return mesh().thisDb().lookupObject<solidPolyMesh>
            (
                "solidPolyMesh"
            ).volToPointInterp();
        }

        //- Insert coefficient contribution from normal derivative terms
        void insertCoeffsNorm
        (
            const solidPolyMesh& solidMesh,
            const surfaceScalarField& muf,
            const GeometricField<vector, fvPatchField, volMesh>& U,
            Field<vector>& blockB,
            BlockLduMatrix<vector>& blockM
        );

        //- Insert coefficient contribution from tangential derivative terms
        void insertCoeffsTang
        (
            const solidPolyMesh& solidMesh,
            const surfaceScalarField& muf,
            const GeometricField<vector, fvPatchField, volMesh>& U,
            Field<vector>& blockB,
            BlockLduMatrix<vector>& blockM
        );

        //- Insert coefficient contribution from boundary conditions
        // WE WILL NOT ADD BCS HERE
        // void insertCoeffsBc
        // (
        //     const solidPolyMesh& solidMesh,
        //     const surfaceScalarField& muf,
        //     const GeometricField<vector, fvPatchField, volMesh>& U,
        //     Field<vector>& blockB,
        //     BlockLduMatrix<vector>& blockM
        // );

        //- Disallow default bitwise copy construct
        pointGaussLsBlockLaplacianScheme
        (
            const pointGaussLsBlockLaplacianScheme&
        );

        //- Disallow default bitwise assignment
        void operator=(const pointGaussLsBlockLaplacianScheme&);


public:

    //- Runtime type information
    TypeName("pointGaussLeastSquares");


    // Constructors

        //- Construct from Istream
        pointGaussLsBlockLaplacianScheme(const fvMesh& mesh, Istream& is);

    // Destructor

        virtual ~pointGaussLsBlockLaplacianScheme()
        {}


    // Member Functions

        // Implicit laplacian with implicit non-orthogonal correction
        tmp<BlockLduMatrix<vector> > fvmBlockLaplacian
        (
            const solidPolyMesh& solidMesh,
            const surfaceScalarField& muf,
            GeometricField<vector, fvPatchField, volMesh>& U,
            Field<vector>& blockB
        );
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
