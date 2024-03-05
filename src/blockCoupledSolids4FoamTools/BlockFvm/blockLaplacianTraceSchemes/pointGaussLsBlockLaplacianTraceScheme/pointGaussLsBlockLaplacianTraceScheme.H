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
    Foam::fv::pointGaussLsDivSigmaScheme

Description
    Second order discretisation of
    laplacianTrace(mu, U) == div(lambda*I*tr(gradU))

    The normal derivatives are discretised using central differencing, and the
    tangential derivatives are discretised using the face Gauss theorem where
    the vertex values are interpolated from the surrounding cell centres using
    a linear least squares fit.

SourceFiles
    pointGaussLsBlockLaplacianTraceScheme.C
    pointGaussLsBlockLaplacianTraceSchemeTang.C

Author
    Philip Cardiff UCD

\*---------------------------------------------------------------------------*/

#ifndef pointGaussLsBlockLaplacianTraceScheme_H
#define pointGaussLsBlockLaplacianTraceScheme_H

#include "blockLaplacianTraceScheme.H"
#include "newLeastSquaresVolPointInterpolation.H"
#include "solidPolyMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

/*---------------------------------------------------------------------------*\
               Class pointGaussLsBlockLaplacianTraceScheme Declaration
\*---------------------------------------------------------------------------*/

class pointGaussLsBlockLaplacianTraceScheme
:
    public fv::blockLaplacianTrace
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
        pointGaussLsBlockLaplacianTraceScheme
        (
            const pointGaussLsBlockLaplacianTraceScheme&
        );

        //- Disallow default bitwise assignment
        void operator=(const pointGaussLsBlockLaplacianTraceScheme&);


public:

    //- Runtime type information
    TypeName("pointGaussLeastSquaresLaplacianTrace");


    // Constructors

        //- Construct from Istream
        pointGaussLsBlockLaplacianTraceScheme
        (
            const fvMesh& mesh, Istream& is
        );

    // Destructor

        virtual ~pointGaussLsBlockLaplacianTraceScheme()
        {}


    // Member Functions

        // Implicit laplacian transpose i.e. div(lambda*I*tr(gradU))
        tmp<BlockLduMatrix<vector> > fvmBlockLaplacianTrace
        (
            const solidPolyMesh& solidMesh,
            const surfaceScalarField& lambdaf,
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
