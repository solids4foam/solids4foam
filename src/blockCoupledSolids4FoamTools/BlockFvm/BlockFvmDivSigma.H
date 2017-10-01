/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

InNamespace
    Foam::BlockFvm

Description
    Implicitly discretised the divergence of sigma operator including implicit
    treatment of the non-orthogonal correction component of the laplacian
    operator.

    The divergence of sigma contains three distinct operators:
        laplacian;
        laplacian transpose;
        laplacian trace.

    Returns a BlockLduMatrix.

    Note: not templated.

SourceFiles
    BlockFvmDivSigma.C

\*---------------------------------------------------------------------------*/


#ifndef BlockFvmDivSigma_H
#define BlockFvmDivSigma_H

#include "volFieldsFwd.H"
#include "surfaceFieldsFwd.H"
//#include "pointFieldsFwd.H"
#include "blockLduSolvers.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    // Forward declaration
    class solidPolyMesh;

/*---------------------------------------------------------------------------*\
                     Namespace fvm functions Declaration
\*---------------------------------------------------------------------------*/

namespace BlockFvm
{
    // Divergence of sigma stress tensor operator
    tmp<BlockLduMatrix<vector> >
    divSigma
    (
        const solidPolyMesh& solidMesh,
        const surfaceScalarField& mu,
        const surfaceScalarField& lambda,
        GeometricField<vector, fvPatchField, volMesh>& U,
        Field<vector>& blockMSource
    );

    // Laplacian operator with implicit non-orthogonal correction
    tmp<BlockLduMatrix<vector> >
    laplacian
    (
        const solidPolyMesh& solidMesh,
        const surfaceScalarField& mu,
        GeometricField<vector, fvPatchField, volMesh>& U,
        Field<vector>& blockMSource
    );

    // Laplacian transpose operator i.e. div(grad(U).T())
    tmp<BlockLduMatrix<vector> >
    laplacianTranspose
    (
        const solidPolyMesh& solidMesh,
        const surfaceScalarField& mu,
        GeometricField<vector, fvPatchField, volMesh>& U,
        Field<vector>& blockMSource
    );

    // Laplacian trace operator i.e. div(tr(gradU))
    tmp<BlockLduMatrix<vector> >
    laplacianTrace
    (
        const solidPolyMesh& solidMesh,
        const surfaceScalarField& lambda,
        GeometricField<vector, fvPatchField, volMesh>& U,
        Field<vector>& blockMSource
    );
} // End namespace BlockFvm

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Operator not templated yet
//#ifdef NoRepository
//#   include "BlockFvmDivSigma.C"
//#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
