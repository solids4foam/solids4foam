/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
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
    newSubsetMotionSolverFvMesh

Description

SourceFiles
    newSubsetMotionSolverFvMesh.C

\*---------------------------------------------------------------------------*/

#ifndef newSubsetMotionSolverFvMesh_H
#define newSubsetMotionSolverFvMesh_H

#include "dynamicFvMesh.H"
#include "newFvMeshSubset.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class motionSolver;

/*---------------------------------------------------------------------------*\
                   Class newSubsetMotionSolverFvMesh Declaration
\*---------------------------------------------------------------------------*/

class newSubsetMotionSolverFvMesh
:
    public dynamicFvMesh
{
    // Private data

        //- Moving mesh dictionary
        dictionary movingMeshCoeffs_;

        //- Subset mesh
        newFvMeshSubset subsetMesh_;

        //- Pointer to mesh motion
        autoPtr<motionSolver> motionPtr_;

        //- Motion under-relaxation coefficient
        scalar alpha_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        newSubsetMotionSolverFvMesh(const newSubsetMotionSolverFvMesh&);

        //- Disallow default bitwise assignment
        void operator=(const newSubsetMotionSolverFvMesh&);


public:

    //- Runtime type information
    TypeName("newSubsetMotionSolverFvMesh");


    // Constructors

        //- Construct from IOobject
        newSubsetMotionSolverFvMesh(const IOobject& io);


    // Destructor

        virtual ~newSubsetMotionSolverFvMesh();


    // Member Functions

        //- Return reference to subset mesh
        const newFvMeshSubset& subsetMesh() const;

        //- Update the mesh for both mesh motion and topology change
        virtual bool update();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
