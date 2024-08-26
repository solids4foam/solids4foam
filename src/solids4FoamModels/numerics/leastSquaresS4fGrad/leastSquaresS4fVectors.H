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
    Foam::leastSquaresS4fVectors

Description
    Least-squares gradient scheme vectors

SourceFiles
    leastSquaresS4fVectors.C

\*---------------------------------------------------------------------------*/

#ifndef leastSquaresS4fVectors_H
#define leastSquaresS4fVectors_H

#include "MeshObject.H"
#include "fvMesh.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class leastSquaresS4fVectors Declaration
\*---------------------------------------------------------------------------*/

class leastSquaresS4fVectors
:
    public MeshObject<fvMesh, MoveableMeshObject, leastSquaresS4fVectors>
{
    // Private data

        //- Use boundary face values for patches which do not "fix the value"
        //  i.e., Neumann type boundaries like traction conditions
        //  Defaults to true
        Switch useNeumannBoundaryFaceValues_;

        //- Least-squares gradient vectors
        surfaceVectorField pVectors_;
        surfaceVectorField nVectors_;


    // Private Member Functions

        //- Construct Least-squares gradient vectors
        void calcLeastSquaresVectors();


public:

    // Declare name of the class and its debug switch
    TypeName("leastSquaresS4fVectors");


    // Constructors

        //- Construct given an fvMesh and the useNeumannBoundaryFaceValues flag
        explicit leastSquaresS4fVectors
        (
            const fvMesh&,
            const Switch& useNeumannBoundaryFaceValues = true
        );


    //- Destructor
    virtual ~leastSquaresS4fVectors();


    // Member functions

        //- Return reference to owner least square vectors
        const surfaceVectorField& pVectors() const
        {
            return pVectors_;
        }

        //- Return reference to neighbour least square vectors
        const surfaceVectorField& nVectors() const
        {
            return nVectors_;
        }

        //- Delete the least square vectors when the mesh moves
        virtual bool movePoints();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
