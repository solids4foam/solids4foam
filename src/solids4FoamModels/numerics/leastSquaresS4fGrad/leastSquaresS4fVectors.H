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
#include "vectorList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class leastSquaresS4fVectors Declaration
\*---------------------------------------------------------------------------*/

class leastSquaresS4fVectors
:
#ifdef OPENFOAM_NOT_EXTEND
    public MeshObject<fvMesh, MoveableMeshObject, leastSquaresS4fVectors>
#else
    public MeshObject<fvMesh, leastSquaresS4fVectors>
#endif
{
    // Private data

        //- Use boundary face values for patches
        const boolList useBoundaryFaceValues_;

        //- Least-squares gradient vectors
        mutable surfaceVectorField pVectors_;
        mutable surfaceVectorField nVectors_;

        //- Cells whose face stencil is rank-deficient and which instead use
        //  the wider point-cell (vertex-neighbour) stencil
        mutable labelList wideCells_;

        //- For each cell in wideCells_, the cells sharing at least one point
        mutable labelListList wideStencil_;

        //- For each cell in wideCells_, the least-squares vector of each
        //  cell in wideStencil_
        mutable List<vectorList> wideVectors_;


    // Private Member Functions

        //- Construct Least-squares gradient vectors
        void calcLeastSquaresVectors() const;

        //- Replace the face stencil by the wider point-cell stencil in cells
        //  where the face stencil does not span three directions
        //  Takes the assembled (not yet inverted) face-stencil dd tensor
        void calcWideStencilVectors(const symmTensorField& dd) const;


public:

    // Declare name of the class and its debug switch
    TypeName("leastSquaresS4fVectors");


    // Constructors

        //- Construct given an fvMesh and the useNeumannBoundaryFaceValues flag
        explicit leastSquaresS4fVectors
        (
            const fvMesh& mesh,
            const boolList& useBoundaryFaceValues
        );

        //- Construct given the object name, fvMesh and the
        //  useNeumannBoundaryFaceValues flag
        leastSquaresS4fVectors
        (
            const word& objName,
            const fvMesh& mesh,
            const boolList& useBoundaryFaceValues
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

        //- Return a const-access to the boundary face values for patches
        const boolList& useBoundaryFaceValues() const
        {
            return useBoundaryFaceValues_;
        }

        //- Return the cells using the wider point-cell stencil
        const labelList& wideCells() const
        {
            return wideCells_;
        }

        //- Return the point-cell stencil of each cell in wideCells()
        const labelListList& wideStencil() const
        {
            return wideStencil_;
        }

        //- Return the least-squares vectors of each cell in wideCells()
        const List<vectorList>& wideVectors() const
        {
            return wideVectors_;
        }

        //- Delete the least square vectors when the mesh moves
#ifdef OPENFOAM_NOT_EXTEND
        virtual bool movePoints();
#else
        virtual bool movePoints() const;

        //- Update after topo change
        //  Delete the least square vectors when mesh changes
        virtual bool updateMesh(const mapPolyMesh&) const;
#endif
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
