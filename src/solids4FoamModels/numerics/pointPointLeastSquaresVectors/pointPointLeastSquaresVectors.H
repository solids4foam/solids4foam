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
    Foam::pointPointLeastSquaresVectors

Description
    Least-squares gradient vectors from points to the point-points neighbours,
    as used in vertex-centred finite volume schemes.

SourceFiles
    pointPointLeastSquaresVectors.C

Author
    Philip Cardiff, UCD.

\*---------------------------------------------------------------------------*/

#ifndef pointPointLeastSquaresVectors_H
#define pointPointLeastSquaresVectors_H

#include "MeshObject.H"
#include "fvMesh.H"
#include "vectorList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class mapPolyMesh;

/*---------------------------------------------------------------------------*\
                Class pointPointLeastSquaresVectors Declaration
\*---------------------------------------------------------------------------*/

class pointPointLeastSquaresVectors
:
#ifdef OPENFOAM_NOT_EXTEND
    public MeshObject<fvMesh, MoveableMeshObject, pointPointLeastSquaresVectors>
#else
    public MeshObject<fvMesh, pointPointLeastSquaresVectors>
#endif
{
    // Private data

        //- Least-squares gradient vectors
        //  There is a vector for each point-neighbour of aa point
        mutable autoPtr< List<vectorList> > vectorsPtr_;


    // Private member functions

        //- Construct Least-squares gradient vectors
        void makeLeastSquaresVectors() const;


public:

    // Declare name of the class and its debug switch
    TypeName("pointPointLeastSquaresVectors");


    // Constructors

        //- Construct given an fvMesh
        explicit pointPointLeastSquaresVectors(const fvMesh&);


    // Destructor

        virtual ~pointPointLeastSquaresVectors();


    // Member functions

        //- Return reference to point least square vectors
        const List<vectorList>& vectors() const;

        //- Update after mesh motion:
        //  Delete the least square vectors when the mesh moves
#ifdef OPENFOAM_NOT_EXTEND
        virtual bool movePoints();
#else
        virtual bool movePoints() const;
#endif

        //- Update after topo change:
        //  Delete the least square vectors when mesh changes
        virtual bool updateMesh(const mapPolyMesh&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
