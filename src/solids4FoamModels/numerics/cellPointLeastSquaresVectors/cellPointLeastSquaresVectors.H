/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
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
    Foam::cellPointLeastSquaresVectors

Description
    Least-squares gradient vectors from cell-centres to the points of the cell,
    as used in vertex-centred finite volume schemes.

SourceFiles
    cellPointLeastSquaresVectors.C

Author
    Philip Cardiff, UCD.

\*---------------------------------------------------------------------------*/

#ifndef cellPointLeastSquaresVectors_H
#define cellPointLeastSquaresVectors_H

#include "MeshObject.H"
#include "fvMesh.H"
#include "vectorList.H"
//#include "tolerancesSwitch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class mapPolyMesh;

/*---------------------------------------------------------------------------*\
                Class cellPointLeastSquaresVectors Declaration
\*---------------------------------------------------------------------------*/

class cellPointLeastSquaresVectors
:
#ifdef OPENFOAMESIORFOUNDATION
    public MeshObject<fvMesh, MoveableMeshObject, cellPointLeastSquaresVectors>
#else
    public MeshObject<fvMesh, cellPointLeastSquaresVectors>
#endif
{
    // Private data

        //- Least-squares gradient vectors
        //  There is a vector for each point in each cell
        mutable autoPtr< List<vectorList> > vectorsPtr_;


    // Private member functions

        //- Construct Least-squares gradient vectors
        void makeLeastSquaresVectors() const;


public:

    // Declare name of the class and its debug switch
    TypeName("cellPointLeastSquaresVectors");


    // Constructors

        //- Construct given an fvMesh
        explicit cellPointLeastSquaresVectors(const fvMesh&);


    // Destructor

        virtual ~cellPointLeastSquaresVectors();


    // Member functions

        //- Return reference to cell least square vectors
        const List<vectorList>& vectors() const;

        //- Update after mesh motion:
        //  Delete the least square vectors when the mesh moves
#ifdef OPENFOAMESIORFOUNDATION
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
