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
    Foam::stabLeastSquaresVectors

Description
    Least-squares gradient scheme vectors

SourceFiles
    stabLeastSquaresVectors.C

\*---------------------------------------------------------------------------*/

#ifndef stabLeastSquaresVectors_H
#define stabLeastSquaresVectors_H

#include "MeshObject.H"
#include "fvMesh.H"
#include "surfaceFieldsFwd.H"
#include "labelPair.H"
#include "tolerancesSwitch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class mapPolyMesh;

/*---------------------------------------------------------------------------*\
                     Class stabLeastSquaresVectors Declaration
\*---------------------------------------------------------------------------*/

class stabLeastSquaresVectors
:
    public MeshObject<fvMesh, stabLeastSquaresVectors>
{
    // Private data

        //- Least-squares gradient vectors
        mutable surfaceVectorField* pVectorsPtr_;
        mutable surfaceVectorField* nVectorsPtr_;


    // Private member functions

        //- Construct Least-squares gradient vectors
        void makeLeastSquaresVectors() const;


public:

    // Declare name of the class and its debug switch
    TypeName("stabLeastSquaresVectors");


    // Static data members

        //- Point merge tolerance
        static const debug::tolerancesSwitch smallDotProdTol_;


    // Constructors

        //- Construct given an fvMesh
        explicit stabLeastSquaresVectors(const fvMesh&);


    // Destructor

        virtual ~stabLeastSquaresVectors();


    // Member functions

        //- Return reference to owner least square vectors
        const surfaceVectorField& pVectors() const;

        //- Return reference to neighbour least square vectors
        const surfaceVectorField& nVectors() const;


        //- Update after mesh motion:
        //  Delete the least square vectors when the mesh moves
        virtual bool movePoints() const;

        //- Update after topo change:
        //  Delete the least square vectors when mesh changes
        virtual bool updateMesh(const mapPolyMesh&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
