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
    Foam::leastSquaresReconstruction

Description
    Mesh object that owns field-specific least-squares reconstruction schemes.

SourceFiles
    leastSquaresReconstruction.C

Author
    Ivan Batistic, UCD.
    Philip Cardiff, UCD.

\*---------------------------------------------------------------------------*/

#ifndef leastSquaresReconstruction_H
#define leastSquaresReconstruction_H

#include "MeshObject.H"
#include "leastSquaresScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
              Class leastSquaresReconstruction Declaration
\*---------------------------------------------------------------------------*/

class leastSquaresReconstruction
:
#ifdef FOAMEXTEND
    public MeshObject<fvMesh, leastSquaresReconstruction>
#else
    public MeshObject
    <
        fvMesh,
        UpdateableMeshObject,
        leastSquaresReconstruction
    >
#endif
{
    // Private data

        //- Displacement reconstruction scheme pointer
        mutable autoPtr<leastSquaresScheme> displacementSchemePtr_;

        //- Pressure reconstruction scheme pointer
        mutable autoPtr<leastSquaresScheme> pressureSchemePtr_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        leastSquaresReconstruction(const leastSquaresReconstruction&);

        //- Disallow default bitwise assignment
        void operator=(const leastSquaresReconstruction&);


public:

    //- Runtime type information
    TypeName("leastSquaresReconstruction");


    // Constructors

        //- Construct from mesh
        explicit leastSquaresReconstruction(const fvMesh& mesh);


    //- Destructor
    virtual ~leastSquaresReconstruction();


    // Member Functions

        //- Return the reconstruction scheme for the given field role
        const leastSquaresScheme& scheme
        (
            const word& fieldRole,
            const boolList& includePatchInStencils,
            const dictionary& dict
        ) const;

        //- Clear geometry-dependent data in constructed schemes
        void clear() const;

#ifdef FOAMEXTEND
        //- Clear geometry-dependent data after mesh motion
        virtual bool movePoints() const;

        //- Delete constructed schemes after a topology change
        virtual bool updateMesh(const mapPolyMesh& mpm) const;

        //- Delete from the registry without dispatching through an fvMesh
        //  whose derived lifetime has ended during polyMesh destruction
        virtual bool deleteObject() const;
#else
        //- Clear geometry-dependent data after mesh motion
        virtual bool movePoints();

        //- Delete constructed schemes after a topology change
        virtual void updateMesh(const mapPolyMesh& mpm);
#endif
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
