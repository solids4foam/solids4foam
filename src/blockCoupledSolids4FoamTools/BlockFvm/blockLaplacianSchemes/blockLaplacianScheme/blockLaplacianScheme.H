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
    Foam::fv::blockLaplacian

Description
    Abstract base class for block coupled of lapacian schemes.

SourceFiles
    blockLaplacian.C

Author
    Philip Cardiff UCD

\*---------------------------------------------------------------------------*/

#ifndef blockLaplacian_H
#define blockLaplacian_H

#include "tmp.H"
#include "volFieldsFwd.H"
#include "surfaceFieldsFwd.H"
#include "linear.H"
#include "correctedSnGrad.H"
#include "typeInfo.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class Type>
class fvMatrix;

class fvMesh;
class solidPolyMesh;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

/*---------------------------------------------------------------------------*\
                           Class blockLaplacian Declaration
\*---------------------------------------------------------------------------*/

class blockLaplacian
:
    public refCount
{

protected:

    // Protected data

        const fvMesh& mesh_;

    // Private Member Functions

        //- Disallow copy construct
        blockLaplacian(const blockLaplacian&);

        //- Disallow default bitwise assignment
        void operator=(const blockLaplacian&);


public:

    // //- Runtime type information
    virtual const word& type() const = 0;

    // Declare run-time constructor selection tables

        declareRunTimeSelectionTable
        (
            tmp,
            blockLaplacian,
            Istream,
            (const fvMesh& mesh, Istream& schemeData),
            (mesh, schemeData)
        );


    // Constructors

        //- Construct from mesh and Istream
        blockLaplacian(const fvMesh& mesh, Istream& is)
        :
            mesh_(mesh)
        {}


    // Selectors

        //- Return a pointer to a new blockLaplacian created on freestore
        static tmp<blockLaplacian> New
        (
            const fvMesh& mesh,
            Istream& schemeData
        );


    // Destructor

        virtual ~blockLaplacian();


    // Member Functions

        //- Return mesh reference
        const fvMesh& mesh() const
        {
            return mesh_;
        }

        //- Implicit laplacian with implicit non-orthgonal correction
        virtual tmp<BlockLduMatrix<vector> > fvmBlockLaplacian
        (
            const solidPolyMesh& solidMesh,
            const surfaceScalarField& muf,
            GeometricField<vector, fvPatchField, volMesh>& U,
            Field<vector>& blockB
        ) = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
