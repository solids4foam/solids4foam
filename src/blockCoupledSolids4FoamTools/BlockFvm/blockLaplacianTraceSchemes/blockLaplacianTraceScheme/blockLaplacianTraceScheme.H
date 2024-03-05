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
    Foam::fv::blockLaplacianTrace

Description
    Abstract base class for block coupled of laplacian trace schemes.

SourceFiles
    blockLaplacianTrace.C

Author
    Philip Cardiff UCD

\*---------------------------------------------------------------------------*/

#ifndef blockLaplacianTrace_H
#define blockLaplacianTrace_H

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
                           Class blockLaplacianTrace Declaration
\*---------------------------------------------------------------------------*/

class blockLaplacianTrace
:
    public refCount
{

protected:

    // Protected data

        const fvMesh& mesh_;

    // Private Member Functions

        //- Disallow copy construct
        blockLaplacianTrace(const blockLaplacianTrace&);

        //- Disallow default bitwise assignment
        void operator=(const blockLaplacianTrace&);


public:

    // //- Runtime type information
    virtual const word& type() const = 0;

    // Declare run-time constructor selection tables

        declareRunTimeSelectionTable
        (
            tmp,
            blockLaplacianTrace,
            Istream,
            (const fvMesh& mesh, Istream& schemeData),
            (mesh, schemeData)
        );


    // Constructors

        //- Construct from mesh and Istream
        blockLaplacianTrace(const fvMesh& mesh, Istream& is)
        :
            mesh_(mesh)
        {}


    // Selectors

        //- Return a pointer to a new blockLaplacianTrace created on
        //  freestore
        static tmp<blockLaplacianTrace> New
        (
            const fvMesh& mesh,
            Istream& schemeData
        );


    // Destructor

        virtual ~blockLaplacianTrace();


    // Member Functions

        //- Return mesh reference
        const fvMesh& mesh() const
        {
            return mesh_;
        }

        //- Implicit laplacian trace i.e. div(lambda*I*tr(gradU))
        virtual tmp<BlockLduMatrix<vector> > fvmBlockLaplacianTrace
        (
            const solidPolyMesh& solidMesh,
            const surfaceScalarField& lambdaf,
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
