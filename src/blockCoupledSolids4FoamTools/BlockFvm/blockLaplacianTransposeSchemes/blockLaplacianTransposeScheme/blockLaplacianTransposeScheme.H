/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     |
    \\  /    A nd           | For copyright notice see file Copyright
     \\/     M anipulation  |
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
    Foam::fv::blockLaplacianTranspose

Description
    Abstract base class for block coupled of lapacian transpose schemes.

SourceFiles
    blockLaplacianTranspose.C

Author
    Philip Cardiff UCD

\*---------------------------------------------------------------------------*/

#ifndef blockLaplacianTranspose_H
#define blockLaplacianTranspose_H

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
                           Class blockLaplacianTranspose Declaration
\*---------------------------------------------------------------------------*/

class blockLaplacianTranspose
:
    public refCount
{

protected:

    // Protected data

        const fvMesh& mesh_;

    // Private Member Functions

        //- Disallow copy construct
        blockLaplacianTranspose(const blockLaplacianTranspose&);

        //- Disallow default bitwise assignment
        void operator=(const blockLaplacianTranspose&);


public:

    // //- Runtime type information
    virtual const word& type() const = 0;

    // Declare run-time constructor selection tables

        declareRunTimeSelectionTable
        (
            tmp,
            blockLaplacianTranspose,
            Istream,
            (const fvMesh& mesh, Istream& schemeData),
            (mesh, schemeData)
        );


    // Constructors

        //- Construct from mesh and Istream
        blockLaplacianTranspose(const fvMesh& mesh, Istream& is)
        :
            mesh_(mesh)
        {}


    // Selectors

        //- Return a pointer to a new blockLaplacianTranspose created on
        //  freestore
        static tmp<blockLaplacianTranspose> New
        (
            const fvMesh& mesh,
            Istream& schemeData
        );


    // Destructor

        virtual ~blockLaplacianTranspose();


    // Member Functions

        //- Return mesh reference
        const fvMesh& mesh() const
        {
            return mesh_;
        }

        //- Implicit laplacian transpose i.e. div(gradU.T())
        virtual tmp<BlockLduMatrix<vector> > fvmBlockLaplacianTranspose
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
