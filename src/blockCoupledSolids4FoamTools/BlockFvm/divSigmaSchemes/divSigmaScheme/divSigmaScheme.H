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
    Foam::fv::divSigmaScheme

Description
    Abstract base class for divergence of sigma schemes.

SourceFiles
    divSigmaScheme.C

Author
    Philip Cardiff UCD

\*---------------------------------------------------------------------------*/

#ifndef divSigmaScheme_H
#define divSigmaScheme_H

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
                           Class divSigmaScheme Declaration
\*---------------------------------------------------------------------------*/

class divSigmaScheme
:
    public refCount
{

protected:

    // Protected data

        const fvMesh& mesh_;

    // Private Member Functions

        //- Disallow copy construct
        divSigmaScheme(const divSigmaScheme&);

        //- Disallow default bitwise assignment
        void operator=(const divSigmaScheme&);


public:

    // //- Runtime type information
    virtual const word& type() const = 0;

    // Declare run-time constructor selection tables

        declareRunTimeSelectionTable
        (
            tmp,
            divSigmaScheme,
            Istream,
            (const fvMesh& mesh, Istream& schemeData),
            (mesh, schemeData)
        );


    // Constructors

        //- Construct from mesh and Istream
        divSigmaScheme(const fvMesh& mesh, Istream& is)
        :
            mesh_(mesh)
        {}


    // Selectors

        //- Return a pointer to a new divSigmaScheme created on freestore
        static tmp<divSigmaScheme> New
        (
            const fvMesh& mesh,
            Istream& schemeData
        );


    // Destructor

        virtual ~divSigmaScheme();


    // Member Functions

        //- Return mesh reference
        const fvMesh& mesh() const
        {
            return mesh_;
        }

        virtual tmp<BlockLduMatrix<vector> > fvmDivSigma
        (
            const solidPolyMesh& solidMesh,
            const surfaceScalarField& muf,
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
