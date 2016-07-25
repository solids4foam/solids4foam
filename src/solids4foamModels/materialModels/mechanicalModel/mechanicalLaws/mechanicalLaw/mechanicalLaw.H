/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     |
    \\  /    A nd           | For copyright notice see file Copyright
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Class
    mechanicalLaw

Description
    Material mechanical for solids.

SourceFiles
    mechanicalLaw.C
    newmechanicalLaw.C

\*---------------------------------------------------------------------------*/

#ifndef mechanicalLaw_H
#define mechanicalLaw_H

#include "IOdictionary.H"
#include "typeInfo.H"
#include "runTimeSelectionTables.H"
#include "volFields.H"
#include "tmp.H"
#include "autoPtr.H"
//#include "plasticityStressReturn.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class mechanicalLaw Declaration
\*---------------------------------------------------------------------------*/

class mechanicalLaw
{
    // Private data

        //- Name
        const word name_;

        //- Const reference to mesh
        const fvMesh& mesh_;

        //- Const reference to law dictionary
        const dictionary& dict_;

        // Current materials index field
        //     1 for the current materials
        //     0 for other materials
        mutable volScalarField* curMaterialPtr_;

        // Current material index i.e. is the current material number 0, 1,
        // 2, etc.
        mutable label curMatIndex_;


    // Private Member Functions

        // Make current material field
        void calcCurMaterial() const;

        //- Disallow copy construct
        mechanicalLaw(const mechanicalLaw&);

        //- Disallow default bitwise assignment
        void operator=(const mechanicalLaw&);


protected:

        //- Return reference to mesh
        const fvMesh& mesh() const
        {
            return mesh_;
        }


public:

    //- Runtime type information
    TypeName("mechanicalLaw");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            mechanicalLaw,
            dictionary,
            (
                const word name,
                const fvMesh& mesh,
                const dictionary& dict
            ),
            (name, mesh, dict)
        );


    // Selectors

        //- Return a reference to the selected mechanical model
        static autoPtr<mechanicalLaw> New
        (
            const word& name,
            const fvMesh& mesh,
            const dictionary& dict
        );


    // Constructors

        //- Construct from dictionary
        mechanicalLaw
        (
            const word& name,
            const fvMesh& mesh,
            const dictionary& dict
        );


    // Destructor

        virtual ~mechanicalLaw()
        {}


    // Member Functions

        //- Return name
        const word& name() const
        {
            return name_;
        }

        //- Return density
        virtual tmp<volScalarField> rho() const = 0;

        //- Return the implicit stiffness
        //  This is the diffusivity for the Laplacian term
        virtual tmp<volScalarField> impK() const = 0;

        //- Return the implicit stiffness surface field
        //  This is the diffusivity for the Laplacian term
        virtual tmp<surfaceScalarField> impKf() const;

        //- Calculate the stress
        virtual void correct(volSymmTensorField& sigma) = 0;

        //- Update stress based on plasticity stress return
        //virtual void correct(surfaceSymmTensorField& sigma) = 0;

        //- Return material residual i.e. a measured of how convergence of
        //  the material model
        virtual scalar residual();

        //- Update the yield stress
        //- called at end of time-step
        virtual void updateYieldStress()
        {}

        // Return reference to current material field
        const volScalarField& curMaterial() const;

        // Return reference to current material field
        volScalarField& curMaterial();

        // Return reference to current material index
        label& curMaterialIndex() const
        {
            return curMatIndex_;
        }

        //- Set current material index; this is only needed for multiMaterial
        //  analyses
        virtual void setMaterialIndex(label curMaterialIndex)
        {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
