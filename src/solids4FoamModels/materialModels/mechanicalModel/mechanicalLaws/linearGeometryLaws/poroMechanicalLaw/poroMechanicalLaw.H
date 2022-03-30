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
    poroMechanicalLaw

Description
    Generalised poro-coupled mechanical law, where the effective stress is
    defined by any other run-time selectable mechanical law, to which a pore
    stress contribution is added.

    An optional Biot's coefficient can be given to scale the pressure term.

SourceFiles
    poroMechanicalLaw.C

Author
    Philip Cardiff, UCD. All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef poroMechanicalLaw_H
#define poroMechanicalLaw_H

#include "mechanicalLaw.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class poroMechanicalLaw Declaration
\*---------------------------------------------------------------------------*/

class poroMechanicalLaw
:
    public mechanicalLaw
{
    // Private Data

        //- Run-time selectable mechanical law to define the effective stress
        autoPtr<mechanicalLaw> effectiveStressMechLawPtr_;

        //- Biot's coefficient
        //  Optional coefficient to scalar the pressure contribution to the
        //  total stress. Defaults to 1.0.
        const dimensionedScalar b_;

        //- Name of pressure field (defaults to "p")
        const word pName_;

        //- Name of object registry (region) that contains the pressure field
        //  Defaults to checking "solid" then "region0"
        const word pRegion_;

        //- Initial (residual) poro-pressure field
        //  Defaults to zero if not specified
        const volScalarField p0_;

        //- Initial (residual) stress surface field
        mutable surfaceScalarField* p0fPtr_;

    // Private Member Functions

        //- Make the initial pore-pressure surface field
        void makeP0f() const;

        //- Return a reference to the initial pore-pressure surface field
        const surfaceScalarField& p0f() const;

        //- Lookup and return a reference to the pressure field
        const volScalarField& lookupPressureField() const;

        //- Disallow default bitwise copy construct
        poroMechanicalLaw(const poroMechanicalLaw&);

        //- Disallow default bitwise assignment
        void operator=(const poroMechanicalLaw&);

public:

    //- Runtime type information
    TypeName("poroMechanicalLaw");

    // Static data members


    // Constructors

        //- Construct from dictionary
        poroMechanicalLaw
        (
            const word& name,
            const fvMesh& mesh,
            const dictionary& dict,
            const nonLinearGeometry::nonLinearType& nonLinGeom
        );


    // Destructor

        virtual ~poroMechanicalLaw();


    // Member Functions

        //- Return density
        virtual tmp<volScalarField> rho() const;

        //- Return the implicit stiffness
        //  This is the diffusivity for the Laplacian term
        virtual tmp<volScalarField> impK() const;

        //- Calculate the stress
        virtual void correct(volSymmTensorField& sigma);

        //- Calculate the stress
        virtual void correct(surfaceSymmTensorField& sigma);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
