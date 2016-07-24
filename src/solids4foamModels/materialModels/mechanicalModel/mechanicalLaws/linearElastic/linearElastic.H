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
    linearElastic

Description
    Hookean linear elastic mechanical law.

SourceFiles
    linearElastic.C

Author
    Philip Cardiff UCD
\*---------------------------------------------------------------------------*/

#ifndef linearElastic_H
#define linearElastic_H

#include "mechanicalLaw.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class linearElastic Declaration
\*---------------------------------------------------------------------------*/

class linearElastic
:
    public mechanicalLaw
{
    // Private data

        //- Density
        dimensionedScalar rho_;

        //- Modulus of elasticity
        dimensionedScalar E_;

        //- Poisson's ratio
        dimensionedScalar nu_;

        //- First Lame parameter
        dimensionedScalar lambda_;

        //- Second Lame parameter (shear modulus)
        dimensionedScalar mu_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        linearElastic(const linearElastic&);

        //- Disallow default bitwise assignment
        void operator=(const linearElastic&);

public:

    //- Runtime type information
    TypeName("linearElastic");

    // Static data members


    // Constructors

        //- Construct from dictionary
        linearElastic
        (
            const word& name,
            const fvMesh& mesh,
            const dictionary& dict
        );


    // Destructor

        virtual ~linearElastic();


    // Member Functions

        //- Return density
        virtual tmp<volScalarField> rho() const;

        //- Return modulus of elasticity
        virtual tmp<volScalarField> E() const;

        //- Return Poisson's ratio
        virtual tmp<volScalarField> nu() const;

        //- Return density
        virtual tmp<volScalarField> rho(scalar t) const
        {
            return rho();
        }

        //- Return modulus of elasticity
        virtual tmp<volScalarField> E(scalar t) const
        {
            return E();
        }

        virtual tmp<volScalarField> E(const volScalarField& epsEq) const
        {
            return E();
        }

        //- Return Poisson's ratio
        virtual tmp<volScalarField> nu(scalar t) const
        {
            return nu();
        }

        //- Return creep compliance
        virtual tmp<volScalarField> J(scalar t) const
        {
            notImplemented(type() + "::J(scalar t)");

            return 1.0/E(t);
        }

        //- Return modulus of plasticity
        virtual tmp<volScalarField> Ep() const;

        virtual tmp<volScalarField> Ep(const volScalarField& sigmaEq) const
        {
            return Ep();
        }

        //- Return yield stress
        virtual tmp<volScalarField> sigmaY() const;

        //- Return yield stress given cellID
        // virtual scalar sigmaY
        // (
        //     const scalar epsilonPEq, const label cellID
        // ) const;

        //- Correct the rheological model
        // virtual void correct()
        // {}

        //- Calculate the stress
        virtual void correct(volSymmTensorField& sigma);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
