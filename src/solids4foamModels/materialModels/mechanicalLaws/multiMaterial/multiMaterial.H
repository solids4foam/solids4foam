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
    multiMaterial

Description
    Zoned multi-material mechanical controlled by an indicator field.

SourceFiles
    multiMaterial.C

\*---------------------------------------------------------------------------*/

#ifndef multiMaterial_H
#define multiMaterial_H

#include "mechanicalLaw.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class multiMaterial Declaration
\*---------------------------------------------------------------------------*/

class multiMaterial
:
    public mechanicalLaw,
    public PtrList<mechanicalLaw>
{
    // Private data

        //- Material indicator field
        volScalarField materials_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        multiMaterial(const multiMaterial&);

        //- Disallow default bitwise assignment
        void operator=(const multiMaterial&);


        //- Calculate indicator field given index
        tmp<scalarField> indicator(const label index) const;

        //- Calculate indicator field given index
        scalar indicator(const label index, const label cellID) const;

public:

    //- Runtime type information
    TypeName("multiMaterial");

    // Static data members


    // Constructors

        //- Construct from dictionary
        multiMaterial
        (
            const word& name,
            const fvMesh& mesh,
            const dictionary& dict
        );


    // Destructor

        virtual ~multiMaterial();


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

        //- Return material residual i.e. a measured of how convergence of
        //  the material model
        virtual scalar residual();

        //- Update the yield stress, called at end of time-step
        virtual void updateYieldStress();

        //- Rate of dissipated energy from plasticity
        //virtual tmp<volScalarField> plasticDissipationRate() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
