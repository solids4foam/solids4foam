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
    poroLinearElastic

Description
    Poro-elastic mechanical law where the total stress is the sum of the solid
    skeleton stress, governed by Hooke's law, and the poro-pressure.

    More reading at:

    T.Tang, O.Hededal & P.Cardiff (2014). On finite volume method implementation
    of poro-elasto-plasticity soil model. International Journal for
    Numerical and Analytical Methods in Geomechanics, 10.1002/nag.2361.

    T.Tang & O.Hededal (2014). Simulation of pore pressure accumulation under
    cyclic loading using finite volume method. Proceedings of the 8th European
    Conference on Numerical Methods in Geotechnical Engineering (NUMGE14),
    Volume 2, Pages 1301-1306.

SourceFiles
    poroLinearElastic.C

Author
    Tian Tang, DTU. All rights reserved.
    Philip Cardiff, UCD. All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef poroLinearElastic_H
#define poroLinearElastic_H

#include "mechanicalLaw.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class poroLinearElastic Declaration
\*---------------------------------------------------------------------------*/

class poroLinearElastic
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
        poroLinearElastic(const poroLinearElastic&);

        //- Disallow default bitwise assignment
        void operator=(const poroLinearElastic&);

public:

    //- Runtime type information
    TypeName("poroLinearElastic");

    // Static data members


    // Constructors

        //- Construct from dictionary
        poroLinearElastic
        (
            const word& name,
            const fvMesh& mesh,
            const dictionary& dict
        );


    // Destructor

        virtual ~poroLinearElastic();


    // Member Functions

        //- Return density
        virtual tmp<volScalarField> rho() const;

        //- Return the implicit stiffness
        //  This is the diffusivity for the Laplacian term
        virtual tmp<volScalarField> impK() const;

        //- Provide access to mu for the coupledUnsNonLinLinearElasticSolid
        virtual const dimensionedScalar& mu() const;

        //- Provide access to lambda for the coupledUnsNonLinLinearElasticSolid
        virtual const dimensionedScalar& lambda() const;

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
