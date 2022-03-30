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
    anisotropicBiotElastic

Description
    Anisotropic elastic mechanical law ported from Tian Tang's minigeotechfoam
    at https://bitbucket.org/tiantang/minigeotechfoam.

    Note: to include pore pressure, you can use this law with the general
    poroMechanicalLaw wrapper mechanical law.

    More reading at:

    T.Tang, O.Hededal & P.Cardiff (2014). On finite volume method implementation
    of poro-elasto-plasticity soil model. International Journal for
    Numerical and Analytical Methods in Geomechanics, 10.1002/nag.2361.

    T.Tang & O.Hededal (2014). Simulation of pore pressure accumulation under
    cyclic loading using finite volume method. Proceedings of the 8th European
    Conference on Numerical Methods in Geotechnical Engineering (NUMGE14),
    Volume 2, Pages 1301-1306.

SourceFiles
    anisotropicBiotElastic.C

Author
    Original author: Tian Tang, DTU. All rights reserved.
    Ported to solids4foam by: Philip Cardiff, UCD. All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef anisotropicBiotElastic_H
#define anisotropicBiotElastic_H

#include "mechanicalLaw.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class anisotropicBiotElastic Declaration
\*---------------------------------------------------------------------------*/

class anisotropicBiotElastic
:
    public mechanicalLaw
{
    // Private Data

        //- Flag to indicate a 2-D model
        const bool model2d_;

        //- Density
        const dimensionedScalar rho_;

        //- Elastic stiffness parameters
        scalar A11_;
        scalar A22_;
        scalar A33_;
        scalar A44_;
        scalar A55_;
        scalar A66_;
        scalar A12_;
        scalar A21_;
        scalar A31_;
        scalar A23_;

        //- Total strain
        volSymmTensorField epsilon_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        anisotropicBiotElastic(const anisotropicBiotElastic&);

        //- Disallow default bitwise assignment
        void operator=(const anisotropicBiotElastic&);

public:

    //- Runtime type information
    TypeName("anisotropicBiotElastic");

    // Static data members


    // Constructors

        //- Construct from dictionary
        anisotropicBiotElastic
        (
            const word& name,
            const fvMesh& mesh,
            const dictionary& dict,
            const nonLinearGeometry::nonLinearType& nonLinGeom
        );


    // Destructor

        virtual ~anisotropicBiotElastic();


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
