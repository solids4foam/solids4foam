/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     3.2
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
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
    robinInterfaceHelper

Description
    Heloer class for apply Robin boundary conditions for the fluid-solid
    interaction simulations. 

Author
    Philip Cardiff, UCD.

SourceFiles
    robinInterfaceHelper.C

\*---------------------------------------------------------------------------*/

#ifndef robinInterfaceHelper_H
#define robinInterfaceHelper_H

#include "dictionary.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class robinInterfaceHelper Declaration
\*---------------------------------------------------------------------------*/

class robinInterfaceHelper
{
    // Private data

        //- List of patches which have a Robin condition
        labelList patchIDs_;

        //- Velocity ddt scheme
        word ddtScheme_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        robinInterfaceHelper(const robinInterfaceHelper&);

        //- Disallow default bitwise assignment
        void operator=(const robinInterfaceHelper&);

public:

    //- Runtime type information
    TypeName("robinInterfaceHelper");


    // Constructors

        //- Construct from the velocity and pressure fields
        robinInterfaceHelper
        (
            const volVectorField& U,
            const volScalarField& p
        );


    // Destructor

        virtual ~robinInterfaceHelper()
        {}


    // Member Functions

        // Access

            //- Const access to patchIDs
            const labelList& patchIDs() const
            {
                return patchIDs_;
            }


        // Evaluation

            //- Reset flux if mesh changed
            void setInterfaceFluxToZero
            (
                surfaceScalarField& phi
            ) const;

            //- Update phi for Robin FSI boundary conditions
            void updateInterface
            (
                surfaceScalarField& phi,
                surfaceScalarField& rAU
            ) const;

            //- Correct phi after solving pressure equation
            void correctInterface
            (
                const volScalarField& p,
                surfaceScalarField& phi,
                surfaceScalarField& rAU
            ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
