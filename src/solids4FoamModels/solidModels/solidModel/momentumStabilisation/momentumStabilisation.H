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
    momentumStabilisation

Description
    Run-time selectable stabilisation (diffusion) term for the solid mechanics
    momentum equation.

Author
    Philip Cardiff, UCD.

SourceFiles
    momentumStabilisation.C

\*---------------------------------------------------------------------------*/

#ifndef momentumStabilisation_H
#define momentumStabilisation_H

#include "dictionary.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class momentumStabilisation Declaration
\*---------------------------------------------------------------------------*/

class momentumStabilisation
{
    // Private data

        //- Copy of settings dict
        const dictionary dict_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        momentumStabilisation(const momentumStabilisation&);

        //- Disallow default bitwise assignment
        void operator=(const momentumStabilisation&);

public:

    //- Runtime type information
    TypeName("momentumStabilisation");


    // Constructors

        //- Construct from dictionary
        momentumStabilisation
        (
            const dictionary& dict
        );


    // Destructor

        virtual ~momentumStabilisation()
        {}


    // Member Functions

        // Access

            //- Const access to the dict
            const dictionary& dict() const
            {
                return dict_;
            }

        // Evaluation

            //- Stabilisation term
            tmp<volVectorField> stabilisation
            (
                const volVectorField& vf,
                const volTensorField& gradVf,
                const volScalarField& gamma
            ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
