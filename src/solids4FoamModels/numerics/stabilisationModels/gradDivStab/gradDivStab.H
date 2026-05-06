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
    gradDivStab

Description
    Face-based stabilisation model based on grad(div(field)).

    This form primarily targets volumetric or compressive oscillations in the
    field while retaining the same cached face-field interface as the other
    stabilisation models.

SourceFiles
    gradDivStab.C
    gradDivStabTemplates.H

Author
    Philip Cardiff, UCD.

\*---------------------------------------------------------------------------*/

#ifndef gradDivStab_H
#define gradDivStab_H

#include "stabilisationModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class gradDivStab Declaration
\*---------------------------------------------------------------------------*/

class gradDivStab
:
    public stabilisationModel
{
    // Private data

        //- Scale factor used for globally scaling the amount of stabilisation
        const scalar scaleFactor_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        gradDivStab
        (
            const gradDivStab&
        );

        //- Disallow default bitwise assignment
        void operator=(const gradDivStab&);

public:

    //- Runtime type information
    TypeName("gradDiv");

    // Constructors

        //- Construct from dictionary
        gradDivStab
        (
            const fvMesh& mesh,
            const dictionary& dict,
            const dimensionSet& dims
        );


    // Destructor

        virtual ~gradDivStab();


    // Member Functions

        // Access

            //- Re-calculate the vector stabilisation
            // Provide the primary variable field and optionally its gradient
            virtual void updateVector
            (
                const volVectorField& field,
                const volTensorField* gradPtr = nullptr
            ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
