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
    generalisedEvenOrderLaplacianStab

Description
    Face-based stabilisation model based on repeated Laplacian application.

    The model generalises the Laplacian-based stabilisation family to higher
    even orders through the user-selected Laplacian power.

SourceFiles
    generalisedEvenOrderLaplacianStab.C
    generalisedEvenOrderLaplacianStabTemplates.H

Author
    Philip Cardiff, UCD. All rights reserved

\*---------------------------------------------------------------------------*/

#ifndef generalisedEvenOrderLaplacianStab_H
#define generalisedEvenOrderLaplacianStab_H

#include "stabilisationModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class generalisedEvenOrderLaplacianStab Declaration
\*---------------------------------------------------------------------------*/

class generalisedEvenOrderLaplacianStab
:
    public stabilisationModel
{
    // Private data

        //- Scale factor used for globally scaling the amount of stabilisation
        const scalar scaleFactor_;

        //- Laplacian power, e.g. 1 == laplacianStab, 2 == JSTStab, ..
        const label laplacianPower_;   // m ≥ 1


    // Private Member Functions

        //- Template calculation of the stabilisation
        template<class GeoField, class SurfaceField>
        void computeDiffStencil
        (
            const GeoField& field,
            SurfaceField& faceFlux,
            scalar scaleFactor,
            const label laplacianPower   // m ≥ 1
        ) const;

        //- Disallow default bitwise copy construct
        generalisedEvenOrderLaplacianStab
        (
            const generalisedEvenOrderLaplacianStab&
        );

        //- Disallow default bitwise assignment
        void operator=(const generalisedEvenOrderLaplacianStab&);

public:

    //- Runtime type information
    TypeName("generalisedEvenOrderLaplacian");

    // Constructors

        //- Construct from dictionary
        generalisedEvenOrderLaplacianStab
        (
            const fvMesh& mesh,
            const dictionary& dict,
            const dimensionSet& dims
        );


    // Destructor

        virtual ~generalisedEvenOrderLaplacianStab();


    // Member Functions

        // Access

            //- Re-calculate the scalar stabilisation
            //  Provide the primary variable field and optionally its gradient
            virtual void updateScalar
            (
                const volScalarField& field,
                const volVectorField* gradPtr = nullptr
            ) const;

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

#ifdef NoRepository
#   include "generalisedEvenOrderLaplacianStabTemplates.H"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
