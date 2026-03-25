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
    alphaStab

Description
    Face-based stabilisation model based on Nishikawa's alpha scheme.

    The stabilisation is formed from an extrapolated jump across each face and
    is intended to damp grid-scale oscillations while remaining consistent for
    smooth solutions.

SourceFiles
    alphaStab.C
    alphaStabTemplates.H

Author
    Philip Cardiff, UCD. All rights reserved

\*---------------------------------------------------------------------------*/

#ifndef alphaStab_H
#define alphaStab_H

#include "stabilisationModel.H"
#include "lookupSolidModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class alphaStab Declaration
\*---------------------------------------------------------------------------*/

class alphaStab
:
    public stabilisationModel
{
    // Private data

        //- Scale factor used for globally scaling the amount of stabilisation
        const scalar scaleFactor_;


    // Private Member Functions

        //- Template calculation of the stabilisation
        template<class GeoField, class GradField, class SurfaceField>
        void computeDiffStencil
        (
            const GeoField& field,
            const GradField& gradField,
            SurfaceField& faceFlux,
            scalar scaleFactor
        ) const;

        void computeDiffStencilHighOrderScalar
        (
            const volScalarField& field,
            const volVectorField& gradField,
            surfaceScalarField& faceFlux,
            scalar scaleFactor
        ) const;

        void computeDiffStencilHighOrderVector
        (
            const volVectorField& field,
            const volTensorField& gradField,
            surfaceVectorField& faceFlux,
            scalar scaleFactor
        ) const;

        //- Disallow default bitwise copy construct
        alphaStab(const alphaStab&);

        //- Disallow default bitwise assignment
        void operator=(const alphaStab&);

public:

    //- Runtime type information
    TypeName("alpha");

    // Constructors

        //- Construct from dictionary
        alphaStab
        (
            const fvMesh& mesh,
            const dictionary& dict,
            const dimensionSet& dims
        );


    // Destructor

        virtual ~alphaStab();


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

            //- Provide an approximate implicit Jacobian contribution
            virtual const fvVectorMatrix& vectorJacobian
            (
                const volVectorField& field,
                const surfaceScalarField* gammaPtr = nullptr,
                const bool rebuild = false
            ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#ifdef NoRepository
#   include "alphaStabTemplates.H"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
