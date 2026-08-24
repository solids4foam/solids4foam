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
    laplacianStab

Description
    Face-based stabilisation model based on the Laplacian of the primary
    field.

    This is the simplest Laplacian-based stabilisation in the directory and
    provides the common cached face-field and approximate Jacobian interface.

SourceFiles
    laplacianStab.C
    laplacianStabTemplates.H

Author
    Philip Cardiff, UCD.

\*---------------------------------------------------------------------------*/

#ifndef laplacianStab_H
#define laplacianStab_H

#include "stabilisationModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class laplacianStab Declaration
\*---------------------------------------------------------------------------*/

class laplacianStab
:
    public stabilisationModel
{
    // Private Member Functions

        //- Template calculation of the stabilisation
        template<class GeoField, class SurfaceField>
        void computeDiffStencil
        (
            const GeoField& field,
            SurfaceField& faceFlux,
            scalar scaleFactor
        ) const;

        //- Disallow default bitwise copy construct
        laplacianStab
        (
            const laplacianStab&
        );

        //- Disallow default bitwise assignment
        void operator=(const laplacianStab&);

public:

    //- Runtime type information
    TypeName("laplacian");

    // Constructors

        //- Construct from dictionary
        laplacianStab
        (
            const fvMesh& mesh,
            const dictionary& dict,
            const dimensionSet& dims
        );


    // Destructor

        virtual ~laplacianStab();


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

            //- Provide an implicit approximation of the stabilisation
            //  This is an approximate Jacobian of the stabilisation term
            //  Provide the primary variable field and optionally a surface
            //  diffusivity field can be optionally provided
            virtual const fvScalarMatrix& scalarJacobian
            (
                const volScalarField& field,
                const surfaceScalarField* gammaPtr = nullptr,
                const bool rebuild = false
            ) const;

            //- Provide an implicit approximation of the stabilisation
            //  This is an approximate Jacobian of the stabilisation term
            //  Provide the primary variable field and optionally a surface
            //  diffusivity field can be optionally provided
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
#   include "laplacianStabTemplates.H"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
