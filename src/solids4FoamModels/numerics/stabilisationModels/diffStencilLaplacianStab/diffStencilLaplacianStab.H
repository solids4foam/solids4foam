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
    diffStencilLaplacianStab

Description
    Face-based stabilisation model based on the difference between compact-
    stencil and extended-stencil gradient/Laplacian approximations.

    The stabilisation is formed from the mismatch between the two operators
    and is intended to damp grid-scale oscillations while remaining
    consistent for smooth solutions.

    On uniform Cartesian meshes with a centred cell gradient and linear
    interpolation, its spectral shape is proportional to sum_q a_q^2 and
    coincides with the directionally split m=2 spectral shape. The generalised
    directionally split family is not implemented here.

SourceFiles
    diffStencilLaplacianStab.C
    diffStencilLaplacianStabTemplates.H

Author
    Philip Cardiff, UCD.

\*---------------------------------------------------------------------------*/

#ifndef diffStencilLaplacianStab_H
#define diffStencilLaplacianStab_H

#include "stabilisationModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class diffStencilLaplacianStab Declaration
\*---------------------------------------------------------------------------*/

class diffStencilLaplacianStab
:
    public stabilisationModel
{
    // Private data

        //- Scale factor used for globally scaling the amount of stabilisation
        const scalar scaleFactor_;

        //- Scale factor used for globally scaling the approximate Jacobian
        //  This defaults to 1.0
        const scalar scaleFactorJacobian_;


    // Private Member Functions

        //- Return the operator-specific reference-mode normalisation
        //  1/(4*rStar) for the raw diffStencil/Rhie-Chow response 4*rStar
        scalar diffStencilNormalisation() const;

        //- Template calculation of the stabilisation
        template<class GeoField, class GradField, class SurfaceField>
        void computeDiffStencil
        (
            const GeoField& field,
            const GradField& gradField,
            SurfaceField& faceFlux,
            scalar scaleFactor
        ) const;

        //- Disallow default bitwise copy construct
        diffStencilLaplacianStab(const diffStencilLaplacianStab&);

        //- Disallow default bitwise assignment
        void operator=(const diffStencilLaplacianStab&);

public:

    //- Runtime type information
    TypeName("diffStencilLaplacian");

    // Constructors

        //- Construct from dictionary
        diffStencilLaplacianStab
        (
            const fvMesh& mesh,
            const dictionary& dict,
            const dimensionSet& dims
        );


    // Destructor

        virtual ~diffStencilLaplacianStab();


    // Member Functions

        // Access

            //- Return true: spectral normalisation is supported
            virtual bool supportsSpectralNormalisation() const
            {
                return true;
            }

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
#   include "diffStencilLaplacianStabTemplates.H"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
