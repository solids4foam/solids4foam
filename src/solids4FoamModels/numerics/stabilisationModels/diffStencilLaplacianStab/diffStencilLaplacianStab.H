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
    Face-based stabilisation model based on the difference between compact
    and extended stencil Laplacian/gradient approximations.

    The method constructs a stabilisation flux (surface field) proportional
    to the difference between a small-stencil face-normal gradient and a
    larger-stencil reconstructed gradient. This selectively damps grid-scale
    oscillatory modes while maintaining consistency for smooth solutions.

    The stabilisation term scales with the square of the mesh spacing (O(h^2))
    under grid refinement. Consequently, it does not alter the formal accuracy
    of a second-order discretisation. However, when applied to higher-order
    schemes, the leading truncation error introduced by this stabilisation
    is second order, and therefore the overall method may be limited to
    second-order convergence unless a higher-order variant of the
    stabilisation is employed.

    The stabilisation may be retrieved as a face flux contribution (surface
    field).

    A user-defined scale factor controls the magnitude of the stabilisation
    and can be specified at run time via dictionary entries.

    The implementation is designed for high performance: internal storage
    is preallocated and updated in-place to avoid repeated memory allocation
    during frequent residual evaluations in nonlinear or segregated solvers.

SourceFiles
    diffStencilLaplacianStab.C
    diffStencilLaplacianStabTemplates.H

Author
    Philip Cardiff, UCD. All rights reserved

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
        //  This defaults to scaleFactor
        const scalar scaleFactorJacobian_;

        //- Pointer to the scalar face stabilisation
        mutable autoPtr<surfaceScalarField> faceScalarPtr_;

        //- Pointer to the vector face stabilisation
        mutable autoPtr<surfaceVectorField> faceVectorPtr_;

        //- Pointer to the scalar approximate Jacobian
        mutable autoPtr<fvScalarMatrix> scalarJacobianPtr_;

        //- Pointer to the vector approximate Jacobian
        mutable autoPtr<fvVectorMatrix> vectorJacobianPtr_;

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

            //- Does this stabilisation model support scalar fields?
            virtual bool supportsScalar() const
            {
                return true;
            }

            //- Does this stabilisation model support vector fields?
            virtual bool supportsVector() const
            {
                return true;
            }

            //- Provide read-only access to the scalar face stabilisation
            const surfaceScalarField& faceScalar() const
            {
                if (faceScalarPtr_.empty())
                {
                    FatalErrorInFunction
                        << "Pointer not yet set: the updateScalar(...) "
                        << "function must be called first"
                        << type() << abort(FatalError);
                }

                return faceScalarPtr_();
            }

            //- Provide read-only access to the vector face stabilisation
            const surfaceVectorField& faceVector() const
            {
                if (faceVectorPtr_.empty())
                {
                    FatalErrorInFunction
                        << "Pointer not yet set: the updateVector(...) "
                        << "function must be called first"
                        << type() << abort(FatalError);
                }

                return faceVectorPtr_();
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
                const surfaceScalarField* gammaPtr = nullptr
            ) const;

            //- Provide an implicit approximation of the stabilisation
            //  This is an approximate Jacobian of the stabilisation term
            //  Provide the primary variable field and optionally a surface
            //  diffusivity field can be optionally provided
            virtual const fvVectorMatrix& vectorJacobian
            (
                const volVectorField& field,
                const surfaceScalarField* gammaPtr = nullptr
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
