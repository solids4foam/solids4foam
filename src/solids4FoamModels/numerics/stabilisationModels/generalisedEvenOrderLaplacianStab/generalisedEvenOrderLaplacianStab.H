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
    Face-based stabilisation model based on the Laplacian^laplacianPower. This
    is a generalisation of the 4th order term of Jameson, Schmidt and Turkel.

    This selectively damps grid-scale oscillatory modes while maintaining
    consistency for smooth solutions.

    For similar damping of the checkboard frequencyies, set the scaleFactor as:

        scaleFactor = (4d)^{1 - m}

        where d is the dimension (1-D, 2-D, 3-D), and m == laplacianPower

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

        //- Pointer to the scalar face stabilisation
        mutable autoPtr<surfaceScalarField> faceScalarPtr_;

        //- Pointer to the vector face stabilisation
        mutable autoPtr<surfaceVectorField> faceVectorPtr_;

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
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#ifdef NoRepository
#   include "generalisedEvenOrderLaplacianStabTemplates.H"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
