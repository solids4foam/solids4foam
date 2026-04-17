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
    Foam::stabilisationModel

Description
    Base class for runtime-selectable scalar and vector stabilisation models.

    Derived classes cache face-based stabilisation fields and may also provide
    cached cell-centred fields and approximate Jacobians.

SourceFiles
    stabilisationModel.C

Author
    Philip Cardiff, UCD.

\*---------------------------------------------------------------------------*/

#ifndef stabilisationModel_H
#define stabilisationModel_H

#include "typeInfo.H"
#include "runTimeSelectionTables.H"
#include "dictionary.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "fvMatrices.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class stabilisationModel Declaration
\*---------------------------------------------------------------------------*/

class stabilisationModel
{
    // Private data

        //- Reference to the mesh
        const fvMesh& mesh_;

        //- Stabilisation dictionary
        const dictionary dict_;

        //- Dimensioned of the stabilisation
        const dimensionSet dims_;

        //- Scale factor used for globally scaling the amount of stabilisation
        const scalar scaleFactor_;

        //- Scale factor used for globally scaling the approximate Jacobian
        //  This defaults to 1.0
        const scalar scaleFactorJacobian_;

        //- Pointer to the scalar face stabilisation
        mutable autoPtr<surfaceScalarField> faceScalarPtr_;

        //- Pointer to the vector face stabilisation
        mutable autoPtr<surfaceVectorField> faceVectorPtr_;

        //- Pointer to the scalar approximate Jacobian
        mutable autoPtr<fvScalarMatrix> scalarJacobianPtr_;

        //- Pointer to the vector approximate Jacobian
        mutable autoPtr<fvVectorMatrix> vectorJacobianPtr_;

        //- Pointer to the cell scalar stabilisation
        mutable autoPtr<volScalarField> cellScalarPtr_;

        //- Pointer to the cell vector stabilisation
        mutable autoPtr<volVectorField> cellVectorPtr_;

        //- Pointer to the reciprocal of the deltaCoeffs squared
        mutable autoPtr<surfaceScalarField> h2Ptr_;

    protected:

        //- Invalidate the cached cell scalar stabilisation
        void clearCellScalarCache() const
        {
            cellScalarPtr_.clear();
        }

        //- Invalidate the cached cell vector stabilisation
        void clearCellVectorCache() const
        {
            cellVectorPtr_.clear();
        }

        //- Check that an optional face diffusivity field is non-negative
        void checkGamma(const surfaceScalarField* gammaPtr) const;


    // Private Member Functions

        //- Disallow copy construct
        stabilisationModel(const stabilisationModel&);

        //- Disallow default bitwise assignment
        void operator=(const stabilisationModel&);


public:

    //- Runtime type information
    TypeName("stabilisationModel");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            stabilisationModel,
            stabModel,
            (
                const fvMesh& mesh,
                const dictionary& dict,
                const dimensionSet& dims
            ),
            (mesh, dict, dims)
        );


    // Selectors

        //- Create a new stabilisation model
        static autoPtr<stabilisationModel> New
        (
            const fvMesh& mesh,
            const dictionary& dict,
            const dimensionSet& dims
        );


    // Constructors

        //- Construct from dictionary
        stabilisationModel
        (
            const fvMesh& mesh,
            const dictionary& dict,
            const dimensionSet& dims
        );


    // Destructor

        virtual ~stabilisationModel() = default;


    // Member Functions

        // Access

            //- Return const reference to the mesh
            const fvMesh& mesh() const
            {
                return mesh_;
            }

            //- Return const reference to dict
            const dictionary& dict() const
            {
                return dict_;
            }

            //- Return const reference to dimensions
            const dimensionSet& dims() const
            {
                return dims_;
            }

            //- Return the scale factor
            scalar scaleFactor() const
            {
                return scaleFactor_;
            }

            //- Return the scale factor for Jacobian
            scalar scaleFactorJacobian() const
            {
                return scaleFactorJacobian_;
            }

            //- Provide access to the scalar face stabilisation pointer
            autoPtr<surfaceScalarField>& faceScalarPtr() const
            {
                return faceScalarPtr_;
            };

            //- Provide read-only access to the scalar face stabilisation
            const surfaceScalarField& faceScalar() const
            {
                if (faceScalarPtr_.empty())
                {
                    FatalErrorInFunction
                        << "Pointer not yet set: the updateScalar(...) "
                        << "function must be called first for "
                        << type() << abort(FatalError);
                }

                return faceScalarPtr_();
            }

            //- Provide access to the vector face stabilisation pointer
            autoPtr<surfaceVectorField>& faceVectorPtr() const
            {
                return faceVectorPtr_;
            };

            //- Provide read-only access to the vector face stabilisation
            const surfaceVectorField& faceVector() const
            {
                if (faceVectorPtr_.empty())
                {
                    FatalErrorInFunction
                        << "Pointer not yet set: the updateVector(...) "
                        << "function must be called first for "
                        << type() << abort(FatalError);
                }

                return faceVectorPtr_();
            }

            //- Provide read-only access to the cell scalar stabilisation
            //  If rebuild is true, re-calculate it from the cached face field
            //  An optional face diffusivity can be included inside the div
            const volScalarField& cellScalar
            (
                const surfaceScalarField* gammaPtr = nullptr,
                const bool rebuild = false
            ) const;

            //- Provide read-only access to the cell vector stabilisation
            //  If rebuild is true, re-calculate it from the cached face field
            //  An optional face diffusivity can be included inside the div
            const volVectorField& cellVector
            (
                const surfaceScalarField* gammaPtr = nullptr,
                const bool rebuild = false
            ) const;

            //- Provide access to the scalarJacobianPtr
            autoPtr<fvScalarMatrix>& scalarJacobianPtr() const
            {
                return scalarJacobianPtr_;
            }

            //- Provide access to the vectorJacobianPtr
            autoPtr<fvVectorMatrix>& vectorJacobianPtr() const
            {
                return vectorJacobianPtr_;
            }

            //- Provide read-only access to the h2 field
            const surfaceScalarField& h2() const
            {
                if (h2Ptr_.empty())
                {
                    h2Ptr_.set
                    (
                        new surfaceScalarField
                        (
                            "h2",
                            1.0/sqr(mesh().deltaCoeffs())
                        )
                    );
                }

                return h2Ptr_();
            }

        // Edit

            //- Re-calculate the scalar stabilisation
            //  Provide the primary variable field and optionally its gradient
            virtual void updateScalar
            (
                const volScalarField& field,
                const volVectorField* gradPtr = nullptr
            ) const
            {
                FatalErrorInFunction
                    << "scalar stabilisation not supported by "
                    << type() << exit(FatalError);

                // Keep the compiler happy
                std::abort();
            }

            //- Re-calculate the vector stabilisation
            //  Provide the primary variable field and optionally its gradient
            virtual void updateVector
            (
                const volVectorField& field,
                const volTensorField* gradPtr = nullptr
            ) const
            {
                FatalErrorInFunction
                    << "vector stabilisation not supported by "
                    << type() << exit(FatalError);

                // Keep the compiler happy
                std::abort();
            }

            //- Provide an implicit approximation of the stabilisation
            //  This is an approximate Jacobian of the stabilisation term
            //  Provide the primary variable field and optionally a surface
            //  diffusivity field can be optionally provided
            virtual const fvScalarMatrix& scalarJacobian
            (
                const volScalarField& field,
                const surfaceScalarField* gammaPtr = nullptr,
                const bool rebuild = false
            ) const
            {
                FatalErrorInFunction
                    << "scalarJacobian not supported by "
                    << type() << exit(FatalError);

                // Keep the compiler happy
                std::abort();
            }

            //- Provide an implicit approximation of the stabilisation
            //  This is an approximate Jacobian of the stabilisation term
            //  Provide the primary variable field and optionally a surface
            //  diffusivity field can be optionally provided
            virtual const fvVectorMatrix& vectorJacobian
            (
                const volVectorField& field,
                const surfaceScalarField* gammaPtr = nullptr,
                const bool rebuild = false
            ) const
            {
                FatalErrorInFunction
                    << "vectorJacobian not supported by "
                    << type() << exit(FatalError);

                // Keep the compiler happy
                std::abort();
            }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
