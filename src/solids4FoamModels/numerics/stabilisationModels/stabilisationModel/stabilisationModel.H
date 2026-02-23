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
    Class for calculating stabilisation terms for vector and scalar geometric
    fields.

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

            //- Does this stabilisation model support scalar fields?
            virtual bool supportsScalar() const
            {
                return false;
            }

            //- Does this stabilisation model support vector fields?
            virtual bool supportsVector() const
            {
                return false;
            }

            //- Provide read-only access to the scalar face stabilisation
            virtual const surfaceScalarField& faceScalar() const
            {
                FatalErrorInFunction
                    << "Scalar stabilisation not supported by "
                    << type() << exit(FatalError);

                // Keep the compiler happy
                std::abort();
            }

            //- Provide read-only access to the vector face stabilisation
            virtual const surfaceVectorField& faceVector() const
            {
                FatalErrorInFunction
                    << "Vector stabilisation not supported by "
                    << type() << exit(FatalError);

                // Keep the compiler happy
                std::abort();
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
            // Provide the primary variable field and optionally its gradient
            virtual void updateVector
            (
                const volVectorField& field,
                const volTensorField* gradPtr = nullptr
            ) const
            {
                FatalErrorInFunction
                    << "scalar stabilisation not supported by "
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
