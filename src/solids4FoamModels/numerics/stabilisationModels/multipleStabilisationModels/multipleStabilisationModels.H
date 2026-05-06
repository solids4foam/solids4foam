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
    multipleStabilisationModels

Description
    Composite stabilisation model that combines two or more stabilisation
    models at run time.

    Each sub-model is constructed from a named sub-dictionary. Their face
    contributions are summed and an optional outer scaleFactor (default 1.0)
    is applied to the accumulated result.

    The Jacobian is delegated to a single user-nominated sub-model specified
    via the jacobianModel entry.

    Example dictionary:

    \verbatim
    momentum
    {
        type            multiple;
        // scaleFactor 1.0;  // optional outer scale; defaults to 1.0
        models          ( stab1  stab2 );
        jacobianModel   stab1;

        stab1
        {
            type        laplacian;
            scaleFactor 0.1;
        }

        stab2
        {
            type        gradDiv;
            scaleFactor 0.05;
        }
    }
    \endverbatim

SourceFiles
    multipleStabilisationModels.C

Author
    Philip Cardiff, UCD.

\*---------------------------------------------------------------------------*/

#ifndef multipleStabilisationModels_H
#define multipleStabilisationModels_H

#include "stabilisationModel.H"
#include "PtrList.H"
#include "wordList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class multipleStabilisationModels Declaration
\*---------------------------------------------------------------------------*/

class multipleStabilisationModels
:
    public stabilisationModel
{
    // Private data

        //- Ordered list of sub-model names matching models_
        wordList modelNames_;

        //- Name of the sub-model that provides the Jacobian
        word jacobianModelName_;

        //- List of sub-models
        PtrList<stabilisationModel> models_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        multipleStabilisationModels(const multipleStabilisationModels&);

        //- Disallow default bitwise assignment
        void operator=(const multipleStabilisationModels&);


public:

    //- Runtime type information
    TypeName("multiple");


    // Constructors

        //- Construct from dictionary
        multipleStabilisationModels
        (
            const fvMesh& mesh,
            const dictionary& dict,
            const dimensionSet& dims
        );


    // Destructor

        virtual ~multipleStabilisationModels();


    // Member Functions

        //- Re-calculate the scalar stabilisation by accumulating all
        //  sub-model contributions
        virtual void updateScalar
        (
            const volScalarField& field,
            const volVectorField* gradPtr = nullptr
        ) const;

        //- Re-calculate the vector stabilisation by accumulating all
        //  sub-model contributions
        virtual void updateVector
        (
            const volVectorField& field,
            const volTensorField* gradPtr = nullptr
        ) const;

        //- Delegate the scalar Jacobian to the nominated sub-model
        virtual const fvScalarMatrix& scalarJacobian
        (
            const volScalarField& field,
            const surfaceScalarField* gammaPtr = nullptr,
            const bool rebuild = false
        ) const;

        //- Delegate the vector Jacobian to the nominated sub-model
        virtual const fvVectorMatrix& vectorJacobian
        (
            const volVectorField& field,
            const surfaceScalarField* gammaPtr = nullptr,
            const bool rebuild = false
        ) const;

        //- Return true only if all sub-models support high-order residual
        //  calculation
        virtual bool supportsHighOrderResidual() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
