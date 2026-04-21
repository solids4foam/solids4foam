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

\*---------------------------------------------------------------------------*/

#include "multipleStabilisationModels.H"
#include "addToRunTimeSelectionTable.H"
#include "compatibilityFunctions.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(multipleStabilisationModels, 0);
    addToRunTimeSelectionTable
    (
        stabilisationModel, multipleStabilisationModels, stabModel
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::multipleStabilisationModels::multipleStabilisationModels
(
    const fvMesh& mesh,
    const dictionary& dict,
    const dimensionSet& dims
)
:
    stabilisationModel(mesh, dict, dims),
    modelNames_(dict.lookup("models")),
    jacobianModelName_(dict.lookup("jacobianModel")),
    models_(modelNames_.size())
{
    // Validate that jacobianModelName_ is one of the listed models
    if (findIndex(modelNames_, jacobianModelName_) == -1)
    {
        FatalErrorInFunction
            << "jacobianModel '" << jacobianModelName_ << "' not found in "
            << "models list " << modelNames_
            << exit(FatalError);
    }

    // Construct each sub-model from its named sub-dictionary
    forAll(modelNames_, i)
    {
        models_.set
        (
            i,
            stabilisationModel::New
            (
                mesh,
                dict.subDict(modelNames_[i]),
                dims
            )
        );
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::multipleStabilisationModels::~multipleStabilisationModels()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::multipleStabilisationModels::updateScalar
(
    const volScalarField& field,
    const volVectorField* gradPtr
) const
{
    clearCellScalarCache();

    // Update all sub-models; a fatal error propagates if any do not support
    // scalar stabilisation
    forAll(models_, i)
    {
        models_[i].updateScalar(field, gradPtr);
    }

    // Initialise the composite face field from the first sub-model
    faceScalarPtr().reset
    (
        new surfaceScalarField
        (
            IOobject
            (
                "faceStabilisation(" + field.name() + ")",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            models_[0].faceScalar()
        )
    );

    // Accumulate contributions from remaining sub-models
    for (label i = 1; i < models_.size(); i++)
    {
        autoPtrRef(faceScalarPtr()) += models_[i].faceScalar();
    }

    // Apply the outer scale factor (defaults to 1.0 if not specified)
    autoPtrRef(faceScalarPtr()) *= scaleFactor();
}


void Foam::multipleStabilisationModels::updateVector
(
    const volVectorField& field,
    const volTensorField* gradPtr
) const
{
    clearCellVectorCache();

    // Update all sub-models; a fatal error propagates if any do not support
    // vector stabilisation
    forAll(models_, i)
    {
        models_[i].updateVector(field, gradPtr);
    }

    // Initialise the composite face field from the first sub-model
    faceVectorPtr().reset
    (
        new surfaceVectorField
        (
            IOobject
            (
                "faceStabilisation(" + field.name() + ")",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            models_[0].faceVector()
        )
    );

    // Accumulate contributions from remaining sub-models
    for (label i = 1; i < models_.size(); i++)
    {
        autoPtrRef(faceVectorPtr()) += models_[i].faceVector();
    }

    // Apply the outer scale factor (defaults to 1.0 if not specified)
    autoPtrRef(faceVectorPtr()) *= scaleFactor();
}


const Foam::fvScalarMatrix& Foam::multipleStabilisationModels::scalarJacobian
(
    const volScalarField& field,
    const surfaceScalarField* gammaPtr,
    const bool rebuild
) const
{
    const label idx = findIndex(modelNames_, jacobianModelName_);
    return models_[idx].scalarJacobian(field, gammaPtr, rebuild);
}


const Foam::fvVectorMatrix& Foam::multipleStabilisationModels::vectorJacobian
(
    const volVectorField& field,
    const surfaceScalarField* gammaPtr,
    const bool rebuild
) const
{
    const label idx = findIndex(modelNames_, jacobianModelName_);
    return models_[idx].vectorJacobian(field, gammaPtr, rebuild);
}


bool Foam::multipleStabilisationModels::supportsHighOrderResidual() const
{
    forAll(models_, i)
    {
        if (!models_[i].supportsHighOrderResidual())
        {
            return false;
        }
    }
    return true;
}


// ************************************************************************* //
