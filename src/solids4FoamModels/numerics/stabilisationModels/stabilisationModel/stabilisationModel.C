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

InClass
    Foam::stabilisationModel

\*---------------------------------------------------------------------------*/

#include "stabilisationModel.H"
#include "fvc.H"
#include "compatibilityFunctions.H"

namespace Foam
{
namespace
{
label readReferenceNyquistDirections
(
    const fvMesh& mesh,
    const dictionary& dict,
    const bool normalise
)
{
    if (!normalise)
    {
        return 0;
    }

    if (!dict.found("referenceNyquistDirections"))
    {
        FatalIOErrorInFunction(dict)
            << "referenceNyquistDirections is required when normalise is true"
            << exit(FatalIOError);
    }

    const label referenceNyquistDirections = readLabel
    (
        dict.lookup("referenceNyquistDirections")
    );

    if (referenceNyquistDirections <= 0)
    {
        FatalIOErrorInFunction(dict)
            << "referenceNyquistDirections must be greater than zero, found "
            << referenceNyquistDirections << exit(FatalIOError);
    }

    const label maxReferenceNyquistDirections =
        min(label(3), mesh.nSolutionD());

    if (referenceNyquistDirections > maxReferenceNyquistDirections)
    {
        FatalIOErrorInFunction(dict)
            << "referenceNyquistDirections must be in the valid range 1 .. "
            << maxReferenceNyquistDirections
            << " = min(3, mesh.nSolutionD()), found "
            << referenceNyquistDirections << exit(FatalIOError);
    }

    return referenceNyquistDirections;
}
}
}

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(stabilisationModel, 0);
    defineRunTimeSelectionTable(stabilisationModel, stabModel);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::stabilisationModel::stabilisationModel
(
    const fvMesh& mesh,
    const dictionary& dict,
    const dimensionSet& dims
)
:
    mesh_(mesh),
    dict_(dict),
    dims_(dims),
    scaleFactor_(dict.lookupOrDefault<scalar>("scaleFactor", 1.0)),
    scaleFactorJacobian_
    (
        dict.lookupOrDefault<scalar>("scaleFactorJacobian", 1.0)
    ),
    normalise_(dict.lookupOrDefault<Switch>("normalise", false)),
    referenceNyquistDirections_
    (
        readReferenceNyquistDirections(mesh, dict, normalise_)
    ),
    faceScalarPtr_(),
    faceVectorPtr_(),
    scalarJacobianPtr_(),
    vectorJacobianPtr_(),
    cellScalarPtr_(),
    cellVectorPtr_(),
    h2Ptr_()
{
    if (!normalise_ && dict.found("referenceNyquistDirections"))
    {
        WarningInFunction
            << "Ignoring referenceNyquistDirections because normalise is false"
            << endl;
    }
}


// * * * * * * * * * * * * * * * Member functions * * * * * * * * * * * * * * //


Foam::autoPtr<Foam::stabilisationModel> Foam::stabilisationModel::New
(
    const fvMesh& mesh,
    const dictionary& dict,
    const dimensionSet& dims
)
{
    const word modelType(dict.lookup("type"));

    Info<< "Selecting stabilisation model " << modelType << endl;

#if (OPENFOAM >= 2112)
    auto* ctorPtr = stabModelConstructorTable(modelType);

    if (!ctorPtr)
    {
        FatalIOErrorInLookup
        (
            dict,
            "stabilisationModel",
            modelType,
            *stabModelConstructorTablePtr_
        )   << exit(FatalIOError);
    }
#else
    stabModelConstructorTable::iterator cstrIter =
        stabModelConstructorTablePtr_->find(modelType);

    if (cstrIter == stabModelConstructorTablePtr_->end())
    {
        FatalErrorIn("stabilisationModel::New(...)")
            << "Unknown stabilisation model type " << modelType
            << nl << nl
            << "Valid stabilisation model types are:" << nl
            << stabModelConstructorTablePtr_->toc()
            << exit(FatalError);
    }

    auto* ctorPtr = cstrIter();
#endif

    autoPtr<stabilisationModel> modelPtr(ctorPtr(mesh, dict, dims));

    if (modelPtr->normalise() && !modelPtr->supportsSpectralNormalisation())
    {
        FatalIOErrorInFunction(dict)
            << "Spectral normalisation is not supported by stabilisation "
            << "model " << modelType << exit(FatalIOError);
    }

    return modelPtr;
}


Foam::scalar Foam::stabilisationModel::nyquistNormalisation
(
    const label paperOrder
) const
{
    if (!normalise_)
    {
        return 1.0;
    }

    scalar referenceResponse = 1.0;
    const scalar directionalResponse = 4.0*referenceNyquistDirections_;

    for (label i = 0; i < paperOrder; ++i)
    {
        referenceResponse *= directionalResponse;
    }

    return 1.0/referenceResponse;
}


void Foam::stabilisationModel::writeSpectralNormalisationInfo
(
    const string& responseDescription,
    const scalar factor
) const
{
    if (normalise_)
    {
        Info<< "Spectral normalisation: model " << word(dict_.lookup("type"))
            << ", referenceNyquistDirections "
            << referenceNyquistDirections_
            << ", " << responseDescription.c_str()
            << ", factor " << factor << endl;
    }
}


void Foam::stabilisationModel::writeSpectralNormalisationInfo
(
    const label paperOrder,
    const bool usesIdealLaplacianStencil
) const
{
    string responseDescription("paper order ");
    responseDescription += Foam::name(paperOrder);

    if (usesIdealLaplacianStencil)
    {
        responseDescription += ", using ideal m=1 snGrad stencil";
    }

    writeSpectralNormalisationInfo
    (
        responseDescription,
        nyquistNormalisation(paperOrder)
    );
}


void Foam::stabilisationModel::checkGamma
(
    const surfaceScalarField* gammaPtr
) const
{
    if (gammaPtr != nullptr)
    {
        const scalar minGamma = gMin(*gammaPtr);

        if (minGamma < 0.0)
        {
            FatalErrorInFunction
                << "Negative gamma field " << gammaPtr->name()
                << " supplied to stabilisation model " << type()
                << ": min(gamma) = " << minGamma
                << abort(FatalError);
        }
    }
}


const Foam::volScalarField& Foam::stabilisationModel::cellScalar
(
    const surfaceScalarField* gammaPtr,
    const bool rebuild
) const
{
    checkGamma(gammaPtr);

    if (cellScalarPtr_.empty() || rebuild)
    {
        if (faceScalarPtr_.empty())
        {
            FatalErrorInFunction
                << "Pointer not yet set: the updateScalar(...) "
                << "function must be called first for "
                << type() << abort(FatalError);
        }

        tmp<volScalarField> tCellScalar
        (
            gammaPtr == nullptr
          ? fvc::div(mesh().magSf()*faceScalar())
          : fvc::div((*gammaPtr)*mesh().magSf()*faceScalar())
        );
        tmpRef(tCellScalar).rename
        (
            "cellStabilisation(" + faceScalar().name() + ")"
        );
        cellScalarPtr_.reset(tCellScalar.ptr());
    }

    return cellScalarPtr_();
}


const Foam::volVectorField& Foam::stabilisationModel::cellVector
(
    const surfaceScalarField* gammaPtr,
    const bool rebuild
) const
{
    checkGamma(gammaPtr);

    if (cellVectorPtr_.empty() || rebuild)
    {
        if (faceVectorPtr_.empty())
        {
            FatalErrorInFunction
                << "Pointer not yet set: the updateVector(...) "
                << "function must be called first for "
                << type() << abort(FatalError);
        }

        tmp<volVectorField> tCellVector
        (
            gammaPtr == nullptr
          ? fvc::div(mesh().magSf()*faceVector())
          : fvc::div((*gammaPtr)*mesh().magSf()*faceVector())
        );
        tmpRef(tCellVector).rename
        (
            "cellStabilisation(" + faceVector().name() + ")"
        );
        cellVectorPtr_.reset(tCellVector.ptr());
    }

    return cellVectorPtr_();
}


// ************************************************************************* //
