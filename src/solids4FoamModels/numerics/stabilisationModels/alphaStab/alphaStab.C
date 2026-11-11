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

#include "alphaStab.H"
#include "addToRunTimeSelectionTable.H"
#include "compatibilityFunctions.H"
#include "fvmLaplacian.H"
#include "lookupSolidModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(alphaStab, 0);
    addToRunTimeSelectionTable
    (
        stabilisationModel, alphaStab, stabModel
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from dictionary
Foam::alphaStab::alphaStab
(
    const fvMesh& mesh,
    const dictionary& dict,
    const dimensionSet& dims
)
:
    stabilisationModel(mesh, dict, dims),
    scaleFactor_(readScalar(dict.lookup("scaleFactor")))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::alphaStab::~alphaStab()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::alphaStab::initFaceScalar(const volScalarField& field) const
{
    if (faceScalarPtr().empty())
    {
        faceScalarPtr().set
        (
            new surfaceScalarField
            (
                IOobject
                (
                    "faceStabilisation(" + field.name() + ")",
                    mesh().time().timeName(),
                    mesh(),
                    IOobject::READ_IF_PRESENT,
                    IOobject::NO_WRITE
                ),
                mesh(),
                dimensionedScalar("0", dims(), 0.0)
            )
        );

        // This is an oriented field
#ifdef OPENFOAM_COM
        autoPtrRef(faceScalarPtr()).setOriented(true);
#endif
    }
}


void Foam::alphaStab::initFaceVector(const volVectorField& field) const
{
    if (faceVectorPtr().empty())
    {
        faceVectorPtr().set
        (
            new surfaceVectorField
            (
                IOobject
                (
                    "faceStabilisation(" + field.name() + ")",
                    mesh().time().timeName(),
                    mesh(),
                    IOobject::READ_IF_PRESENT,
                    IOobject::NO_WRITE
                ),
                mesh(),
                dimensionedVector("0", dims(), vector::zero)
            )
        );

        // This is an oriented field
#ifdef OPENFOAM_COM
        autoPtrRef(faceVectorPtr()).setOriented(true);
#endif
    }
}


void Foam::alphaStab::updateScalar
(
    const volScalarField& p,
    const volVectorField* gradPtr
) const
{
    if (gradPtr == nullptr)
    {
        FatalErrorInFunction
            << "grad(" << p.name() << ") must be provided with this "
            << "stabilisation method" << exit(FatalError);
    }

    // A high-order solid model on this mesh supplies its pressure scheme
    if (foundSolidModel(p.mesh()))
    {
        const solidModel& solMod = lookupSolidModel(p.mesh());

        if (solMod.highOrderResidual())
        {
            updateScalarHighOrder(p, solMod.pressureLeastSquares());
            return;
        }
    }

    clearCellScalarCache();

    // If required, initialise the face stabilisation field
    initFaceScalar(p);

    // Update the second-order stabilisation
    computeDiffStencil(p, *gradPtr, autoPtrRef(faceScalarPtr()), scaleFactor_);
}


void Foam::alphaStab::updateVector
(
    const volVectorField& p,
    const volTensorField* gradPtr
) const
{
    if (gradPtr == nullptr)
    {
        FatalErrorInFunction
            << "grad(" << p.name() << ") must be provided with this "
            << "stabilisation method" << exit(FatalError);
    }

    // A high-order solid model on this mesh supplies its displacement scheme
    if (foundSolidModel(p.mesh()))
    {
        const solidModel& solMod = lookupSolidModel(p.mesh());

        if (solMod.highOrderResidual())
        {
            updateVectorHighOrder(p, solMod.displacementLeastSquares());
            return;
        }
    }

    clearCellVectorCache();

    // If required, initialise the face stabilisation field
    initFaceVector(p);

    // Update the second-order stabilisation
    computeDiffStencil(p, *gradPtr, autoPtrRef(faceVectorPtr()), scaleFactor_);
}


void Foam::alphaStab::updateScalarHighOrder
(
    const volScalarField& field,
    const leastSquaresScheme& reconstruction
) const
{
    clearCellScalarCache();

    // If required, initialise the face stabilisation field
    initFaceScalar(field);

    // Update the high-order stabilisation
    computeDiffStencilHighOrder
    (
        field,
        autoPtrRef(faceScalarPtr()),
        scaleFactor_,
        reconstruction
    );
}


void Foam::alphaStab::updateVectorHighOrder
(
    const volVectorField& field,
    const leastSquaresScheme& reconstruction
) const
{
    clearCellVectorCache();

    // If required, initialise the face stabilisation field
    initFaceVector(field);

    // Update the high-order stabilisation
    computeDiffStencilHighOrder
    (
        field,
        autoPtrRef(faceVectorPtr()),
        scaleFactor_,
        reconstruction
    );
}


const Foam::fvScalarMatrix& Foam::alphaStab::scalarJacobian
(
    const volScalarField& field,
    const surfaceScalarField* gammaPtr,
    const bool rebuild
) const
{
    checkGamma(gammaPtr);

    // Scalar Jacobian taken from diffStencilLaplacianStab, as for the vector
    // Jacobian below
    if (scalarJacobianPtr().empty() || rebuild)
    {
        if (gammaPtr == nullptr)
        {
            scalarJacobianPtr().reset
            (
                new fvScalarMatrix
                (
                    scaleFactorJacobian()
                   *fvm::laplacian(field, "laplacian(" + field.name() + ")")
                )
            );
        }
        else
        {
            const word schemeName
            (
                "laplacian(" + gammaPtr->name() + "," + field.name() + ")"
            );

            scalarJacobianPtr().reset
            (
                new fvScalarMatrix
                (
                    scaleFactorJacobian()
                   *fvm::laplacian(*gammaPtr, field, schemeName)
                )
            );
        }
    }

    return scalarJacobianPtr()();
}


const Foam::fvVectorMatrix& Foam::alphaStab::vectorJacobian
(
    const volVectorField& field,
    const surfaceScalarField* gammaPtr,
    const bool rebuild
) const
{
    checkGamma(gammaPtr);

    // Vector Jacobian taken from diffStencilLaplacianStab
    // Residual can use alpha stabilisation but Jacobian is build using
    // implicit Laplacian. Alpha implicit Laplacian will be similar so this
    // is ok workaround
    if (vectorJacobianPtr().empty() || rebuild)
    {
        if (gammaPtr == nullptr)
        {
            vectorJacobianPtr().reset
            (
                new fvVectorMatrix
                (
                    scaleFactorJacobian()
                   *fvm::laplacian(field, "laplacian(" + field.name() + ")")
                )
            );
        }
        else
        {
            const word schemeName
            (
                "laplacian(" + gammaPtr->name() + "," + field.name() + ")"
            );

            vectorJacobianPtr().reset
            (
                new fvVectorMatrix
                (
                    scaleFactorJacobian()
                   *fvm::laplacian(*gammaPtr, field, schemeName)
                )
            );
        }
    }

    return vectorJacobianPtr()();
}

// ************************************************************************* //
