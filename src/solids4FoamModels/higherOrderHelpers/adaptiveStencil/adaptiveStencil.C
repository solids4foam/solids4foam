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

#include "adaptiveStencil.H"
#include "fvMesh.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "compatibilityFunctions.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::adaptiveStencil::adaptiveStencil
(
    const dictionary& schemeDict,
    const label nTerms
)
:
    active_(schemeDict.found("adaptiveStencilCoeffs")),
    maxConditionNumber_(GREAT),
    extraCellsPerPass_(max(label(1), nTerms/2)),
    maxExtraCells_(3*nTerms),
    normalise_(true)
{
    if (active_)
    {
        const dictionary& dict = schemeDict.subDict("adaptiveStencilCoeffs");

        maxConditionNumber_ = readScalar(dict.lookup("maxConditionNumber"));
        extraCellsPerPass_ =
            dict.lookupOrDefault<label>
            (
                "extraCellsPerPass",
                extraCellsPerPass_
            );
        maxExtraCells_ =
            dict.lookupOrDefault<label>("maxExtraCells", maxExtraCells_);
        normalise_ = dict.lookupOrDefault<Switch>("normalise", normalise_);

        if (maxConditionNumber_ <= 1.0)
        {
            FatalErrorInFunction
                << "maxConditionNumber should be larger than 1, but it is "
                << "set to " << maxConditionNumber_
                << abort(FatalError);
        }

        if (extraCellsPerPass_ < 1 || maxExtraCells_ < 1)
        {
            FatalErrorInFunction
                << "extraCellsPerPass and maxExtraCells should be positive, "
                << "but extraCellsPerPass is set to " << extraCellsPerPass_
                << " and maxExtraCells to " << maxExtraCells_
                << abort(FatalError);
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::adaptiveStencil::conditionEstimate
(
    const Eigen::MatrixXd& Q,
    const Eigen::DiagonalMatrix<double, Eigen::Dynamic>& W
) const
{
    // The least-squares system is underdetermined
    if (Q.cols() < Q.rows())
    {
        return GREAT;
    }

    // Weighted least-squares matrix
    Eigen::MatrixXd Qhat = Q*W.diagonal().cwiseSqrt().asDiagonal();

    if (normalise_)
    {
        // Scale each polynomial term to unit norm
        for (Eigen::Index termI = 0; termI < Qhat.rows(); ++termI)
        {
            const scalar termNorm = Qhat.row(termI).norm();

            // A polynomial term is not seen by the stencil
            if (termNorm < VSMALL)
            {
                return GREAT;
            }

            Qhat.row(termI) /= termNorm;
        }
    }

    // Diagonal of R from the column-pivoted QR decomposition has
    // non-increasing magnitudes, and the ratio of the largest and the smallest
    // is a lower bound of the condition number
    const Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr(Qhat.transpose());

    const Eigen::VectorXd diagR = qr.matrixQR().diagonal().cwiseAbs();

    const scalar estimate = diagR.maxCoeff()/(diagR.minCoeff() + VSMALL);

    // Guard against not-a-number and overflow
    if (!(estimate < GREAT))
    {
        return GREAT;
    }

    return estimate;
}


void Foam::adaptiveStencil::writeCellField
(
    const fvMesh& mesh,
    const word& name,
    const scalarField& values
)
{
    volScalarField field
    (
        IOobject
        (
            name,
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("0", dimless, 0.0),
        "zeroGradient"
    );

    primitiveFieldRef(field) = values;
    field.correctBoundaryConditions();

    Info<< "Writing " << field.name() << " to time = "
        << mesh.time().value() << endl;
    field.write();
}


void Foam::adaptiveStencil::writeFaceField
(
    const fvMesh& mesh,
    const word& name,
    const scalarField& values
)
{
    surfaceScalarField field
    (
        IOobject
        (
            name,
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("0", dimless, 0.0)
    );

    scalarField& fieldI = primitiveFieldRef(field);

    forAll(fieldI, faceI)
    {
        fieldI[faceI] = values[faceI];
    }

    forAll(mesh.boundaryMesh(), patchI)
    {
        const polyPatch& pp = mesh.boundaryMesh()[patchI];
        scalarField& patchField = boundaryFieldRef(field)[patchI];

        forAll(patchField, i)
        {
            patchField[i] = values[pp.start() + i];
        }
    }

    Info<< "Writing " << field.name() << " to time = "
        << mesh.time().value() << endl;
    field.write();
}


// ************************************************************************* //
