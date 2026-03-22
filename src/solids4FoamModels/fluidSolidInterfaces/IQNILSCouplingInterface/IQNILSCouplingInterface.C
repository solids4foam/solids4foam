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

#include "IQNILSCouplingInterface.H"
#include "addToRunTimeSelectionTable.H"
#include "RectangularMatrix.H"
#include <limits>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace fluidSolidInterfaces
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(IQNILSCouplingInterface, 0);
addToRunTimeSelectionTable
(
    fluidSolidInterface, IQNILSCouplingInterface, dictionary
);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

label IQNILSCouplingInterface::couplingReuse() const
{
    return couplingReuse_;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

IQNILSCouplingInterface::IQNILSCouplingInterface
(
    Time& runTime,
    const word& region
)
:
    fluidSolidInterface(typeName, runTime, region),
    relaxationFactor_
    (
        fsiProperties().lookupOrAddDefault<scalar>("relaxationFactor", 0.01)
    ),
    couplingReuse_(fsiProperties().lookupOrAddDefault<int>("couplingReuse", 0)),
    minSignificant_
    (
        fsiProperties().lookupOrAddDefault<scalar>("minSignificant", 1.0e-12)
    ),
    relMinSignificant_
    (
        fsiProperties().lookupOrAddDefault<scalar>
        (
            "relMinSignificant", -1.0
        )
    ),
    normalizeCouplingColumns_
    (
        fsiProperties().lookupOrAddDefault<bool>
        (
            "normalizeCouplingColumns", false
        )
    ),
    qrSolveTolerance_
    (
        fsiProperties().lookupOrAddDefault<scalar>
        (
            "qrSolveTolerance", 1.0e-10
        )
    ),
    qrFilterTolerance_
    (
        fsiProperties().lookupOrAddDefault<scalar>
        (
            "qrFilterTolerance", -1.0
        )
    ),
    combinedCouplingSystem_
    (
        fsiProperties().lookupOrAddDefault<bool>
        (
            "combinedCouplingSystem", false
        )
    ),
    preciceStyleCouplingQR_
    (
        fsiProperties().lookupOrAddDefault<bool>
        (
            "preciceStyleCouplingQR", false
        )
    ),
    reusePreviousStepModes_
    (
        fsiProperties().lookupOrAddDefault<bool>
        (
            "reusePreviousStepModes", false
        )
    ),
    minPreviousStepModesForReuse_
    (
        fsiProperties().lookupOrAddDefault<label>
        (
            "minPreviousStepModesForReuse", 1
        )
    ),
    maxReuseUpdateNormRatio_
    (
        fsiProperties().lookupOrAddDefault<scalar>
        (
            "maxReuseUpdateNormRatio", -1.0
        )
    ),
    reorthogonalizeCouplingColumns_
    (
        fsiProperties().lookupOrAddDefault<bool>
        (
            "reorthogonalizeCouplingColumns", false
        )
    ),
    residualSumPreconditioning_
    (
        fsiProperties().lookupOrAddDefault<bool>
        (
            "residualSumPreconditioning", false
        )
    ),
    residualPreconditioningEpsilon_
    (
        fsiProperties().lookupOrAddDefault<scalar>
        (
            "residualPreconditioningEpsilon", SMALL
        )
    ),
    predictSolid_(fsiProperties().lookupOrAddDefault<bool>("predictSolid", true)),
    fluidPatchesPointsV_(nGlobalPatches()),
    fluidPatchesPointsW_(nGlobalPatches()),
    solidPatchesFacesTractionW_(nGlobalPatches()),
    fluidPatchesPointsT_(nGlobalPatches()),
    previousStepPatchesPointsV_(nGlobalPatches()),
    previousStepPatchesPointsW_(nGlobalPatches()),
    previousStepPatchesFacesTractionW_(nGlobalPatches()),
    previousNonEmptyStepPatchesPointsV_(nGlobalPatches()),
    previousNonEmptyStepPatchesPointsW_(nGlobalPatches()),
    previousNonEmptyStepPatchesFacesTractionW_(nGlobalPatches()),
    residualPointScale_(nGlobalPatches())
{}


void IQNILSCouplingInterface::removeCouplingMode
(
    const label interfaceI,
    const label modeI
)
{
    for
    (
        label i = modeI;
        i < fluidPatchesPointsT_[interfaceI].size() - 1;
        i++
    )
    {
        fluidPatchesPointsT_[interfaceI][i] =
            fluidPatchesPointsT_[interfaceI][i + 1];

        fluidPatchesPointsV_[interfaceI][i] =
            fluidPatchesPointsV_[interfaceI][i + 1];

        fluidPatchesPointsW_[interfaceI][i] =
            fluidPatchesPointsW_[interfaceI][i + 1];

        solidPatchesFacesTractionW_[interfaceI][i] =
            solidPatchesFacesTractionW_[interfaceI][i + 1];
    }

    fluidPatchesPointsT_[interfaceI].remove();
    fluidPatchesPointsV_[interfaceI].remove();
    fluidPatchesPointsW_[interfaceI].remove();
    solidPatchesFacesTractionW_[interfaceI].remove();
}


void IQNILSCouplingInterface::removeOldCouplingModes
(
    const label interfaceI,
    const label minTimeIndexToKeep
)
{
    while
    (
        fluidPatchesPointsT_[interfaceI].size()
     && fluidPatchesPointsT_[interfaceI][0] < minTimeIndexToKeep
    )
    {
        removeCouplingMode(interfaceI, 0);
    }
}


void IQNILSCouplingInterface::clearCouplingModes(const label interfaceI)
{
    fluidPatchesPointsT_[interfaceI].clear();
    fluidPatchesPointsV_[interfaceI].clear();
    fluidPatchesPointsW_[interfaceI].clear();
    solidPatchesFacesTractionW_[interfaceI].clear();
}


void IQNILSCouplingInterface::cacheCurrentStepModes()
{
    if (!(reusePreviousStepModes_ && couplingReuse() == 0))
    {
        return;
    }

    forAll(fluid().globalPatches(), interfaceI)
    {
        if (fluidPatchesPointsV_[interfaceI].size() > 0)
        {
            previousNonEmptyStepPatchesPointsV_[interfaceI].clear();
            previousNonEmptyStepPatchesPointsW_[interfaceI].clear();
            previousNonEmptyStepPatchesFacesTractionW_[interfaceI].clear();

            forAll(fluidPatchesPointsV_[interfaceI], modeI)
            {
                previousNonEmptyStepPatchesPointsV_[interfaceI].append
                (
                    fluidPatchesPointsV_[interfaceI][modeI]
                );
                previousNonEmptyStepPatchesPointsW_[interfaceI].append
                (
                    fluidPatchesPointsW_[interfaceI][modeI]
                );
                previousNonEmptyStepPatchesFacesTractionW_[interfaceI].append
                (
                    solidPatchesFacesTractionW_[interfaceI][modeI]
                );
            }
        }

        previousStepPatchesPointsV_[interfaceI].clear();
        previousStepPatchesPointsW_[interfaceI].clear();
        previousStepPatchesFacesTractionW_[interfaceI].clear();

        forAll(fluidPatchesPointsV_[interfaceI], modeI)
        {
            previousStepPatchesPointsV_[interfaceI].append
            (
                fluidPatchesPointsV_[interfaceI][modeI]
            );
            previousStepPatchesPointsW_[interfaceI].append
            (
                fluidPatchesPointsW_[interfaceI][modeI]
            );
            previousStepPatchesFacesTractionW_[interfaceI].append
            (
                solidPatchesFacesTractionW_[interfaceI][modeI]
            );
        }

        clearCouplingModes(interfaceI);
    }
}


bool IQNILSCouplingInterface::usePreviousStepModesInSolve
(
    const label interfaceI
) const
{
    return
        reusePreviousStepModes_
     && couplingReuse() == 0
     && outerCorr() == 1
     && previousStepPatchesPointsV_[interfaceI].size()
        >= minPreviousStepModesForReuse_;
}


bool IQNILSCouplingInterface::usePreviousNonEmptyStepModesInSolve
(
    const label interfaceI
) const
{
    return
        reusePreviousStepModes_
     && couplingReuse() == 0
     && outerCorr() == 1
     && previousStepPatchesPointsV_[interfaceI].size() == 0
     && previousNonEmptyStepPatchesPointsV_[interfaceI].size()
        >= minPreviousStepModesForReuse_;
}


label IQNILSCouplingInterface::nSolveCouplingModes(const label interfaceI) const
{
    label nModes = fluidPatchesPointsV_[interfaceI].size();

    if (usePreviousStepModesInSolve(interfaceI))
    {
        nModes += previousStepPatchesPointsV_[interfaceI].size();
    }
    else if (usePreviousNonEmptyStepModesInSolve(interfaceI))
    {
        nModes += previousNonEmptyStepPatchesPointsV_[interfaceI].size();
    }

    return nModes;
}


void IQNILSCouplingInterface::solveCouplingMode
(
    const label interfaceI,
    const label newestModeI,
    vectorField& v,
    vectorField& w
) const
{
    const label nActiveModes = fluidPatchesPointsV_[interfaceI].size();

    if (newestModeI < nActiveModes)
    {
        const label modeI = nActiveModes - 1 - newestModeI;
        v = fluidPatchesPointsV_[interfaceI][modeI];
        w = fluidPatchesPointsW_[interfaceI][modeI];
    }
    else
    {
        const label previousNewestModeI = newestModeI - nActiveModes;
        const bool usePreviousStepModes =
            usePreviousStepModesInSolve(interfaceI);
        const label nPreviousModes =
            usePreviousStepModes
          ? previousStepPatchesPointsV_[interfaceI].size()
          : previousNonEmptyStepPatchesPointsV_[interfaceI].size();
        const label modeI = nPreviousModes - 1 - previousNewestModeI;

        if (usePreviousStepModes)
        {
            v = previousStepPatchesPointsV_[interfaceI][modeI];
            w = previousStepPatchesPointsW_[interfaceI][modeI];
        }
        else
        {
            v = previousNonEmptyStepPatchesPointsV_[interfaceI][modeI];
            w = previousNonEmptyStepPatchesPointsW_[interfaceI][modeI];
        }
    }
}


void IQNILSCouplingInterface::solveSecondaryTractionMode
(
    const label interfaceI,
    const label newestModeI,
    vectorField& wTraction
) const
{
    const label nActiveModes = solidPatchesFacesTractionW_[interfaceI].size();

    if (newestModeI < nActiveModes)
    {
        const label modeI = nActiveModes - 1 - newestModeI;
        wTraction = solidPatchesFacesTractionW_[interfaceI][modeI];
    }
    else
    {
        const label previousNewestModeI = newestModeI - nActiveModes;
        const bool usePreviousStepModes =
            usePreviousStepModesInSolve(interfaceI);
        const label nPreviousModes =
            usePreviousStepModes
          ? previousStepPatchesFacesTractionW_[interfaceI].size()
          : previousNonEmptyStepPatchesFacesTractionW_[interfaceI].size();
        const label modeI = nPreviousModes - 1 - previousNewestModeI;

        if (usePreviousStepModes)
        {
            wTraction = previousStepPatchesFacesTractionW_[interfaceI][modeI];
        }
        else
        {
            wTraction =
                previousNonEmptyStepPatchesFacesTractionW_[interfaceI][modeI];
        }
    }
}


label IQNILSCouplingInterface::combinedCouplingSystemSize()
{
    label size = 0;

    forAll(fluid().globalPatches(), interfaceI)
    {
        size += fluidZonesPointsDispls()[interfaceI].size();
    }

    return size;
}


label IQNILSCouplingInterface::nCombinedCouplingModes()
{
    label cols = 0;

    forAll(fluid().globalPatches(), interfaceI)
    {
        cols = max(cols, nSolveCouplingModes(interfaceI));
    }

    return cols;
}


void IQNILSCouplingInterface::combinedCouplingMode
(
    const label newestModeI,
    vectorField& v,
    vectorField& w
)
{
    const label totalSize = combinedCouplingSystemSize();

    v.setSize(totalSize);
    v = Zero;
    w.setSize(totalSize);
    w = Zero;

    label offset = 0;

    forAll(fluid().globalPatches(), interfaceI)
    {
        const label localSize = fluidZonesPointsDispls()[interfaceI].size();

        if (newestModeI < nSolveCouplingModes(interfaceI))
        {
            vectorField localV;
            vectorField localW;
            solveCouplingMode(interfaceI, newestModeI, localV, localW);

            forAll(localV, pointI)
            {
                v[offset + pointI] = localV[pointI];
                w[offset + pointI] = localW[pointI];
            }
        }

        offset += localSize;
    }
}


vectorField IQNILSCouplingInterface::combinedResidual()
{
    vectorField residual(combinedCouplingSystemSize(), Zero);

    label offset = 0;

    forAll(fluid().globalPatches(), interfaceI)
    {
        const vectorField localResidual
        (
            fluidZonesPointsDispls()[interfaceI]
          - solidZonesPointsDispls()[interfaceI]
        );

        forAll(localResidual, pointI)
        {
            residual[offset + pointI] = localResidual[pointI];
        }

        offset += localResidual.size();
    }

    return residual;
}


void IQNILSCouplingInterface::applyCombinedCorrection
(
    const vectorField& correction
)
{
    label offset = 0;

    forAll(fluid().globalPatches(), interfaceI)
    {
        fluidZonesPointsDisplsPrev()[interfaceI] =
            fluidZonesPointsDispls()[interfaceI];

        fluidZonesPointsDispls()[interfaceI] =
            solidZonesPointsDispls()[interfaceI];

        forAll(fluidZonesPointsDispls()[interfaceI], pointI)
        {
            fluidZonesPointsDispls()[interfaceI][pointI] +=
                correction[offset + pointI];
        }

        offset += fluidZonesPointsDispls()[interfaceI].size();
    }
}


void IQNILSCouplingInterface::filterRelativeCouplingModes(const label interfaceI)
{
    if
    (
        relMinSignificant_ <= 0.0
     || fluidPatchesPointsV_[interfaceI].size() < 2
    )
    {
        return;
    }

    const word patchName =
        fluidMesh().boundary()
        [
            fluid().globalPatches()[interfaceI].patch().index()
        ].name();

    bool removedMode = true;

    while (removedMode && fluidPatchesPointsV_[interfaceI].size() > 1)
    {
        removedMode = false;

        const label cols = fluidPatchesPointsV_[interfaceI].size();
        DynamicList<vectorField> Q(cols);

        for (label i = 0; i < cols; i++)
        {
            vectorField v(fluidPatchesPointsV_[interfaceI][cols - 1 - i]);
            const scalar origNorm = Foam::sqrt(sum(v & v));
            scalarField r;
            scalar orthNorm = 0.0;
            const label status =
                orthogonalizeCouplingColumn(v, r, orthNorm, Q);

            // Keep the newest column and only filter older information.
            if
            (
                i > 0
             && origNorm > VSMALL
             && (status < 0 || orthNorm < relMinSignificant_*origNorm)
            )
            {
                Info<< "Removing IQN-ILS mode (" << patchName
                    << "): ";

                if (status < 0)
                {
                    Info<< "repeated orthogonalization failed";
                }
                else
                {
                    Info<< "|v_orth| = " << orthNorm
                        << " < relMinSignificant*|v| = "
                        << relMinSignificant_*origNorm;
                }

                Info<< endl;

                removeCouplingMode(interfaceI, cols - 1 - i);
                removedMode = true;
                break;
            }

            if (status >= 0)
            {
                Q.append(v);
            }
        }
    }
}


label IQNILSCouplingInterface::orthogonalizeCouplingColumn
(
    vectorField& v,
    scalarField& r,
    scalar& rho,
    const DynamicList<vectorField>& qBasis
) const
{
    static const scalar theta = 1.0/0.7;

    const label colNum = qBasis.size();
    scalar rho0 = Foam::sqrt(sum(v & v));
    rho = rho0;
    r.setSize(colNum + 1);
    r = 0.0;

    label nPasses = 0;

    while (true)
    {
        vectorField u(v.size(), Zero);
        scalarField s(colNum, 0.0);

        for (label j = 0; j < colNum; j++)
        {
            s[j] = sum(qBasis[j] & v);
            u += qBasis[j]*s[j];
        }

        for (label j = 0; j < colNum; j++)
        {
            r[j] += s[j];
        }

        v -= u;
        rho = Foam::sqrt(sum(v & v));
        const scalar normCoefficients =
            Foam::sqrt(Foam::sum(sqr(s)));

        nPasses++;

        if (rho <= std::numeric_limits<scalar>::min())
        {
            return -1;
        }

        if (rho*theta <= rho0 + normCoefficients)
        {
            if (nPasses >= 4)
            {
                return -1;
            }

            rho0 = rho;
        }
        else
        {
            break;
        }
    }

    v /= rho;
    r[colNum] = rho;

    return nPasses;
}


label IQNILSCouplingInterface::buildPreciceStyleCouplingQR
(
    const label interfaceI,
    DynamicList<vectorField>& solveQ,
    DynamicList<vectorField>& solveW,
    RectangularMatrix<scalar>& R
) const
{
    const label storedCols = nSolveCouplingModes(interfaceI);
    const word patchName =
        fluidMesh().boundary()
        [
            fluid().globalPatches()[interfaceI].patch().index()
        ].name();

    solveQ.clear();
    solveW.clear();
    R = RectangularMatrix<scalar>(storedCols, storedCols, 0.0);

    for (label i = 0; i < storedCols; i++)
    {
        vectorField v;
        vectorField w;
        solveCouplingMode(interfaceI, i, v, w);
        applyResidualPreconditioning(interfaceI, v);
        const scalar origNorm = Foam::sqrt(sum(v & v));

        if (normalizeCouplingColumns_ && origNorm > VSMALL)
        {
            v /= origNorm;
            w /= origNorm;
        }

        scalarField r;
        scalar rho = 0.0;
        const label status =
            orthogonalizeCouplingColumn(v, r, rho, solveQ);

        if
        (
            i > 0
         && qrFilterTolerance_ > 0.0
         && origNorm > VSMALL
         && rho < qrFilterTolerance_*origNorm
        )
        {
            Info<< "Ignoring IQN-ILS QR2 column (" << patchName
                << "): |v_orth| = " << rho
                << " < qrFilterTolerance*|v| = "
                << qrFilterTolerance_*origNorm << endl;
            continue;
        }

        if (status < 0)
        {
            Info<< "Ignoring IQN-ILS QR2 column (" << patchName
                << "): repeated orthogonalization failed" << endl;
            continue;
        }

        const label col = solveQ.size();
        solveQ.append(v);
        solveW.append(w);

        for (label row = 0; row <= col; row++)
        {
            R[row][col] = r[row];
        }
    }

    return solveQ.size();
}


label IQNILSCouplingInterface::buildPreciceStyleCombinedCouplingQR
(
    DynamicList<vectorField>& solveQ,
    DynamicList<vectorField>& solveW,
    RectangularMatrix<scalar>& R
)
{
    const label storedCols = nCombinedCouplingModes();

    solveQ.clear();
    solveW.clear();
    R = RectangularMatrix<scalar>(storedCols, storedCols, 0.0);

    for (label i = 0; i < storedCols; i++)
    {
        vectorField v;
        vectorField w;
        combinedCouplingMode(i, v, w);
        applyResidualPreconditioningCombined(v);
        const scalar origNorm = Foam::sqrt(sum(v & v));

        if (normalizeCouplingColumns_ && origNorm > VSMALL)
        {
            v /= origNorm;
            w /= origNorm;
        }

        scalarField r;
        scalar rho = 0.0;
        const label status =
            orthogonalizeCouplingColumn(v, r, rho, solveQ);

        if
        (
            i > 0
         && qrFilterTolerance_ > 0.0
         && origNorm > VSMALL
         && rho < qrFilterTolerance_*origNorm
        )
        {
            Info<< "Ignoring IQN-ILS combined QR2 column: |v_orth| = "
                << rho << " < qrFilterTolerance*|v| = "
                << qrFilterTolerance_*origNorm << endl;
            continue;
        }

        if (status < 0)
        {
            Info<< "Ignoring IQN-ILS combined QR2 column: "
                << "repeated orthogonalization failed" << endl;
            continue;
        }

        const label col = solveQ.size();
        solveQ.append(v);
        solveW.append(w);

        for (label row = 0; row <= col; row++)
        {
            R[row][col] = r[row];
        }
    }

    return solveQ.size();
}


void IQNILSCouplingInterface::filterCouplingModes(const label interfaceI)
{
    const word patchName =
        fluidMesh().boundary()
        [
            fluid().globalPatches()[interfaceI].patch().index()
        ].name();

    while
    (
        fluidPatchesPointsV_[interfaceI].size()
     && (
            3*fluidPatchesPointsV_[interfaceI][0].size()
          < fluidPatchesPointsV_[interfaceI].size()
        )
    )
    {
        Info<< "Removing oldest IQN-ILS mode (" << patchName
            << "): too many columns for available rows" << endl;

        removeCouplingMode(interfaceI, 0);
    }

    filterRelativeCouplingModes(interfaceI);

    if
    (
        minSignificant_ <= 0.0
     || fluidPatchesPointsV_[interfaceI].size() < 2
    )
    {
        return;
    }

    while (fluidPatchesPointsV_[interfaceI].size() > 1)
    {
        const label cols = fluidPatchesPointsV_[interfaceI].size();
        DynamicList<vectorField> Q(cols);
        scalar minDiag = GREAT;
        label minDiagI = -1;

        for (label i = 0; i < cols; i++)
        {
            Q.append(fluidPatchesPointsV_[interfaceI][cols - 1 - i]);
        }

        for (label i = 0; i < cols; i++)
        {
            const scalar diag = Foam::sqrt(sum(Q[i] & Q[i]));

            if (diag < minDiag)
            {
                minDiag = diag;
                minDiagI = i;
            }

            if (diag > VSMALL)
            {
                Q[i] /= diag;
            }

            for (label j = i + 1; j < cols; j++)
            {
                Q[j] -= sum(Q[i] & Q[j])*Q[i];
            }
        }

        if (minDiag >= minSignificant_)
        {
            break;
        }

        Info<< "Removing IQN-ILS mode (" << patchName
            << "): |diag(R)| = " << minDiag
            << " < minSignificant = " << minSignificant_ << endl;

        removeCouplingMode(interfaceI, cols - 1 - minDiagI);
    }
}


void IQNILSCouplingInterface::updateResidualPreconditioning()
{
    if (!residualSumPreconditioning_)
    {
        return;
    }

    forAll(fluid().globalPatches(), interfaceI)
    {
        const vectorField& res = residuals()[interfaceI];
        scalarField& scale = residualPointScale_[interfaceI];
        scale.setSize(res.size());

        forAll(res, pointI)
        {
            const scalar magnitude = Foam::mag(res[pointI]);
            scale[pointI] =
                1.0/(residualPreconditioningEpsilon_ + magnitude);
        }
    }
}


void IQNILSCouplingInterface::applyResidualPreconditioning
(
    const label interfaceI,
    vectorField& v
) const
{
    if (!residualSumPreconditioning_)
    {
        return;
    }

    const scalarField& scale = residualPointScale_[interfaceI];

    if (scale.size() != v.size())
    {
        return;
    }

    forAll(v, pointI)
    {
        v[pointI] *= scale[pointI];
    }
}


void IQNILSCouplingInterface::applyResidualPreconditioningCombined
(vectorField& v) const
{
    if (!residualSumPreconditioning_)
    {
        return;
    }

    label offset = 0;

    forAll(fluid().globalPatches(), interfaceI)
    {
        const scalarField& scale = residualPointScale_[interfaceI];

        forAll(scale, pointI)
        {
            const label index = offset + pointI;

            if (index < v.size())
            {
                v[index] *= scale[pointI];
            }
        }

        offset += scale.size();
    }
}


bool IQNILSCouplingInterface::usePreciceStyleCouplingQR
(
    const label interfaceI
) const
{
    return
        preciceStyleCouplingQR_
     && fluidPatchesPointsV_[interfaceI].size() > 0;
}


bool IQNILSCouplingInterface::useCombinedPreciceStyleCouplingQR() const
{
    if (!preciceStyleCouplingQR_)
    {
        return false;
    }

    forAll(fluid().globalPatches(), interfaceI)
    {
        if (fluidPatchesPointsV_[interfaceI].size() > 0)
        {
            return true;
        }
    }

    return false;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool IQNILSCouplingInterface::evolve()
{
    initializeFields();

    updateInterpolatorAndGlobalPatches();

    scalar residualNorm = 0;

    // Check if coupling switch needs to be updated
    if (!coupled())
    {
        updateCoupled();
    }

    if (predictSolid_ && coupled())
    {
        updateForce();
        solid().evolve();

        residualNorm =
            updateResidual();

        updateResidualPreconditioning();
    }

    do
    {
        outerCorr()++;

        // Transfer the displacement from the solid to the fluid
        updateDisplacement();

        // Move the fluid mesh
        moveFluidMesh();

        // Solve fluid
        fluid().evolve();

        if (coupled())
        {
            // Transfer the force from the fluid to the solid
            updateForce();

            // Solve solid
            solid().evolve();

            // Calculate the FSI residual
            residualNorm = updateResidual();

            updateResidualPreconditioning();
        }
        else
        {
            residualNorm = 0.0;
        }

        // Optional: write residuals to file
        if (writeResidualsToFile() && Pstream::master())
        {
            residualFile()
                << runTime().value() << " "
                << outerCorr() << " "
                << residualNorm << endl;
        }
    }
    while (residualNorm > outerCorrTolerance() && outerCorr() < nOuterCorr());

    cacheCurrentStepModes();

    solid().updateTotalFields();

    // Optional: correct fluid mesh to avoid build-up of interface position
    // errors
    if (additionalMeshCorrection())
    {
        // Transfer the displacement from the solid to the fluid, where we will
        // use no relaxation; in that way, we can force the solid and fluid
        // interfaces to stay aligned
        forAll(fluid().globalPatches(), interfaceI)
        {
            fluidZonesPointsDisplsPrev()[interfaceI] =
                fluidZonesPointsDispls()[interfaceI];

            fluidZonesPointsDispls()[interfaceI] += residuals()[interfaceI];
        }

        // Move the fluid mesh
        moveFluidMesh();
    }

    return 0;
}


void IQNILSCouplingInterface::updateDisplacement()
{
    Info<< nl << "Time = " << fluid().runTime().timeName()
        << ", iteration: " << outerCorr() << endl;

    if (outerCorr() == 1)
    {
        forAll(fluid().globalPatches(), interfaceI)
        {
            Info<< "Modes before clean-up ("
                << fluidMesh().boundary()
                   [
                       fluid().globalPatches()[interfaceI].patch().index()
                   ].name()
                << "): " << fluidPatchesPointsT_[interfaceI].size();

            removeOldCouplingModes
            (
                interfaceI,
                fluid().runTime().timeIndex() - couplingReuse()
            );

            Info<< ", modes after clean-up ("
                << fluidMesh().boundary()
                   [
                       fluid().globalPatches()[interfaceI].patch().index()
                   ].name()
                << "): " << fluidPatchesPointsT_[interfaceI].size();

            if (usePreviousStepModesInSolve(interfaceI))
            {
                Info<< ", cached previous-step modes available: "
                    << previousStepPatchesPointsV_[interfaceI].size();
            }
            else if (usePreviousNonEmptyStepModesInSolve(interfaceI))
            {
                Info<< ", cached previous non-empty step modes available: "
                    << previousNonEmptyStepPatchesPointsV_[interfaceI].size();
            }

            Info<< endl;
        }
    }
    else if (outerCorr() == 2)
    {
        // Set reference in the first coupling iteration
        forAll(fluid().globalPatches(), interfaceI)
        {
            solidZonesPointsDisplsRef()[interfaceI] =
                solidZonesPointsDispls()[interfaceI];

            fluidZonesPointsDisplsRef()[interfaceI] =
                fluidZonesPointsDispls()[interfaceI];

            solidZonesTractionsRef()[interfaceI] =
                solidZonesTractions()[interfaceI];
        }
    }
    else
    {
        forAll(fluid().globalPatches(), interfaceI)
        {
            // Use consecutive iteration differences for the IQN-ILS secants.
            fluidPatchesPointsV_[interfaceI].append
            (
                vectorField
                (
                    (
                        solidZonesPointsDispls()[interfaceI]
                      - fluidZonesPointsDispls()[interfaceI]
                    )
                  - (
                        solidZonesPointsDisplsRef()[interfaceI]
                      - fluidZonesPointsDisplsRef()[interfaceI]
                    )
                )
            );

            fluidPatchesPointsW_[interfaceI].append
            (
                vectorField
                (
                    solidZonesPointsDispls()[interfaceI]
                  - solidZonesPointsDisplsRef()[interfaceI]
                )
            );

            solidPatchesFacesTractionW_[interfaceI].append
            (
                vectorField
                (
                    solidZonesTractions()[interfaceI]
                  - solidZonesTractionsRef()[interfaceI]
                )
            );

            fluidPatchesPointsT_[interfaceI].append
            (
                fluid().runTime().timeIndex()
            );

            solidZonesPointsDisplsRef()[interfaceI] =
                solidZonesPointsDispls()[interfaceI];

            fluidZonesPointsDisplsRef()[interfaceI] =
                fluidZonesPointsDispls()[interfaceI];

            solidZonesTractionsRef()[interfaceI] =
                solidZonesTractions()[interfaceI];
        }
    }


    forAll(fluid().globalPatches(), interfaceI)
    {
        filterCouplingModes(interfaceI);
    }

    if (combinedCouplingSystem_)
    {
        const label storedCols = nCombinedCouplingModes();

        if (storedCols > 0)
        {
            DynamicList<vectorField> solveQ(storedCols);
            DynamicList<vectorField> solveW(storedCols);
            DynamicList<vectorField> qBasis(storedCols);
            vectorField solveResidual(combinedResidual());
            RectangularMatrix<scalar> preciceR(0, 0, 0.0);
            label cols = 0;

            if (useCombinedPreciceStyleCouplingQR())
            {
                cols =
                    buildPreciceStyleCombinedCouplingQR
                    (
                        solveQ,
                        solveW,
                        preciceR
                    );
            }
            else
            {
                for (label i = 0; i < storedCols; i++)
                {
                    vectorField v;
                    vectorField w;
                    combinedCouplingMode(i, v, w);
                    const scalar origNorm = Foam::sqrt(sum(v & v));

                    if (normalizeCouplingColumns_ && origNorm > VSMALL)
                    {
                        v /= origNorm;
                        w /= origNorm;
                    }

                    vectorField vOrth(v);
                    forAll(qBasis, qI)
                    {
                        vOrth -= sum(qBasis[qI] & vOrth)*qBasis[qI];
                    }

                    const scalar orthNorm = Foam::sqrt(sum(vOrth & vOrth));

                    if
                    (
                        qrFilterTolerance_ > 0.0
                     && i > 0
                     && origNorm > VSMALL
                     && orthNorm < qrFilterTolerance_*origNorm
                    )
                    {
                        Info<< "Ignoring IQN-ILS combined solve column: |v_orth| = "
                            << orthNorm << " < qrFilterTolerance*|v| = "
                            << qrFilterTolerance_*origNorm << endl;
                        continue;
                    }

                    solveQ.append(v);
                    solveW.append(w);

                    if (orthNorm > VSMALL)
                    {
                        vOrth /= orthNorm;
                    }

                    qBasis.append(vOrth);
                }

                cols = solveQ.size();
            }

            if (cols > 0)
            {
                Field<scalar> coeffs(cols, 0.0);
                RectangularMatrix<scalar> R(cols, cols, 0.0);
                RectangularMatrix<scalar> Rcolsum(1, cols);

                if (useCombinedPreciceStyleCouplingQR())
                {
                    for (label i = 0; i < cols; i++)
                    {
                        coeffs[i] = -sum(solveQ[i] & solveResidual);
                    }

                    for (label i = 0; i < cols; i++)
                    {
                        for (label j = i; j < cols; j++)
                        {
                            R[i][j] = preciceR[i][j];
                        }
                    }
                }
                else
                {
                    RectangularMatrix<scalar> C(cols, 1);

                    for (label i = 0; i < cols; i++)
                    {
                        R[i][i] = Foam::sqrt(sum(solveQ[i] & solveQ[i]));
                        solveQ[i] /= R[i][i];

                        for (label j = i + 1; j < cols; j++)
                        {
                            R[i][j] = sum(solveQ[i] & solveQ[j]);
                            solveQ[j] -= R[i][j]*solveQ[i];

                            if (reorthogonalizeCouplingColumns_)
                            {
                                const scalar correction =
                                    sum(solveQ[i] & solveQ[j]);
                                R[i][j] += correction;
                                solveQ[j] -= correction*solveQ[i];
                            }
                        }

                        C[i][0] = sum(solveQ[i] & solveResidual);
                    }

                    for (label i = 0; i < cols; i++)
                    {
                        coeffs[cols - 1 - i] = C[cols - 1 - i][0];
                    }
                }

                for (label j = 0; j < cols; j++)
                {
                    Rcolsum[0][j] = 0.0;

                    for (label i = 0; i < j + 1; i++)
                    {
                        Rcolsum[0][j] += cmptMag(R[i][j]);
                    }
                }

                scalar epsilon = qrSolveTolerance_*max(Rcolsum);

                for (label i = 0; i < cols; i++)
                {
                    if (cmptMag(R[i][i]) > epsilon)
                    {
                        for (label j = i + 1; j < cols; j++)
                        {
                            R[i][j] /= R[i][i];
                        }

                        coeffs[i] /= R[i][i];
                        R[i][i] = 1.0;
                    }
                }

                for (label j = cols - 1; j >= 0; j--)
                {
                    if (cmptMag(R[j][j]) > epsilon)
                    {
                        for (label i = 0; i < j; i++)
                        {
                            coeffs[i] -= coeffs[j]*R[i][j];
                        }
                    }
                    else
                    {
                        coeffs[j] = 0.0;
                    }
                }

                vectorField couplingCorrection
                (
                    combinedCouplingSystemSize(),
                    Zero
                );

                for (label i = 0; i < cols; i++)
                {
                    couplingCorrection += solveW[i]*coeffs[i];
                }

                bool reusedPreviousStepModes = false;

                forAll(fluid().globalPatches(), interfaceI)
                {
                    reusedPreviousStepModes =
                        reusedPreviousStepModes
                     || usePreviousStepModesInSolve(interfaceI);
                }

                if
                (
                    reusedPreviousStepModes
                 && maxReuseUpdateNormRatio_ > 0.0
                )
                {
                    const scalar residualNorm =
                        Foam::sqrt(sum(solveResidual & solveResidual));
                    const scalar correctionNorm =
                        Foam::sqrt(sum(couplingCorrection & couplingCorrection));
                    const scalar maxCorrectionNorm =
                        maxReuseUpdateNormRatio_*residualNorm;

                    if
                    (
                        residualNorm > VSMALL
                     && correctionNorm > maxCorrectionNorm
                    )
                    {
                        const scalar scale =
                            maxCorrectionNorm/max(correctionNorm, VSMALL);

                        Info<< "Limiting reused combined IQN-ILS correction: "
                            << "|du| = " << correctionNorm
                            << " > maxReuseUpdateNormRatio*|r| = "
                            << maxCorrectionNorm
                            << ", scale = " << scale << endl;

                        couplingCorrection *= scale;

                    }
                }

                applyCombinedCorrection(couplingCorrection);

                return;
            }
        }

        forAll(fluid().globalPatches(), interfaceI)
        {
            Info<< "Current fsi under-relaxation factor ("
                << fluidMesh().boundary()
                   [
                       fluid().globalPatches()[interfaceI].patch().index()
                   ].name()
                << "): " << relaxationFactor_ << endl;

            fluidZonesPointsDisplsPrev()[interfaceI] =
                fluidZonesPointsDispls()[interfaceI];

            if ((outerCorr() == 1) && predictor())
            {
                fluidZonesPointsDispls()[interfaceI] += residuals()[interfaceI];
            }
            else
            {
                fluidZonesPointsDispls()[interfaceI] +=
                    relaxationFactor_*residuals()[interfaceI];
            }
        }

        return;
    }

    forAll(fluid().globalPatches(), interfaceI)
    {
        const label storedCols = nSolveCouplingModes(interfaceI);

        if (storedCols > 0)
        {
            // Previoulsy given in the function:
            // updateDisplacementUsingIQNILS();

            // Consider fluidPatchesPointsV as a matrix V
            // with as columns the items
            // in the DynamicList and calculate the QR-decomposition of V
            // with modified Gram-Schmidt
            DynamicList<vectorField> solveQ(storedCols);
            DynamicList<vectorField> solveW(storedCols);
            DynamicList<vectorField> qBasis(storedCols);
            const word patchName =
                fluidMesh().boundary()
                [
                    fluid().globalPatches()[interfaceI].patch().index()
                ].name();
            vectorField solveResidual
            (
                fluidZonesPointsDispls()[interfaceI]
              - solidZonesPointsDispls()[interfaceI]
            );
            RectangularMatrix<scalar> preciceR(0, 0, 0.0);
            label cols = 0;

            if (usePreciceStyleCouplingQR(interfaceI))
            {
                cols =
                    buildPreciceStyleCouplingQR
                    (
                        interfaceI,
                        solveQ,
                        solveW,
                        preciceR
                    );
            }
            else
            {
                for (label i = 0; i < storedCols; i++)
                {
                    vectorField v;
                    vectorField w;
                    solveCouplingMode(interfaceI, i, v, w);
                    const scalar origNorm = Foam::sqrt(sum(v & v));

                    if (normalizeCouplingColumns_ && origNorm > VSMALL)
                    {
                        v /= origNorm;
                        w /= origNorm;
                    }

                    vectorField vOrth(v);
                    forAll(qBasis, qI)
                    {
                        vOrth -= sum(qBasis[qI] & vOrth)*qBasis[qI];
                    }

                    const scalar orthNorm = Foam::sqrt(sum(vOrth & vOrth));

                    if
                    (
                        qrFilterTolerance_ > 0.0
                     && i > 0
                     && origNorm > VSMALL
                     && orthNorm < qrFilterTolerance_*origNorm
                    )
                    {
                        Info<< "Ignoring IQN-ILS solve column (" << patchName
                            << "): |v_orth| = " << orthNorm
                            << " < qrFilterTolerance*|v| = "
                            << qrFilterTolerance_*origNorm << endl;
                        continue;
                    }

                    solveQ.append(v);
                    solveW.append(w);

                    if (orthNorm > VSMALL)
                    {
                        vOrth /= orthNorm;
                    }

                    qBasis.append(vOrth);
                }

                cols = solveQ.size();
            }

            if (cols > 0)
            {
                Field<scalar> coeffs(cols, 0.0);
                RectangularMatrix<scalar> R(cols, cols, 0.0);
                RectangularMatrix<scalar> Rcolsum(1, cols);

                if (usePreciceStyleCouplingQR(interfaceI))
                {
                    for (label i = 0; i < cols; i++)
                    {
                        coeffs[i] = -sum(solveQ[i] & solveResidual);
                    }

                    for (label i = 0; i < cols; i++)
                    {
                        for (label j = i; j < cols; j++)
                        {
                            R[i][j] = preciceR[i][j];
                        }
                    }
                }
                else
                {
                    RectangularMatrix<scalar> C(cols, 1);

                    for (label i = 0; i < cols; i++)
                    {
                        // Normalize column i
                        R[i][i] = Foam::sqrt(sum(solveQ[i] & solveQ[i]));
                        solveQ[i] /= R[i][i];

                        // Orthogonalize columns to the right of column i
                        for (label j = i+1; j < cols; j++)
                        {
                            R[i][j] = sum(solveQ[i] & solveQ[j]);
                            solveQ[j] -= R[i][j]*solveQ[i];

                            if (reorthogonalizeCouplingColumns_)
                            {
                                const scalar correction = sum(solveQ[i] & solveQ[j]);
                                R[i][j] += correction;
                                solveQ[j] -= correction*solveQ[i];
                            }
                        }

                        C[i][0] = sum
                            (
                                solveQ[i]
                              & solveResidual
                            );
                    }

                    for (label i = 0; i < cols; i++)
                    {
                        coeffs[cols-1-i] = C[cols-1-i][0];
                    }
                }

                for (label j = 0; j < cols; j++)
                {
                    Rcolsum[0][j] = 0.0;

                    for (label i = 0; i < j+1; i++)
                    {
                        Rcolsum[0][j] += cmptMag(R[i][j]);
                    }
                }

                scalar epsilon = qrSolveTolerance_*max(Rcolsum);

                for (label i = 0; i < cols; i++)
                {
                    if (cmptMag(R[i][i]) > epsilon)
                    {
                        for (label j = i + 1; j < cols; j++)
                        {
                            R[i][j] /= R[i][i];
                        }

                        coeffs[i] /= R[i][i];
                        R[i][i] = 1.0;
                    }
                }

                for (label j = cols-1; j >= 0; j--)
                {
                    if (cmptMag(R[j][j]) > epsilon)
                    {
                        for (label i = 0; i < j; i++)
                        {
                            coeffs[i] -= coeffs[j]*R[i][j];
                        }
                    }
                    else
                    {
                        coeffs[j] = 0.0;
                    }
                }

                fluidZonesPointsDisplsPrev()[interfaceI] =
                    fluidZonesPointsDispls()[interfaceI];

                vectorField couplingCorrection
                (
                    fluidZonesPointsDispls()[interfaceI].size(),
                    Zero
                );

                for (label i = 0; i < cols; i++)
                {
                    couplingCorrection += solveW[i]*coeffs[i];
                }

                if
                (
                    usePreviousStepModesInSolve(interfaceI)
                 && maxReuseUpdateNormRatio_ > 0.0
                )
                {
                    const scalar residualNorm =
                        Foam::sqrt(sum(solveResidual & solveResidual));
                    const scalar correctionNorm =
                        Foam::sqrt(sum(couplingCorrection & couplingCorrection));
                    const scalar maxCorrectionNorm =
                        maxReuseUpdateNormRatio_*residualNorm;

                    if
                    (
                        residualNorm > VSMALL
                     && correctionNorm > maxCorrectionNorm
                    )
                    {
                        const scalar scale =
                            maxCorrectionNorm/max(correctionNorm, VSMALL);

                        Info<< "Limiting reused IQN-ILS correction ("
                            << patchName << "): |du| = " << correctionNorm
                            << " > maxReuseUpdateNormRatio*|r| = "
                            << maxCorrectionNorm
                            << ", scale = " << scale << endl;

                        couplingCorrection *= scale;
                    }
                }

                fluidZonesPointsDispls()[interfaceI] =
                    solidZonesPointsDispls()[interfaceI];

                fluidZonesPointsDispls()[interfaceI] += couplingCorrection;
            }
            else
            {
                Info<< "Current fsi under-relaxation factor ("
                    << patchName << "): " << relaxationFactor_ << endl;

                fluidZonesPointsDisplsPrev()[interfaceI] =
                    fluidZonesPointsDispls()[interfaceI];

                if ((outerCorr() == 1) && predictor())
                {
                    fluidZonesPointsDispls()[interfaceI] += residuals()[interfaceI];
                }
                else
                {
                    fluidZonesPointsDispls()[interfaceI] +=
                        relaxationFactor_*residuals()[interfaceI];
                }
            }
        }
        else
        {
            // Relax the interface displacement
            Info<< "Current fsi under-relaxation factor ("
                << fluidMesh().boundary()
                   [
                       fluid().globalPatches()[interfaceI].patch().index()
                   ].name()
                << "): " << relaxationFactor_ << endl;

            fluidZonesPointsDisplsPrev()[interfaceI] =
                fluidZonesPointsDispls()[interfaceI];

            if ((outerCorr() == 1) && predictor())
            {
                fluidZonesPointsDispls()[interfaceI] += residuals()[interfaceI];
            }
            else
            {
                fluidZonesPointsDispls()[interfaceI] +=
                    relaxationFactor_*residuals()[interfaceI];
            }
        }
    }

    // Update movingWallPressure boundary conditions, if found
    fluidSolidInterface::updateMovingWallPressureAcceleration();

    // Update elasticWallPressure boundary conditions, if found
    fluidSolidInterface::updateElasticWallPressureAcceleration();

    // Make sure that displacement on all processors is equal to one
    // calculated on master processor
    fluidSolidInterface::syncFluidZonePointsDispl(fluidZonesPointsDispls());
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fluidSolidInterfaces

} // End namespace Foam

// ************************************************************************* //
