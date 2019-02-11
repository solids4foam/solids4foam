/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "IQNILSCouplingInterface.H"
#include "addToRunTimeSelectionTable.H"
#include "RectangularMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace fluidSolidInterfaces
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(IQNILSCouplingInterface, 0);
addToRunTimeSelectionTable
(
    physicsModel, IQNILSCouplingInterface, fluidSolidInteraction
);
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
        fsiProperties().lookupOrDefault<scalar>("relaxationFactor", 0.01)
    ),
    couplingReuse_(fsiProperties().lookupOrDefault<int>("couplingReuse", 0)),
    predictSolid_(fsiProperties().lookupOrDefault<bool>("predictSolid", true)),
    fluidPatchPointsV_(),
    fluidPatchPointsW_(),
    fluidPatchPointsT_()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool IQNILSCouplingInterface::evolve()
{
    initializeFields();

    updateInterpolatorAndGlobalPatches();

    scalar residualNorm = 0;

    if (predictSolid_)
    {
        updateForce();

        solid().evolve();

        residualNorm =
            updateResidual();
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

        // Transfer the force from the fluid to the solid
        updateForce();

        // Solve solid
        solid().evolve();

        // Calculate the FSI residual
        residualNorm = updateResidual();
    }
    while (residualNorm > outerCorrTolerance() && outerCorr() < nOuterCorr());

    solid().updateTotalFields();

    return 0;
}


void IQNILSCouplingInterface::updateDisplacement()
{
    Info<< nl << "Time = " << fluid().runTime().timeName()
        << ", iteration: " << outerCorr() << endl;

    if (outerCorr() == 1)
    {
        // Clean up data from old time steps

        Info<< "Modes before clean-up : " << fluidPatchPointsT_.size();

        while (true)
        {
            if (fluidPatchPointsT_.size())
            {
                if
                (
                    fluid().runTime().timeIndex()-couplingReuse()
                  > fluidPatchPointsT_[0]
                )
                {
                    for (label i = 0; i < fluidPatchPointsT_.size() - 1; i++)
                    {
                        fluidPatchPointsT_[i] = fluidPatchPointsT_[i + 1];
                        fluidPatchPointsV_[i] = fluidPatchPointsV_[i + 1];
                        fluidPatchPointsW_[i] = fluidPatchPointsW_[i + 1];
                    }

                    fluidPatchPointsT_.remove();
                    fluidPatchPointsV_.remove();
                    fluidPatchPointsW_.remove();
                }
                else
                {
                    break;
                }
            }
            else
            {
                break;
            }
        }

        Info<< ", modes after clean-up : "
            << fluidPatchPointsT_.size() << endl;
    }
    else if (outerCorr() == 2)
    {
        // Set reference in the first coupling iteration
        solidZonePointsDisplRef() = solidZonePointsDispl();
        fluidZonePointsDisplRef() = fluidZonePointsDispl();
    }
    else
    {
        // Reference has been set in the first coupling iteration
        fluidPatchPointsV_.append
        (
            (
                solidZonePointsDispl()
              - fluidZonePointsDispl()
            )
          - (
                solidZonePointsDisplRef()
              - fluidZonePointsDisplRef()
            )
        );

        fluidPatchPointsW_.append
        (
            solidZonePointsDispl()
            - solidZonePointsDisplRef()
        );

        fluidPatchPointsT_.append
        (
            fluid().runTime().timeIndex()
        );
    }

    if (fluidPatchPointsT_.size() > 1)
    {
        // Previoulsy given in the function:
        // updateDisplacementUsingIQNILS();

        // Consider fluidPatchPointsV as a matrix V
        // with as columns the items
        // in the DynamicList and calculate the QR-decomposition of V
        // with modified Gram-Schmidt
        label cols = fluidPatchPointsV_.size();
        RectangularMatrix<scalar> R(cols, cols, 0.0);
        RectangularMatrix<scalar> C(cols, 1);
        RectangularMatrix<scalar> Rcolsum(1, cols);
        DynamicList<vectorField> Q;

        for (label i = 0; i < cols; i++)
        {
            Q.append(fluidPatchPointsV_[cols-1-i]);
        }

        for (label i = 0; i < cols; i++)
        {
            // Normalize column i
            R[i][i] = Foam::sqrt(sum(Q[i] & Q[i]));
            Q[i] /= R[i][i];

            // Orthogonalize columns to the right of column i
            for (label j = i+1; j < cols; j++)
            {
                R[i][j] = sum(Q[i] & Q[j]);
                Q[j] -= R[i][j]*Q[i];
            }

            // Project minus the residual vector on the Q
            C[i][0] = sum
                (
                    Q[i]
                  & (
                        fluidZonePointsDispl()
                      - solidZonePointsDispl()
                    )
                );
        }

        // Solve the upper triangular system
        for (label j = 0; j < cols; j++)
        {
            Rcolsum[0][j] = 0.0;
            for (label i = 0; i < j+1; i++)
            {
                Rcolsum[0][j] += cmptMag(R[i][j]);
            }
        }
        scalar epsilon = 1.0E-10*max(Rcolsum);
        for (label i = 0; i < cols; i++)
        {
            if (cmptMag(R[i][i]) > epsilon)
            {
                for (label j = i + 1; j < cols; j++)
                {
                    R[i][j] /= R[i][i];
                }
                C[i][0] /= R[i][i];
                R[i][i] = 1.0;
            }
        }
        for (label j = cols-1; j >= 0; j--)
        {
            if (cmptMag(R[j][j]) > epsilon)
            {
                for (label i = 0; i < j; i++)
                {
                    C[i][0] -= C[j][0]*R[i][j];
                }
            }
            else
            {
                C[j][0] = 0.0;
            }
        }

        fluidZonePointsDisplPrev() = fluidZonePointsDispl();

        fluidZonePointsDispl() = solidZonePointsDispl();

        for (label i = 0; i < cols; i++)
        {
            fluidZonePointsDispl() += fluidPatchPointsW_[i]*C[cols-1-i][0];
        }
    }
    else
    {
        // Relax the interface displacement
        Info<< "Current fsi under-relaxation factor: "
            << relaxationFactor_ << endl;

        fluidZonePointsDisplPrev() = fluidZonePointsDispl();

        if ((outerCorr() == 1) && predictor())
        {
            fluidZonePointsDispl() += residual();
        }
        else
        {
            fluidZonePointsDispl() += relaxationFactor_*residual();
        }
    }

    // Update movingWallPressure boundary conditions, if found
    fluidSolidInterface::updateMovingWallPressureAcceleration();

    // Make sure that displacement on all processors is equal to one
    // calculated on master processor
    fluidSolidInterface::syncFluidZonePointsDispl(fluidZonePointsDispl());
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fluidSolidInterfaces

} // End namespace Foam

// ************************************************************************* //
