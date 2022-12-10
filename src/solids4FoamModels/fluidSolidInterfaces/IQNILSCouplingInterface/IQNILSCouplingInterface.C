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
        fsiProperties().lookupOrDefault<scalar>("relaxationFactor", 0.01)
    ),
    couplingReuse_(fsiProperties().lookupOrDefault<int>("couplingReuse", 0)),
    predictSolid_(fsiProperties().lookupOrDefault<bool>("predictSolid", true)),
    fluidPatchesPointsV_(nGlobalPatches()),
    fluidPatchesPointsW_(nGlobalPatches()),
    fluidPatchesPointsT_(nGlobalPatches())
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
        // Clean up data from old time steps
        forAll(fluid().globalPatches(), interfaceI)
        {
            Info<< "Modes before clean-up ("
                << fluidMesh().boundary()
                   [
                       fluid().globalPatches()[interfaceI].patch().index()
                   ].name()
                << "): " << fluidPatchesPointsT_[interfaceI].size();

            while (true)
            {
                if (fluidPatchesPointsT_[interfaceI].size())
                {
                    if
                    (
                        (fluid().runTime().timeIndex() - couplingReuse())
                      > fluidPatchesPointsT_[interfaceI][0]
                    )
                    {
                        for
                        (
                            label i = 0;
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
                        }

                        fluidPatchesPointsT_[interfaceI].remove();
                        fluidPatchesPointsV_[interfaceI].remove();
                        fluidPatchesPointsW_[interfaceI].remove();
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

            Info<< ", modes after clean-up ("
                << fluidMesh().boundary()
                   [
                       fluid().globalPatches()[interfaceI].patch().index()
                   ].name()
                << "): " << fluidPatchesPointsT_[interfaceI].size() << endl;
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
        }
    }
    else
    {
        forAll(fluid().globalPatches(), interfaceI)
        {
            // Reference has been set in the first coupling iteration
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

            fluidPatchesPointsT_[interfaceI].append
            (
                fluid().runTime().timeIndex()
            );
        }
    }


    forAll(fluid().globalPatches(), interfaceI)
    {
        if (fluidPatchesPointsT_[interfaceI].size() > 1)
        {
            // Previoulsy given in the function:
            // updateDisplacementUsingIQNILS();

            // Consider fluidPatchesPointsV as a matrix V
            // with as columns the items
            // in the DynamicList and calculate the QR-decomposition of V
            // with modified Gram-Schmidt
            label cols = fluidPatchesPointsV_[interfaceI].size();
            RectangularMatrix<scalar> R(cols, cols, 0.0);
            RectangularMatrix<scalar> C(cols, 1);
            RectangularMatrix<scalar> Rcolsum(1, cols);
            DynamicList<vectorField> Q;

            for (label i = 0; i < cols; i++)
            {
                Q.append(fluidPatchesPointsV_[interfaceI][cols-1-i]);
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
                            fluidZonesPointsDispls()[interfaceI]
                          - solidZonesPointsDispls()[interfaceI]
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

            fluidZonesPointsDisplsPrev()[interfaceI] =
                fluidZonesPointsDispls()[interfaceI];

            fluidZonesPointsDispls()[interfaceI] =
                solidZonesPointsDispls()[interfaceI];

            for (label i = 0; i < cols; i++)
            {
                fluidZonesPointsDispls()[interfaceI] +=
                    fluidPatchesPointsW_[interfaceI][i]*C[cols-1-i][0];
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
