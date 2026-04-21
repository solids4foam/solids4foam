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

#include "volStrainRateStab.H"
#include "addToRunTimeSelectionTable.H"
#include "fvc.H"
#include "fvmLaplacian.H"
#include "compatibilityFunctions.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(volStrainRateStab, 0);
    addToRunTimeSelectionTable
    (
        stabilisationModel, volStrainRateStab, stabModel
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::volStrainRateStab::volStrainRateStab
(
    const fvMesh& mesh,
    const dictionary& dict,
    const dimensionSet& dims
)
:
    stabilisationModel(mesh, dict, dims),
    mode_(word(dict.lookup("mode"))),
    tau_
    (
        mode_ == "physicalDamping"
      ? readScalar(dict.lookup("tau"))
      : scalar(1.0)
    )
{
    // Validate mode at construction time for an early, clear error message
    if
    (
        mode_ != "physicalDamping"
     && mode_ != "firstOrderTemporal"
     && mode_ != "secondOrderTemporal"
     && mode_ != "spatioTemporal"
     && mode_ != "h2PhysicalDamping"
    )
    {
        FatalIOErrorInFunction(dict)
            << "Unknown mode: " << mode_ << nl
            << "Valid options are: physicalDamping, firstOrderTemporal, "
            << "secondOrderTemporal, spatioTemporal, h2PhysicalDamping"
            << exit(FatalIOError);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::volStrainRateStab::~volStrainRateStab()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::volStrainRateStab::updateVector
(
    const volVectorField& field,
    const volTensorField* gradPtr
) const
{
    clearCellVectorCache();

    if (gradPtr == nullptr)
    {
        FatalErrorInFunction
            << "A non-null gradient pointer (gradPtr) must be supplied to "
            << type() << "::updateVector; it must also carry an old-time "
            << "level accessible via gradPtr->oldTime()."
            << abort(FatalError);
    }

    // Initialise the face stabilisation field on first call
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

    // Face-interpolated volumetric strain at current and previous time levels.
    // tr(gradD) = div(D) is the volumetric strain (dimensionless for
    // displacement fields in SI units).
    const surfaceScalarField trGradDf
    (
        fvc::interpolate(tr(*gradPtr))
    );
    const surfaceScalarField trGradDOldf
    (
        fvc::interpolate(tr(gradPtr->oldTime()))
    );

    // Volumetric strain increment at faces: O(deltaT) for smooth solutions
    const surfaceScalarField trDiff(trGradDf - trGradDOldf);

    // Face unit normals
    surfaceVectorField n(mesh().Sf() / mesh().magSf());

    tmp<surfaceVectorField> trhs;

    if (mode_ == "physicalDamping")
    {
        // faceVector = scaleFactor * (tau/deltaT) * trDiff * n
        // Approximates scaleFactor * tau * d(tr(gradD))/dt * n.
        // Does NOT vanish as deltaT -> 0; represents permanent physical
        // viscous-like damping of volumetric strain oscillations.
        const scalar deltaT = mesh().time().deltaT().value();

        trhs = tmp<surfaceVectorField>
        (
            scaleFactor() * (tau_ / deltaT) * trDiff * n
        );
    }
    else if (mode_ == "firstOrderTemporal")
    {
        // faceVector = scaleFactor * trDiff * n
        // Since trDiff = O(deltaT), this term vanishes as O(deltaT) with
        // time-step refinement.
        trhs = tmp<surfaceVectorField>
        (
            scaleFactor() * trDiff * n
        );
    }
    else if (mode_ == "secondOrderTemporal")
    {
        // faceVector = scaleFactor * deltaT * trDiff * n
        // Since trDiff = O(deltaT), this term vanishes as O(deltaT^2) with
        // time-step refinement.  scaleFactor carries an implicit [1/s].
        const scalar deltaT = mesh().time().deltaT().value();

        trhs = tmp<surfaceVectorField>
        (
            scaleFactor() * deltaT * trDiff * n
        );
    }
    else if (mode_ == "h2PhysicalDamping")
    {
        // faceVector = scaleFactor * (h^2/m^2 / deltaT) * trDiff * n
        // Vanishes O(h^2) with mesh refinement (like Rhie-Chow); does NOT
        // vanish as deltaT -> 0 (like physicalDamping).  Combines spatial
        // consistency with robust crash prevention at small time steps.
        const scalar deltaT = mesh().time().deltaT().value();
        const dimensionedScalar one("one", dimless/sqr(dimLength), 1.0);
        const surfaceScalarField h2f_nd(one * h2());

        trhs = tmp<surfaceVectorField>
        (
            scaleFactor() * (h2f_nd / deltaT) * trDiff * n
        );
    }
    else // spatioTemporal
    {
        // faceVector = scaleFactor * (h^2/m^2) * trDiff * n
        // Since trDiff = O(deltaT), this term vanishes as O(h^2 * deltaT)
        // with mesh and time-step refinement.  scaleFactor carries an
        // implicit [1/m^2].
        const dimensionedScalar one("one", dimless/sqr(dimLength), 1.0);
        const surfaceScalarField h2f_nd(one * h2());

        trhs = tmp<surfaceVectorField>
        (
            scaleFactor() * h2f_nd * trDiff * n
        );
    }

#ifdef OPENFOAM_COM
    tmpRef(trhs).setOriented(true);
#endif

    autoPtrRef(faceVectorPtr()) = tmpRef(trhs);
}


const Foam::fvVectorMatrix& Foam::volStrainRateStab::vectorJacobian
(
    const volVectorField& field,
    const surfaceScalarField* gammaPtr,
    const bool rebuild
) const
{
    if (vectorJacobianPtr().empty() || rebuild)
    {
        // Scheme name for fvSchemes look-up
        const word schemeName
        (
            gammaPtr
          ? "laplacian(" + gammaPtr->name() + "," + field.name() + ")"
          : "laplacian(" + field.name() + ")"
        );

        if (mode_ == "physicalDamping")
        {
            // Jacobian coefficient: scaleFactorJacobian * tau/deltaT
            // Note: rebuild should be triggered by the caller when deltaT
            // changes between time steps.
            const scalar C
            (
                scaleFactorJacobian() * tau_ / mesh().time().deltaT().value()
            );

            if (gammaPtr)
            {
                vectorJacobianPtr().reset
                (
                    new fvVectorMatrix
                    (
                        C * fvm::laplacian(*gammaPtr, field, schemeName)
                    )
                );
            }
            else
            {
                vectorJacobianPtr().reset
                (
                    new fvVectorMatrix
                    (
                        C * fvm::laplacian(field, schemeName)
                    )
                );
            }
        }
        else if (mode_ == "firstOrderTemporal")
        {
            // Jacobian coefficient: scaleFactorJacobian
            if (gammaPtr)
            {
                vectorJacobianPtr().reset
                (
                    new fvVectorMatrix
                    (
                        scaleFactorJacobian()
                       *fvm::laplacian(*gammaPtr, field, schemeName)
                    )
                );
            }
            else
            {
                vectorJacobianPtr().reset
                (
                    new fvVectorMatrix
                    (
                        scaleFactorJacobian()
                       *fvm::laplacian(field, schemeName)
                    )
                );
            }
        }
        else if (mode_ == "secondOrderTemporal")
        {
            // Jacobian coefficient: scaleFactorJacobian * deltaT
            // Note: rebuild should be triggered by the caller when deltaT
            // changes between time steps.
            const scalar C
            (
                scaleFactorJacobian() * mesh().time().deltaT().value()
            );

            if (gammaPtr)
            {
                vectorJacobianPtr().reset
                (
                    new fvVectorMatrix
                    (
                        C * fvm::laplacian(*gammaPtr, field, schemeName)
                    )
                );
            }
            else
            {
                vectorJacobianPtr().reset
                (
                    new fvVectorMatrix
                    (
                        C * fvm::laplacian(field, schemeName)
                    )
                );
            }
        }
        else if (mode_ == "h2PhysicalDamping")
        {
            // Jacobian coefficient: scaleFactorJacobian / deltaT, with h^2/m^2
            // as face diffusivity matching updateVector.
            // Note: rebuild should be triggered by the caller when deltaT
            // changes between time steps.
            const dimensionedScalar one("one", dimless/sqr(dimLength), 1.0);
            const surfaceScalarField h2f_nd(one * h2());
            const scalar C
            (
                scaleFactorJacobian() / mesh().time().deltaT().value()
            );

            if (gammaPtr)
            {
                vectorJacobianPtr().reset
                (
                    new fvVectorMatrix
                    (
                        fvm::laplacian
                        (
                            C * h2f_nd * (*gammaPtr),
                            field,
                            schemeName
                        )
                    )
                );
            }
            else
            {
                vectorJacobianPtr().reset
                (
                    new fvVectorMatrix
                    (
                        C * fvm::laplacian(h2f_nd, field, schemeName)
                    )
                );
            }
        }
        else // spatioTemporal
        {
            // Jacobian uses h^2/m^2 as the face diffusivity, matching the
            // spatial scaling in updateVector
            const dimensionedScalar one("one", dimless/sqr(dimLength), 1.0);
            const surfaceScalarField h2f_nd(one * h2());

            if (gammaPtr)
            {
                vectorJacobianPtr().reset
                (
                    new fvVectorMatrix
                    (
                        scaleFactorJacobian()
                       *fvm::laplacian
                        (
                            h2f_nd * (*gammaPtr),
                            field,
                            schemeName
                        )
                    )
                );
            }
            else
            {
                vectorJacobianPtr().reset
                (
                    new fvVectorMatrix
                    (
                        scaleFactorJacobian()
                       *fvm::laplacian(h2f_nd, field, schemeName)
                    )
                );
            }
        }
    }

    return vectorJacobianPtr();
}


// ************************************************************************* //
