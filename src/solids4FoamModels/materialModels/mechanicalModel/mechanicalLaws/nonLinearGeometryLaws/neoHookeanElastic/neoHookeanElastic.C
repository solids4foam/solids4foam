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

#include "neoHookeanElastic.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(neoHookeanElastic, 0);
    addToRunTimeSelectionTable
    (
        mechanicalLaw, neoHookeanElastic, nonLinGeomMechLaw
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from dictionary
Foam::neoHookeanElastic::neoHookeanElastic
(
    const word& name,
    const fvMesh& mesh,
    const dictionary& dict,
    const nonLinearGeometry::nonLinearType& nonLinGeom
)
:
    mechanicalLaw(name, mesh, dict, nonLinGeom),
    mu_("mu", dimPressure, 0.0),
    K_("K", dimPressure, 0.0),
    alternatePressureDefinition_
    (
        dict.lookupOrDefault<Switch>("alternatePressureDefinition", false)
    ),
    pressureDisplacement_
    (
        dict.lookupOrDefault<Switch>("pressureDisplacement", false)
    ),
    pressureDisplacementCoeff_
    (
        "pressureDisplacementCoeff",
        dimless,
        dict.lookupOrDefault<scalar>("pressureDisplacementCoeff", -1.0)
    )
{
    // Read mechanical properties
    if
    (
        dict.found("E") && dict.found("nu")
     && !dict.found("mu") && !dict.found("K")
    )
    {
        const dimensionedScalar E = dimensionedScalar(dict.lookup("E"));
        const dimensionedScalar nu = dimensionedScalar(dict.lookup("nu"));

        mu_ = (E/(2.0*(1.0 + nu)));

        if (planeStress())
        {
            K_ = (nu*E/((1.0 + nu)*(1.0 - nu))) + (2.0/3.0)*mu_;
        }
        else if (pressureDisplacement_ && nu.value() > 0.5 - SMALL)
        {
            K_.value() = GREAT;
        }
        else
        {
            K_ = (nu*E/((1.0 + nu)*(1.0 - 2.0*nu))) + (2.0/3.0)*mu_;
        }
    }
    else if
    (
        dict.found("mu") && dict.found("K")
     && !dict.found("E") && !dict.found("nu")
    )
    {
        mu_ = dimensionedScalar(dict.lookup("mu"));
        K_ = dimensionedScalar(dict.lookup("K"));
    }
    else
    {
        FatalErrorIn(type())
            << "Either E and nu or mu and K should be specified"
            << abort(FatalError);
    }

    if (pressureDisplacement_ && pressureDisplacementCoeff_.value() < 0.0)
    {
        dimensionedScalar nu("nu", dimless, 0.5);

        if (K_.value() < GREAT - SMALL)
        {
            nu =
                (3.0*K_ - 2.0*mu_)
               /(2.0*(3.0*K_ + mu_));
        }

        pressureDisplacementCoeff_ = 3.0*nu/(1.0 + nu);
    }

    // Store old F
    F().storeOldTime();
    Ff().storeOldTime();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::neoHookeanElastic::~neoHookeanElastic()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::neoHookeanElastic::impK() const
{
    if (pressureDisplacement_)
    {
        return tmp<volScalarField>
        (
            new volScalarField
            (
                IOobject
                (
                    "impK",
                    mesh().time().timeName(),
                    mesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh(),
                mu_
            )
        );
    }

    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "impK",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            (4.0/3.0)*mu_ + K_ // == 2*mu + lambda
        )
    );
}


Foam::tmp<Foam::volScalarField> Foam::neoHookeanElastic::bulkModulus() const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "K",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            K_
        )
    );
}


Foam::tmp<Foam::volScalarField> Foam::neoHookeanElastic::shearModulus() const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "mu",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            mu_
        )
    );
}


void Foam::neoHookeanElastic::materialTangentField(List<mat66>& matTan) const
{
    // Set the list size
    matTan.resize(mesh().nFaces());

    // Calculate tangent field
    {
        // Lookup current stress and store it as the reference
        const surfaceSymmTensorField& sigmaRef =
            mesh().lookupObject<surfaceSymmTensorField>("sigmaf");

        // Lookup gradient of displacement
        const surfaceTensorField& gradDRef =
            mesh().lookupObject<surfaceTensorField>("grad(D)f");

        // Create fields to be used for perturbations
        surfaceSymmTensorField sigmaPerturb("sigmaPerturb", sigmaRef);
        surfaceTensorField gradDPerturb("gradDPerturb", gradDRef);

        // Small number used for perturbations
        const scalar eps(readScalar(dict().lookup("tangentEps")));

        // Define matrix indices for readability
        const label XX = symmTensor::XX;
        const label YY = symmTensor::YY;
        const label ZZ = symmTensor::ZZ;
        const label XY = symmTensor::XY;
        const label YZ = symmTensor::YZ;
        const label XZ = symmTensor::XZ;

        // For each component of gradD, sequentially apply a perturbation and
        // then calculate the resulting sigma
        for (label cmptI = 0; cmptI < symmTensor::nComponents; cmptI++)
        {
            // Map tensor component to symmTensor
            // We can avoid this is we perturb epsilon directly
            label tensorCmptI = -1;
            if (cmptI == symmTensor::XX)
            {
                tensorCmptI = tensor::XX;
            }
            else if (cmptI == symmTensor::XY)
            {
                tensorCmptI = tensor::XY;
            }
            else if (cmptI == symmTensor::XZ)
            {
                tensorCmptI = tensor::XZ;
            }
            else if (cmptI == symmTensor::YY)
            {
                tensorCmptI = tensor::YY;
            }
            else if (cmptI == symmTensor::YZ)
            {
                tensorCmptI = tensor::YZ;
            }
            else // if (cmptI == symmTensor::ZZ)
            {
                tensorCmptI = tensor::ZZ;
            }

            // Reset gradDPerturb and multiply by 1.0 to avoid it being removed
            // from the object registry
            gradDPerturb = 1.0*gradDRef;

            // Perturb this component of gradD and calculate FPerturb
            gradDPerturb.replace
            (
                tensorCmptI, gradDRef.component(tensorCmptI) + eps
            );

            // Calculate perturbed stress
            const_cast<neoHookeanElastic&>(*this).correct
            (
                sigmaPerturb, gradDPerturb
            );

            // Calculate tangent component
            const surfaceSymmTensorField tangCmpt((sigmaPerturb - sigmaRef)/eps);
            const symmTensorField& tangCmptI = tangCmpt.internalField();

            // Insert tangent component
            forAll(tangCmptI, faceI)
            {
                // Take a reference to the current tangent
                mat66& curMatTan = matTan[faceI];

                // Zero the tangent
                curMatTan.clear();

                curMatTan(XX, cmptI) = tangCmptI[faceI][XX];
                curMatTan(YY, cmptI) = tangCmptI[faceI][YY];
                curMatTan(ZZ, cmptI) = tangCmptI[faceI][ZZ];
                curMatTan(XY, cmptI) = tangCmptI[faceI][XY];
                curMatTan(YZ, cmptI) = tangCmptI[faceI][YZ];
                curMatTan(XZ, cmptI) = tangCmptI[faceI][XZ];
            }

            forAll(tangCmpt.boundaryField(), patchI)
            {
                const symmTensorField& tangCmptP =
                    tangCmpt.boundaryField()[patchI];
                const label start = mesh().boundaryMesh()[patchI].start();

                forAll(tangCmptP, fI)
                {
                    const label faceID = start + fI;

                    // Take a reference to the current tangent
                    mat66& curMatTan = matTan[faceID];

                    // Zero the tangent
                    curMatTan.clear();

                    curMatTan(XX, cmptI) = tangCmptP[fI][XX];
                    curMatTan(YY, cmptI) = tangCmptP[fI][YY];
                    curMatTan(ZZ, cmptI) = tangCmptP[fI][ZZ];
                    curMatTan(XY, cmptI) = tangCmptP[fI][XY];
                    curMatTan(YZ, cmptI) = tangCmptP[fI][YZ];
                    curMatTan(XZ, cmptI) = tangCmptP[fI][XZ];
                }
            }
        }
    }
}


void Foam::neoHookeanElastic::correct(volSymmTensorField& sigma)
{
    // Update the deformation gradient field
    // Note: if true is returned, it means that linearised elasticity was
    // enforced by the solver via the enforceLinear switch
    if (updateF(sigma, mu_, K_))
    {
        return;
    }

    // NOTE [IMPORTANT]:
    // Do NOT write F.T() & F directly: see the comment in
    // StVenantKirchhoffElastic.C
    const volTensorField& F = this->F();
    const volTensorField FT(F.T());

    // Calculate the Jacobian of the deformation gradient
    const volScalarField J(det(F));

    // Calculate the left Cauchy Green strain
    const volSymmTensorField b(symm(F & FT));

    // Calculate the volume preserving left Cauchy Green strain
    const volSymmTensorField bEbar(pow(J, -2.0/3.0)*b);

    // Calculate the deviatoric stress
    const volSymmTensorField s
    (
        pressureDisplacement_
      ? mu_*(b - I)/J
      : mu_*dev(bEbar)
    );

    if (pressureDisplacement_)
    {
        const volScalarField& p = mesh().lookupObject<volScalarField>("p");

        sigma = -pressureDisplacementCoeff_*p*I + s;

        return;
    }

    // Update the hydrostatic stress
    if (alternatePressureDefinition_)
    {
        updateSigmaHyd
        (
            K()*(J - 1.0),
            (4.0/3.0)*mu_ + K_
        );
    }
    else
    {
        updateSigmaHyd
        (
            0.5*K()*(pow(J, 2.0) - 1.0),
            (4.0/3.0)*mu_ + K_
        );
    }

    // Calculate the Cauchy stress
    sigma = (1.0/J)*(sigmaHyd()*I + s);
}


void Foam::neoHookeanElastic::correct(surfaceSymmTensorField& sigma)
{
    // Update the deformation gradient field
    // Note: if true is returned, it means that linearised elasticity was
    // enforced by the solver via the enforceLinear switch
    if (updateF(sigma, mu_, K_))
    {
        return;
    }

    // Update stress
    correctF(sigma, Ff());
}


void Foam::neoHookeanElastic::correct
(
    surfaceSymmTensorField& sigma,
    const surfaceTensorField& gradD
) const
{
    // Update deformation gradient
    const surfaceTensorField F(I + gradD.T());

    // Update stress
    correctF(sigma, F);
}


void Foam::neoHookeanElastic::correctF
(
    surfaceSymmTensorField& sigma,
    const surfaceTensorField& F
) const
{
    // NOTE [IMPORTANT]:
    // Do NOT write F.T() & F directly: see the comment in
    // StVenantKirchhoffElastic.C
    const surfaceTensorField FT(F.T());

    // Calculate the Jacobian of the deformation gradient
    const surfaceScalarField J(det(F));

    // Calculate left Cauchy Green strain tensor
    const surfaceSymmTensorField b(symm(F & FT));

    // Calculate left Cauchy Green strain tensor with volumetric term removed
    const surfaceSymmTensorField bEbar(pow(J, -2.0/3.0)*b);

    // Calculate deviatoric stress
    const surfaceSymmTensorField s
    (
        pressureDisplacement_
      ? mu_*(b - I)/J
      : mu_*dev(bEbar)
    );

    if (pressureDisplacement_)
    {
        const surfaceScalarField& pf =
            mesh().lookupObject<surfaceScalarField>("pf");

        sigma = -pressureDisplacementCoeff_*pf*I + s;

        return;
    }

    // Calculate the Cauchy stress
    sigma = (1.0/J)*(0.5*K_*(pow(J, 2) - 1)*I + s);
}


void Foam::neoHookeanElastic::setRestart()
{
    F().writeOpt() = IOobject::AUTO_WRITE;
    Ff().writeOpt() = IOobject::AUTO_WRITE;
}

// ************************************************************************* //
