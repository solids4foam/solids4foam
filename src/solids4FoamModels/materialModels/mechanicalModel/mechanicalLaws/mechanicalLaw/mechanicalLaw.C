/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     |
    \\  /    A nd           | For copyright notice see file Copyright
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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
    Foam::mechanicalLaw

\*---------------------------------------------------------------------------*/

#include "mechanicalLaw.H"
#include "volFields.H"
#include "fvc.H"
#include "IOdictionary.H"
#include "lookupSolidModel.H"
#include "solidModel.H"
#include "fvm.H"
#include "fvc.H"
#include "zeroGradientFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(mechanicalLaw, 0);
    defineRunTimeSelectionTable(mechanicalLaw, linGeomMechLaw);
    defineRunTimeSelectionTable(mechanicalLaw, nonLinGeomMechLaw);
}

// * * * * * * * * * * *  Private Member Funtcions * * * * * * * * * * * * * //

void Foam::mechanicalLaw::makeMu() const
{
    if (muPtr_.valid())
    {
        FatalErrorIn("void " + type() + "::makeMu()")
            << "pointer already set" << abort(FatalError);
    }

    muPtr_.set
    (
        new volScalarField
        (
            IOobject
            (
                "mu_" + name_,
                mesh().time().timeName(),
                mesh(),
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
            ),
            mesh(),
            dimensionedScalar("0", dimPressure, 0.0)
        )
    );
}


void Foam::mechanicalLaw::makeMuf() const
{
    if (mufPtr_.valid())
    {
        FatalErrorIn("void " + type() + "::makeMuf()")
            << "pointer already set" << abort(FatalError);
    }

    mufPtr_.set
    (
        new surfaceScalarField
        (
            IOobject
            (
                "muf_" + name_,
                mesh().time().timeName(),
                mesh(),
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
            ),
            mesh(),
            dimensionedScalar("0", dimPressure, 0.0)
        )
    );
}


void Foam::mechanicalLaw::makeK() const
{
    if (KPtr_.valid())
    {
        FatalErrorIn("void " + type() + "::makeK()")
            << "pointer already set" << abort(FatalError);
    }

    KPtr_.set
    (
        new volScalarField
        (
            IOobject
            (
                "K_" + name_,
                mesh().time().timeName(),
                mesh(),
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
            ),
            mesh(),
            dimensionedScalar("0", dimPressure, 0.0)
        )
    );
}


void Foam::mechanicalLaw::makeKf() const
{
    if (KfPtr_.valid())
    {
        FatalErrorIn("void " + type() + "::makeKf()")
            << "pointer already set" << abort(FatalError);
    }

    KfPtr_.set
    (
        new surfaceScalarField
        (
            IOobject
            (
                "Kf_" + name_,
                mesh().time().timeName(),
                mesh(),
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
            ),
            mesh(),
            dimensionedScalar("0", dimPressure, 0.0)
        )
    );
}


void Foam::mechanicalLaw::makeSigma0() const
{
    if (sigma0Ptr_.valid())
    {
        FatalErrorIn("void " + type() + "::makeSigma0()")
            << "pointer already set" << abort(FatalError);
    }

    sigma0Ptr_.set
    (
        new volSymmTensorField
        (
            IOobject
            (
                "sigma0",
                mesh().time().timeName(),
                mesh(),
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
            ),
            mesh(),
            dimensionedSymmTensor("0", dimPressure, symmTensor::zero)
        )
    );
}


void Foam::mechanicalLaw::makeEpsilon() const
{
    if (epsilonPtr_.valid())
    {
        FatalErrorIn("void " + type() + "::makeEpsilon()")
            << "pointer already set" << abort(FatalError);
    }

    epsilonPtr_.set
    (
        new volSymmTensorField
        (
            IOobject
            (
                "epsilon_",
                mesh().time().timeName(),
                mesh(),
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
            ),
            mesh(),
            dimensionedSymmTensor("0", dimless, symmTensor::zero)
        )
    );
}


void Foam::mechanicalLaw::makeEpsilonf() const
{
    if (epsilonfPtr_.valid())
    {
        FatalErrorIn("void " + type() + "::makeEpsilonf()")
            << "pointer already set" << abort(FatalError);
    }

    epsilonfPtr_.set
    (
        new surfaceSymmTensorField
        (
            IOobject
            (
                "epsilonf_",
                mesh().time().timeName(),
                mesh(),
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
            ),
            mesh(),
            dimensionedSymmTensor("0", dimless, symmTensor::zero)
        )
    );
}


void Foam::mechanicalLaw::makeF()
{
    if (FPtr_.valid())
    {
        FatalErrorIn("void " + type() + "::makeF()")
            << "pointer already set" << abort(FatalError);
    }

    FPtr_.set
    (
        new volTensorField
        (
            IOobject
            (
                "F_",
                mesh().time().timeName(),
                mesh(),
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
            ),
            mesh(),
            dimensionedTensor("I", dimless, I)
        )
    );
}


void Foam::mechanicalLaw::makeFf()
{
    if (FfPtr_.valid())
    {
        FatalErrorIn("void " + type() + "::makeFf()")
            << "pointer already set" << abort(FatalError);
    }

    FfPtr_.set
    (
        new surfaceTensorField
        (
            IOobject
            (
                "Ff_",
                mesh().time().timeName(),
                mesh(),
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
            ),
            mesh(),
            dimensionedTensor("I", dimless, I)
        )
    );
}


void Foam::mechanicalLaw::makeRelF()
{
    if (relFPtr_.valid())
    {
        FatalErrorIn("void " + type() + "::makeRelF()")
            << "pointer already set" << abort(FatalError);
    }

    relFPtr_.set
    (
        new volTensorField
        (
            IOobject
            (
                "relF_",
                mesh().time().timeName(),
                mesh(),
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
            ),
            mesh(),
            dimensionedTensor("I", dimless, I)
        )
    );
}


void Foam::mechanicalLaw::makeRelFf()
{
    if (relFfPtr_.valid())
    {
        FatalErrorIn("void " + type() + "::makeRelFf()")
            << "pointer already set" << abort(FatalError);
    }

    relFfPtr_.set
    (
        new surfaceTensorField
        (
            IOobject
            (
                "relFf_",
                mesh().time().timeName(),
                mesh(),
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
            ),
            mesh(),
            dimensionedTensor("I", dimless, I)
        )
    );
}


void Foam::mechanicalLaw::makeSigmaHyd()
{
    if (sigmaHydPtr_.valid() || gradSigmaHydPtr_.valid())
    {
        FatalErrorIn("void " + type() + "::makeSigmaHyd()")
            << "pointer already set" << abort(FatalError);
    }

    sigmaHydPtr_.set
    (
        new volScalarField
        (
            IOobject
            (
                "sigmaHyd",
                mesh().time().timeName(),
                mesh(),
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            mesh(),
            dimensionedScalar("zero", dimPressure, 0.0),
            zeroGradientFvPatchScalarField::typeName
        )
    );

    gradSigmaHydPtr_.set
    (
        new volVectorField
        (
            IOobject
            (
                "grad(sigmaHyd)",
                mesh().time().timeName(),
                mesh(),
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
            ),
            mesh(),
            dimensionedVector("zero", dimPressure/dimLength, vector::zero)
        )
    );
}


void Foam::mechanicalLaw::makeSigmaEff()
{
    if (sigmaEffPtr_.valid())
    {
        FatalErrorIn("void " + type() + "::makeSigmaEff()")
            << "pointer already set" << abort(FatalError);
    }

    sigmaEffPtr_.set
    (
        new volSymmTensorField
        (
            IOobject
            (
                "sigmaEff_",
                mesh().time().timeName(),
                mesh(),
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            mesh(),
            dimensionedSymmTensor("zero", dimPressure, symmTensor::zero)
        )
    );
}


// * * * * * * * * * * * * * * Protected Members * * * * * * * * * * * * * * //

bool Foam::mechanicalLaw::planeStress() const
{
    if (mesh_.foundObject<IOdictionary>("mechanicalProperties"))
    {
        return Switch
        (
            mesh_.lookupObject<IOdictionary>
            (
                "mechanicalProperties"
            ).lookup("planeStress")
        );
    }
    else if
    (
        mesh_.objectRegistry::parent().foundObject<objectRegistry>("region0")
    )
    {
        return Switch
        (
            mesh_.objectRegistry::parent().subRegistry
            (
                "region0"
            ).lookupObject<IOdictionary>
            (
                "mechanicalProperties"
            ).lookup("planeStress")
        );
    }
    else
    {
        return Switch
        (
            mesh_.objectRegistry::parent().subRegistry
            (
                "solid"
            ).lookupObject<IOdictionary>
            (
                "mechanicalProperties"
            ).lookup("planeStress")
        );
    }
}


const Foam::volScalarField& Foam::mechanicalLaw::mu() const
{
    if (muPtr_.empty())
    {
        makeMu();
    }

    return muPtr_();
}


Foam::volScalarField& Foam::mechanicalLaw::mu()
{
    if (muPtr_.empty())
    {
        makeMu();
    }

    return muPtr_();
}


const Foam::surfaceScalarField& Foam::mechanicalLaw::muf() const
{
    if (mufPtr_.empty())
    {
        makeMuf();
    }

    return mufPtr_();
}


Foam::surfaceScalarField& Foam::mechanicalLaw::muf()
{
    if (mufPtr_.empty())
    {
        makeMuf();
    }

    return mufPtr_();
}


const Foam::volScalarField& Foam::mechanicalLaw::mu(const dimensionedScalar& mu)
{
    if (muPtr_.empty())
    {
        makeMu();
    }

    // Set constant field, overwriting any existing values
    muPtr_() = mu;

    // Return a reference to the field
    return muPtr_();
}


const Foam::surfaceScalarField& Foam::mechanicalLaw::muf
(
    const dimensionedScalar& mu
)
{
    if (mufPtr_.empty())
    {
        makeMuf();
    }

    // Set constant field, overwriting any existing values
    mufPtr_() = mu;

    // Return a reference to the field
    return mufPtr_();
}


const Foam::volScalarField& Foam::mechanicalLaw::K() const
{
    if (KPtr_.empty())
    {
        makeK();
    }

    return KPtr_();
}


Foam::volScalarField& Foam::mechanicalLaw::K()
{
    if (KPtr_.empty())
    {
        makeK();
    }

    return KPtr_();
}


const Foam::surfaceScalarField& Foam::mechanicalLaw::Kf() const
{
    if (KfPtr_.empty())
    {
        makeKf();
    }

    return KfPtr_();
}


Foam::surfaceScalarField& Foam::mechanicalLaw::Kf()
{
    if (KfPtr_.empty())
    {
        makeKf();
    }

    return KfPtr_();
}


const Foam::volScalarField& Foam::mechanicalLaw::K(const dimensionedScalar& K)
{
    if (KPtr_.empty())
    {
        makeK();
    }

    // Set constant field, overwriting any existing values
    KPtr_() = K;

    // Return a reference to the field
    return KPtr_();
}


const Foam::surfaceScalarField& Foam::mechanicalLaw::Kf
(
    const dimensionedScalar& K
)
{
    if (KfPtr_.empty())
    {
        makeKf();
    }

    // Set constant field, overwriting any existing values
    KfPtr_() = K;

    // Return a reference to the field
    return KfPtr_();
}


const Foam::volScalarField& Foam::mechanicalLaw::impK
(
    const dimensionedScalar& impK
)
{
    if (impKPtr_.empty())
    {
        // Construct the field
        impKPtr_.set
        (
            new volScalarField
            (
                IOobject
                (
                    "impK_" + name_,
                    mesh().time().timeName(),
                    mesh(),
                    IOobject::READ_IF_PRESENT,
                    IOobject::NO_WRITE
                ),
                mesh(),
                impK
            )
        );
    }

    // Set constant field, overwriting any existing values
    impKPtr_() = impK;

    // Return a reference to the field
    return impKPtr_();
}


const Foam::volSymmTensorField& Foam::mechanicalLaw::sigma0() const
{
    if (sigma0Ptr_.empty())
    {
        makeSigma0();
    }

    return sigma0Ptr_();
}


Foam::volSymmTensorField& Foam::mechanicalLaw::sigma0()
{
    if (sigma0Ptr_.empty())
    {
        makeSigma0();
    }

    return sigma0Ptr_();
}


const Foam::volSymmTensorField& Foam::mechanicalLaw::epsilon() const
{
    if (epsilonPtr_.empty())
    {
        makeEpsilon();
    }

    return epsilonPtr_();
}


Foam::volSymmTensorField& Foam::mechanicalLaw::epsilon()
{
    if (epsilonPtr_.empty())
    {
        makeEpsilon();
    }

    return epsilonPtr_();
}


const Foam::surfaceSymmTensorField& Foam::mechanicalLaw::epsilonf() const
{
    if (epsilonfPtr_.empty())
    {
        makeEpsilonf();
    }

    return epsilonfPtr_();
}


Foam::surfaceSymmTensorField& Foam::mechanicalLaw::epsilonf()
{
    if (epsilonfPtr_.empty())
    {
        makeEpsilonf();
    }

    return epsilonfPtr_();
}


Foam::volTensorField& Foam::mechanicalLaw::F()
{
    if (FPtr_.empty())
    {
        makeF();
    }

    return FPtr_();
}


Foam::surfaceTensorField& Foam::mechanicalLaw::Ff()
{
    if (FfPtr_.empty())
    {
        makeFf();
    }

    return FfPtr_();
}


Foam::volTensorField& Foam::mechanicalLaw::relF()
{
    if (relFPtr_.empty())
    {
        makeRelF();
    }

    return relFPtr_();
}


Foam::surfaceTensorField& Foam::mechanicalLaw::relFf()
{
    if (relFfPtr_.empty())
    {
        makeRelFf();
    }

    return relFfPtr_();
}


Foam::volScalarField& Foam::mechanicalLaw::sigmaHyd()
{
    if (sigmaHydPtr_.empty())
    {
        makeSigmaHyd();
    }

    return sigmaHydPtr_();
}


Foam::volVectorField& Foam::mechanicalLaw::gradSigmaHyd()
{
    if (gradSigmaHydPtr_.empty())
    {
        makeSigmaHyd();
    }

    return gradSigmaHydPtr_();
}


Foam::volSymmTensorField& Foam::mechanicalLaw::sigmaEff()
{
    if (sigmaEffPtr_.empty())
    {
        makeSigmaEff();
    }

    return sigmaEffPtr_();
}


void Foam::mechanicalLaw::updateEpsilon()
{
    if (incremental())
    {
        // Lookup gradient of displacement increment
        const volTensorField& gradDD =
            mesh().lookupObject<volTensorField>("grad(DD)");

        epsilon() = epsilon().oldTime() + symm(gradDD);
    }
    else
    {
        // Lookup gradient of displacement
        const volTensorField& gradD =
            mesh().lookupObject<volTensorField>("grad(D)");

        epsilon() = symm(gradD);
    }
}


void Foam::mechanicalLaw::updateEpsilonf()
{
    if (incremental())
    {
        // Lookup gradient of displacement increment
        const surfaceTensorField& gradDD =
            mesh().lookupObject<surfaceTensorField>("grad(DD)f");

        epsilonf() = epsilonf().oldTime() + symm(gradDD);
    }
    else
    {
        // Lookup gradient of displacement
        const surfaceTensorField& gradD =
            mesh().lookupObject<surfaceTensorField>("grad(D)f");

        epsilonf() = symm(gradD);
    }
}


bool Foam::mechanicalLaw::updateF
(
    volSymmTensorField& sigma,
    const dimensionedScalar& mu,
    const dimensionedScalar& K
)
{
    return updateF
    (
        sigma,
        mechanicalLaw::mu(mu),
        mechanicalLaw::K(K)
    );
}


bool Foam::mechanicalLaw::updateF
(
    volSymmTensorField& sigma,
    const volScalarField& mu,
    const volScalarField& K
)
{
    if (curTimeIndex_ != mesh().time().timeIndex())
    {
        curTimeIndex_ = mesh().time().timeIndex();
        warnAboutEnforceLinear_ = true;
    }

    // Check if the mathematical model is in total or updated Lagrangian form
    if (nonLinGeom() == nonLinearGeometry::UPDATED_LAGRANGIAN)
    {
        if (!incremental())
        {
            FatalErrorIn(type() + "::correct(volSymmTensorField& sigma)")
                << "Not implemented for non-incremental updated Lagrangian"
                << abort(FatalError);
        }

        // Lookup gradient of displacement increment
        const volTensorField& gradDD =
            mesh().lookupObject<volTensorField>("grad(DD)");

        // Calculate the relative deformation gradient
        relF() = I + gradDD.T();

        // Update the total deformation gradient
        F() = relF() & F().oldTime();

        if (enforceLinear())
        {
            if (warnAboutEnforceLinear_)
            {
                warnAboutEnforceLinear_ = false;

                WarningIn
                (
                    "void " + type() + "::"
                    "correct(volSymmTensorField& sigma)"
                )   << "Material linearity enforced for stability!" << endl;
            }

            // Calculate stress using Hooke's law
            sigma =
                sigma.oldTime()
              + 2.0*mu*symm(gradDD) + (K - (2.0/3.0)*mu)*tr(gradDD)*I;

            return true;
        }
    }
    else if (nonLinGeom() == nonLinearGeometry::TOTAL_LAGRANGIAN)
    {
        if (incremental())
        {
            // Lookup gradient of displacement increment
            const volTensorField& gradDD =
                mesh().lookupObject<volTensorField>("grad(DD)");

            // Update the total deformation gradient
            // Note: grad is wrt reference configuration
            F() = F().oldTime() + gradDD.T();

            // Update the relative deformation gradient: not needed
            relF() = F() & inv(F().oldTime());

            if (enforceLinear())
            {
                if (warnAboutEnforceLinear_)
                {
                    warnAboutEnforceLinear_ = false;

                    WarningIn
                    (
                        "void " + type() + "::correct(volSymmTensorField& sigma)"
                    )   << "Material linearity enforced for stability!" << endl;
                }

                // Calculate stress using Hooke's law
                sigma =
                    sigma.oldTime()
                  + 2.0*mu*dev(symm(gradDD)) + K*tr(gradDD)*I;

                return true;
            }
        }
        else
        {
            // Lookup gradient of displacement
            const volTensorField& gradD =
                mesh().lookupObject<volTensorField>("grad(D)");

            // Update the total deformation gradient
            F() = I + gradD.T();

            // Update the relative deformation gradient: not needed
            relF() = F() & inv(F().oldTime());

            if (enforceLinear())
            {
                if (warnAboutEnforceLinear_)
                {
                    warnAboutEnforceLinear_ = false;

                    WarningIn
                    (
                        "void " + type() + "::correct(volSymmTensorField& sigma)"
                    )   << "Material linearity enforced for stability!" << endl;
                }

                // Calculate stress using Hooke's law
                sigma = 2.0*mu*dev(symm(gradD)) + K*tr(gradD)*I;

                return true;
            }
        }
    }
    else
    {
        FatalErrorIn
        (
            "void " + type() + "::correct(volSymmTensorField& sigma)"
        )   << "Unknown nonLinGeom type: " << nonLinGeom() << abort(FatalError);
    }

    // linearised elasticity was not enforced
    return false;
}


bool Foam::mechanicalLaw::updateF
(
    surfaceSymmTensorField& sigma,
    const dimensionedScalar& mu,
    const dimensionedScalar& K
)
{
    return updateF
    (
        sigma,
        mechanicalLaw::muf(mu),
        mechanicalLaw::Kf(K)
    );
}


bool Foam::mechanicalLaw::updateF
(
    surfaceSymmTensorField& sigma,
    const surfaceScalarField& mu,
    const surfaceScalarField& K
)
{
    // Check if the mathematical model is in total or updated Lagrangian form
    if (nonLinGeom() == nonLinearGeometry::UPDATED_LAGRANGIAN)
    {
        if (!incremental())
        {
            FatalErrorIn(type() + "::correct(surfaceSymmTensorField& sigma)")
                << "Not implemented for non-incremental updated Lagrangian"
                << abort(FatalError);
        }

        // Lookup gradient of displacement increment
        const surfaceTensorField& gradDD =
            mesh().lookupObject<surfaceTensorField>("grad(DD)f");

        // Update the relative deformation gradient: not needed
        relFf() = I + gradDD.T();

        // Update the total deformation gradient
        Ff() = relFf() & Ff().oldTime();

        if (enforceLinear())
        {
            WarningIn
            (
                "void " + type() + "::correct(surfaceSymmTensorField& sigma)"
            )   << "Material linearity enforced for stability!" << endl;

            // Calculate stress using Hooke's law
            sigma =
                sigma.oldTime() + 2.0*mu*dev(symm(gradDD)) + K*tr(gradDD)*I;

            return true;
        }
    }
    else if (nonLinGeom() == nonLinearGeometry::TOTAL_LAGRANGIAN)
    {
        if (incremental())
        {
            // Lookup gradient of displacement increment
            const surfaceTensorField& gradDD =
                mesh().lookupObject<surfaceTensorField>("grad(DD)f");

            // Update the total deformation gradient
            // Note: grad is wrt reference configuration
            Ff() = Ff().oldTime() + gradDD.T();

            // Update the relative deformation gradient: not needed
            relFf() = Ff() & inv(Ff().oldTime());

            if (enforceLinear())
            {
                WarningIn
                (
                    "void " + type()
                  + "::correct(surfaceSymmTensorField& sigma)"
                )   << "Material linearity enforced for stability!" << endl;

                // Calculate stress using Hooke's law
                sigma =
                    sigma.oldTime()
                  + 2.0*mu*dev(symm(gradDD)) + K*tr(gradDD)*I;

                return true;
            }
        }
        else
        {
            // Lookup gradient of displacement
            const surfaceTensorField& gradD =
                mesh().lookupObject<surfaceTensorField>("grad(D)f");

            // Update the total deformation gradient
            Ff() = I + gradD.T();

            // Update the relative deformation gradient: not needed
            relFf() = Ff() & inv(Ff().oldTime());

            if (enforceLinear())
            {
                WarningIn
                (
                    "void " + type()
                  + "::correct(surfaceSymmTensorField& sigma)"
                )   << "Material linearity enforced for stability!" << endl;

                // Calculate stress using Hooke's law
                sigma = 2.0*mu*dev(symm(gradD)) + K*tr(gradD)*I;

                return true;
            }
        }
    }
    else
    {
        FatalErrorIn
        (
            "void " + type() + "::correct(surfaceSymmTensorField& sigma)"
        )   << "Unknown nonLinGeom type: " << nonLinGeom() << abort(FatalError);
    }

    // linearised elasticity was not enforced
    return false;
}


void Foam::mechanicalLaw::updateSigmaHyd
(
    const volScalarField& sigmaHydExplicit,
    const dimensionedScalar& impK
)
{
    updateSigmaHyd(sigmaHydExplicit, mechanicalLaw::impK(impK));
}


void Foam::mechanicalLaw::updateSigmaHyd
(
    const volScalarField& sigmaHydExplicit,
    const volScalarField& impK
)
{
    if (solvePressureEqn_)
    {
#ifdef OPENFOAMESIORFOUNDATION
        SolverPerformance<scalar>::debug = 0;
#endif

        // Store previous iteration to allow relaxation, if needed
        sigmaHyd().storePrevIter();

        // Lookup the momentum equation inverse diagonal field
        const volScalarField* ADPtr = NULL;
        bool allocatedMemory = false;
        if (baseMesh().foundObject<volScalarField>("DEqnA"))
        {
            if (mesh() == baseMesh())
            {
                ADPtr = &mesh().lookupObject<volScalarField>("DEqnA");
            }
            else
            {
#ifdef OPENFOAMESIORFOUNDATION
                FatalErrorIn("void Foam::mechanicalLaw::updateSigmaHyd(...)")
                    << "Multi-materials are not yet ported for this version of "
                    << "OpenFOAM" << abort(FatalError);
#else
                ADPtr =
                    new volScalarField
                    (
                        baseMesh().lookupObject<mechanicalModel>
                        (
                            "mechanicalProperties"
                        ).solSubMeshes().lookupBaseMeshVolField<scalar>
                        (
                            "DEqnA", mesh()
                        )
                    );
                allocatedMemory = true;
#endif
            }
        }
        else if (baseMesh().foundObject<volScalarField>("DDEqnA"))
        {
            notImplemented
            (
                "DDEqnA in solidModels should instead be called DEqnA"
            );
        }
        else
        {
            FatalErrorIn
            (
                "void " + type() + "updateSigmaHyd(...)\n"
            )   << "Cannot find the DEqnA or DDEqnA field: this should be "
                << "stored in the solidModel" << abort(FatalError);
        }
        const volScalarField& AD = *ADPtr;

        // Pressure diffusivity field
        const surfaceScalarField rDAf
        (
            "rDAf",
            pressureSmoothingScaleFactor_*fvc::interpolate
            (
                impK/AD, "interpolate(" + gradSigmaHyd().name() + ")"
            )
        );
        const dimensionedScalar one("one", dimless, 1.0);

        // Solve pressure laplacian
        // Note: the fvm and fvc laplacian terms cancel at convergence and the
        // laplacian - div(grad) term produce a smoothing/diffusion to quell
        // oscillations
        fvScalarMatrix sigmaHydEqn
        (
            fvm::Sp(one, sigmaHyd())
          - fvm::laplacian(rDAf, sigmaHyd(), "laplacian(rDA,sigmaHyd)")
         ==
            sigmaHydExplicit
          - fvc::div(rDAf*fvc::interpolate(gradSigmaHyd()) & mesh().Sf())
        );

        // Solve the pressure equation
        sigmaHydEqn.solve();

        // Relax the pressure field
        sigmaHyd().relax();

        if (allocatedMemory)
        {
            delete ADPtr;
        }
    }
    else
    {
        // Explicitly calculate hydrostatic stress
        // We use 1.0* to overwritting the field IOobject attributes e.g. its
        // name and writeOpt
        sigmaHyd() = 1.0*sigmaHydExplicit;
    }

    // Update the gradient
    gradSigmaHyd() = fvc::grad(sigmaHyd());
}


const Foam::Switch& Foam::mechanicalLaw::enforceLinear() const
{
    // Lookup the solideModel
    const solidModel& solMod = lookupSolidModel(mesh(), baseMeshRegionName_);

    return solMod.enforceLinear();
}


bool Foam::mechanicalLaw::incremental() const
{
    // Lookup the solideModel
    const solidModel& solMod = lookupSolidModel(mesh(), baseMeshRegionName_);

    return solMod.incremental();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::mechanicalLaw::mechanicalLaw
(
    const word& name,
    const fvMesh& mesh,
    const dictionary& dict,
    const nonLinearGeometry::nonLinearType& nonLinGeom
)
:
    name_(name),
    mesh_(mesh),
    dict_(dict),
    baseMeshRegionName_(),
    nonLinGeom_(nonLinGeom),
    muPtr_(),
    mufPtr_(),
    KPtr_(),
    KfPtr_(),
    epsilonPtr_(),
    epsilonfPtr_(),
    FPtr_(),
    FfPtr_(),
    relFPtr_(),
    relFfPtr_(),
    solvePressureEqn_
    (
        dict.lookupOrDefault<Switch>("solvePressureEqn", false)
    ),
    pressureSmoothingScaleFactor_
    (
        dict.lookupOrDefault<scalar>("pressureSmoothingScaleFactor", 100.0)
    ),
    sigmaHydPtr_(),
    gradSigmaHydPtr_(),
    sigmaEffPtr_(),
    curTimeIndex_(-1),
    warnAboutEnforceLinear_(true)
{
    // Set the base mesh region name
    // For an FSI case, the region will be called solid, else it will be called
    // region0.
    if (mesh.time().foundObject<fvMesh>("solid"))
    {
        baseMeshRegionName_ = "solid";
    }
    else if (mesh.time().foundObject<fvMesh>("region0"))
    {
        baseMeshRegionName_ = "region0";
    }
    else
    {
        FatalErrorIn
        (
            "Foam::mechanicalLaw::mechanicalLaw\n"
            "(\n"
            "    const word& name,\n"
            "    const fvMesh& mesh,\n"
            "    const dictionary& dict\n"
            ")"
        ) << "solid region name not found" << abort(FatalError);
    }

    if (solvePressureEqn_)
    {
        Info<< "    Laplacian equation will be solved for pressure" << nl
            << "    pressureSmoothingScaleFactor: "
            << pressureSmoothingScaleFactor_
            << endl;
    }
}


// * * * * * * * * * * * * * * * Member functions * * * * * * * * * * * * * * //


Foam::dimensionedScalar Foam::mechanicalLaw::rhoScalar() const
{
    return dimensionedScalar(dict_.lookup("rho"));
}


Foam::tmp<Foam::volScalarField> Foam::mechanicalLaw::rho() const
{
    tmp<volScalarField> tresult
    (
        new volScalarField
        (
            IOobject
            (
                "rho",
                mesh().time().timeName(),
                mesh(),
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
            ),
            mesh(),
            rhoScalar(),
            zeroGradientFvPatchScalarField::typeName
        )
    );

#ifdef OPENFOAMESIORFOUNDATION
    tresult.ref().correctBoundaryConditions();
#else
    tresult().correctBoundaryConditions();
#endif

    return tresult;
}


Foam::tmp<Foam::surfaceScalarField> Foam::mechanicalLaw::impKf() const
{
    return fvc::interpolate(impK());
}


Foam::tmp<Foam::Field<Foam::symmTensor4thOrder>>
Foam::mechanicalLaw::materialTangentField() const
{
    // Default to uniform field
    // This function can be overwritten in specific mechanical laws
    tmp<Field<symmTensor4thOrder>> tresult
    (
        new Field<symmTensor4thOrder>(mesh().nFaces(), materialTangent())
    );

    return tresult;
}


void Foam::mechanicalLaw::correct(surfaceSymmTensorField&)
{
    notImplemented
    (
        type() + "::correct(surfaceSymmTensorField&)\n"
        "The correct(surfaceSymmTensorField&) function is not implemented\n"
        " for the " + type() + " mechanical law"
    );
}


Foam::scalar Foam::mechanicalLaw::residual()
{
    // Default to zero; this can be overwritten by any derived mechanical law
    return 0.0;
}


Foam::scalar Foam::mechanicalLaw::newDeltaT()
{
    // Default to a large number
    return mesh_.time().endTime().value();
}


// ************************************************************************* //
