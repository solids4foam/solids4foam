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

#include "pdGuccioneElastic.H"
#include "addToRunTimeSelectionTable.H"
#include "fvc.H"
#include "zeroGradientFvPatchFields.H"
#include "extrapolatedFvPatchFields.H"
#include "eig3.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(pdGuccioneElastic, 0);
    addToRunTimeSelectionTable
    (
        mechanicalLaw, pdGuccioneElastic, nonLinGeomMechLaw
    );
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from dictionary
Foam::pdGuccioneElastic::pdGuccioneElastic
(
    const word& name,
    const fvMesh& mesh,
    const dictionary& dict,
    const nonLinearGeometry::nonLinearType& nonLinGeom
)
:
    mechanicalLaw(name, mesh, dict, nonLinGeom),
    bulkModulus_("K", dimPressure, GREAT),
    // bulkModulus_(dict.lookup("bulkModulus")),
    k_(dict.lookup("k")),
    cf_(readScalar(dict.lookup("cf"))),
    ct_(readScalar(dict.lookup("ct"))),
    cfs_(readScalar(dict.lookup("cfs"))),
    mu_(0.5*k_*(cf_ + cfs_ + ct_)/3),
    // mu_(0.75*(cf_ - 2.0*cfs_ + 2.0*cfs_)*k_*0.5),
    // check: is this the equivalent of linear shear modulus?
    muEff_
    (
        IOobject
        (
            "muEff",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        mu_,
        zeroGradientFvPatchScalarField::typeName
        // extrapolatedFvPatchScalarField::typeName
    ),
    f0_
    (
        IOobject
        (
            "f0",
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
    ),
    f0f_
    (
        IOobject
        (
            "f0f",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        fvc::interpolate(f0_)
    ),
    s0_
    (
        IOobject
        (
            "s0",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedVector("i", dimless, vector(1, 0, 0))
    ),
    s0f_
    (
        IOobject
        (
            "s0f",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedVector("i", dimless, vector(1, 0, 0))
    ),
    n0_
    (
        IOobject
        (
            "n0",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedVector("0", dimless, vector::zero)
    ),
    n0f_
    (
        IOobject
        (
            "n0f",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedVector("0", dimless, vector::zero)
    ),
    R_
    (
        IOobject
        (
            "R",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedTensor("0", dimless, tensor::zero)
    ),
    Rf_
    (
        IOobject
        (
            "Rf",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedTensor("0", dimless, tensor::zero)
    ),
    f0f0_("f0f0", sqr(f0_)),
    f0f0f_("f0f0f", sqr(f0f_)),
    S_
    (
        IOobject
        (
            "S2PK",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedSymmTensor("0", dimPressure, symmTensor::zero)
    ),
    Sf_
    (
        IOobject
        (
            "S2PKf",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedSymmTensor("0", dimPressure, symmTensor::zero)
    ),
    expQf_
    (
        IOobject
        (
            "expQf",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("1", dimless, 1)
    ),
    impKcoeff_
    (
        dict.lookupOrDefault<scalar>
        (
            "impKcoeff",
            1.0
        )
    )
{
    // Check f0 are unit vectors

#ifdef OPENFOAM_NOT_EXTEND
    if (min(mag(mag(f0_.primitiveField()))) < SMALL)
#else
    if (min(mag(mag(f0_.internalField()))) < SMALL)
#endif
    {
        FatalErrorIn
        (
            "pdGuccioneElastic::pdGuccioneElastic()"
        )
            << "At least one f0 vector has a length of zero!"
            << abort(FatalError);
    }

    // Normalise f0
    f0_ /= mag(f0_);
    f0f_ /= mag(f0f_);

    // Re-calculate f0f0
    f0f0_ = sqr(f0_);
    f0f0f_ = sqr(f0f_);

    // Store old F
    F().storeOldTime();
    Ff().storeOldTime();

    // Calculate sheet direction (s0) which is orthogonal to f0. There are
    // infinite vectors which are orthogonal to f0 so we will start with i
    // and remove any component in the f0 direction
    // Remove component in f0 direction
    s0_ = ((I - f0f0_) & s0_);
    s0f_ = ((I - f0f0f_) & s0f_);

    // Check for any vectors with zero magnitude; if found, then use the j
    // direction
    const volScalarField magS0
    (
        "magS0",
        mag(s0_) - dimensionedScalar("SMALL", dimless, SMALL)
    );
    const volScalarField posMagS0("posMagS0", pos(magS0));
    s0_ = posMagS0*s0_ + (1.0 - posMagS0)*((I - f0f0_) & vector(0, 1, 0));

    const surfaceScalarField magS0f
    (
        "magS0f",
        mag(s0f_) - dimensionedScalar("SMALL", dimless, SMALL)
    );
    const surfaceScalarField posMagS0f("posMagS0f", pos(magS0f));
    s0f_ = posMagS0f*s0f_ + (1.0 - posMagS0f)*((I - f0f0f_) & vector(0, 1, 0));

    // Make s0 unit vectors
    s0_ /= mag(s0_) + dimensionedScalar("SMALL", dimless, SMALL);
    s0f_ /= mag(s0f_) + dimensionedScalar("SMALL", dimless, SMALL);

    // Calculate n0 as orthogonal to f0 and s0
    n0_ = f0_ ^ s0_;
    n0_ /= mag(n0_) + dimensionedScalar("SMALL", dimless, SMALL);
    n0f_ = f0f_ ^ s0f_;
    n0f_ /= mag(n0f_) + dimensionedScalar("SMALL", dimless, SMALL);

    // Assign the components of R
    R_.replace(tensor::XX, f0_.component(vector::X));
    R_.replace(tensor::YX, f0_.component(vector::Y));
    R_.replace(tensor::ZX, f0_.component(vector::Z));
    R_.replace(tensor::XY, s0_.component(vector::X));
    R_.replace(tensor::YY, s0_.component(vector::Y));
    R_.replace(tensor::ZY, s0_.component(vector::Z));
    R_.replace(tensor::XZ, n0_.component(vector::X));
    R_.replace(tensor::YZ, n0_.component(vector::Y));
    R_.replace(tensor::ZZ, n0_.component(vector::Z));

    Rf_.replace(tensor::XX, f0f_.component(vector::X));
    Rf_.replace(tensor::YX, f0f_.component(vector::Y));
    Rf_.replace(tensor::ZX, f0f_.component(vector::Z));
    Rf_.replace(tensor::XY, s0f_.component(vector::X));
    Rf_.replace(tensor::YY, s0f_.component(vector::Y));
    Rf_.replace(tensor::ZY, s0f_.component(vector::Z));
    Rf_.replace(tensor::XZ, n0f_.component(vector::X));
    Rf_.replace(tensor::YZ, n0f_.component(vector::Y));
    Rf_.replace(tensor::ZZ, n0f_.component(vector::Z));

    if (dict.lookupOrDefault<Switch>("writeS0N0R", Switch(false)))
    {
        Info<< "Writing s0, n0 and R" << endl;
        s0_.write();
        n0_.write();
        R_.write();
    }

    // Update shear modulus
    calcInitShearModulus();

    // Update effective shear modulus
    calcEffectiveShearModulus();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::pdGuccioneElastic::~pdGuccioneElastic()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::pdGuccioneElastic::impK() const
{
    // Update effective shear modulus first
    // calcEffectiveShearModulus();

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
            impKcoeff_*muEff_
            // mesh()
            // impKcoeff_*mu_
            // 4*mu_
            // (cf_ - 2.0*cfs_ + 2.0*cfs_)*k_
            // (cf_ - 2.0*cfs_ + 2.0*cfs_)*k_ + bulkModulus_
        )
    );
}

Foam::tmp<Foam::volScalarField>
Foam::pdGuccioneElastic::bulkModulus() const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "bulkModulusLaw",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            bulkModulus_
        )
    );
}


void Foam::pdGuccioneElastic::correct(volSymmTensorField& sigma)
{
    // Update the deformation gradient field
    // Note: if true is returned, it means that linearised elasticity was
    // enforced by the solver via the enforceLinear switch
    if (updateF(sigma, mu_, bulkModulus_))
    {
        return;
    }

    // Take a reference to the deformation gradient to make the code easier to
    // read
    const volTensorField& F = this->F();
    const volTensorField FT("FT", F.T());

    // Calculate the Jacobian of the deformation gradient
    // const volScalarField J(det(F)); // J=1

    // Calculate the right Cauchy–Green deformation tensor
    const volSymmTensorField C("C", symm(FT & F));

    // Calculate the Green-Lagrange strain
    const volSymmTensorField E("E", 0.5*(C - I));

    const Switch useLocalCoordSys
    (
        dict().lookupOrDefault<Switch>
        (
            "calculateStressInLocalCoordinateSystem",
            Switch(false)
        )
    );

    if (useLocalCoordSys)
    {
        // Calculate the Green strain in the local coordinate system
        const volTensorField RT("RT", R_.T());
        const volSymmTensorField EStar("EStar", symm(RT & E & R_));

        // Extract the components of EStar
        // Note: EStar is symmetric
        const volScalarField E11("E11", EStar.component(symmTensor::XX));
        const volScalarField E12("E12", EStar.component(symmTensor::XY));
        const volScalarField E13("E13", EStar.component(symmTensor::XZ));
        const volScalarField E22("E22", EStar.component(symmTensor::YY));
        const volScalarField E23("E23", EStar.component(symmTensor::YZ));
        const volScalarField E33("E33", EStar.component(symmTensor::ZZ));

        // Calculate Q
        const volScalarField Q
        (
            "Q",
            cf_*sqr(E11)
          + ct_*(sqr(E22) + sqr(E33) + 2*sqr(E23))
          + cfs_*(2*sqr(E12) + 2*sqr(E13))
        );

        // Calculate the derivative of Q wrt to EStar
        volSymmTensorField dQdEStar
        (
            IOobject
            (
                "dQdEStar",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            dimensionedSymmTensor("0", dimless, symmTensor::zero)
        );

        dQdEStar.replace(symmTensor::XX, 2*cf_*E11);
        dQdEStar.replace(symmTensor::XY, 2*cfs_*E12);
        dQdEStar.replace(symmTensor::XZ, 2*cfs_*E13);
        dQdEStar.replace(symmTensor::YY, 2*ct_*E22);
        dQdEStar.replace(symmTensor::YZ, 2*ct_*E23);
        dQdEStar.replace(symmTensor::ZZ, 2*ct_*E33);

        // Calculate the local 2nd Piola-Kirchhoff stress (without the
        // hydrostatic term)
        S_ = dQdEStar*0.5*k_*exp(Q);

        // Rotate S from the local fibre coordinate system to the global
        // coordinate system
        S_ = symm(R_ & S_ & RT);
    }
    else
    {
        // Calculate E . E
        const volSymmTensorField sqrE("sqrE", symm(E & E));

        // Calculate the invariants of E
        const volScalarField I1("I1", tr(E));
        const volScalarField I2("I2", 0.5*(sqr(tr(E)) - tr(sqrE)));
        const volScalarField I4("I4", E && f0f0_);
        const volScalarField I5("I5", sqrE && f0f0_);

        // Calculate Q
        const volScalarField Q
        (
            ct_*sqr(I1)
          - 2.0*ct_*I2
         + (cf_ - 2.0*cfs_ + ct_)*sqr(I4)
         + 2.0*(cfs_ - ct_)*I5
        );

        // Calculate the derivative of Q wrt to E
        const volSymmTensorField dQdE
        (
            2.0*ct_*E
          + 2.0*(cf_ - 2.0*cfs_ + ct_)*I4*f0f0_
          + 2.0*(cfs_ - ct_)*symm((E & f0f0_) + (f0f0_ & E))
        );

        // Update the 2nd Piola-Kirchhoff stress (without the hydrostatic term)
        S_ = dQdE*0.5*k_*exp(Q);
    }

    // Convert the second Piola-Kirchhoff stress to the Cauchy stress and take
    // the deviatoric component
    // const volSymmTensorField s(symm(F & S_ & FT));
    const volSymmTensorField s("s", dev(symm(F & S_ & FT)));
    // const volSymmTensorField s(dev(J*symm(F & S_ & F.T())));

    // Calculate the hydrostatic stress
    // updateSigmaHyd
    // (
    //     0.5*bulkModulus_*(pow(J, 2.0) - 1.0)/J,
    //     (4.0/3.0)*mu_ + bulkModulus_
    // );

    // Lookup pressure field
    const volScalarField& p =
        mesh().lookupObject<volScalarField>("p");

    // Convert the second Piola-Kirchhoff deviatoric stress to the Cauchy stress
    // and add hydrostatic stress term

    sigma = s - p*I;
    // sigma = s + sigmaHyd()*I;
}


void Foam::pdGuccioneElastic::correct(surfaceSymmTensorField& sigma)
{
    // notImplemented("pdGuccioneElastic::correct(surfaceSymmTensorField& sigma)");

    // Update the deformation gradient field
    // Note: if true is returned, it means that linearised elasticity was
    // enforced by the solver via the enforceLinear switch
    if (updateF(sigma, mu_, bulkModulus_))
    {
        return;
    }

    expQf_.storePrevIter();

    // Take a reference to the deformation gradient to make the code easier to
    // read
    const surfaceTensorField& Ff = this->Ff();
    const surfaceTensorField FfT("FfT", Ff.T());

    // Calculate the Jacobian of the deformation gradient
    // const volScalarField J(det(F)); // J=1

    // Calculate the right Cauchy–Green deformation tensor
    const surfaceSymmTensorField C("C", symm(FfT & Ff));

    // Calculate the Green-Lagrange strain
    const surfaceSymmTensorField E("E", 0.5*(C - I));

    const Switch useLocalCoordSys
    (
        dict().lookupOrDefault<Switch>
        (
            "calculateStressInLocalCoordinateSystem",
            Switch(false)
        )
    );

    if (useLocalCoordSys)
    {
        // Calculate the Green strain in the local coordinate system
        const surfaceTensorField RfT("RfT", Rf_.T());
        const surfaceSymmTensorField EStar("EStar", symm(RfT & E & Rf_));

        // Extract the components of EStar
        // Note: EStar is symmetric
        const surfaceScalarField E11("E11", EStar.component(symmTensor::XX));
        const surfaceScalarField E12("E12", EStar.component(symmTensor::XY));
        const surfaceScalarField E13("E13", EStar.component(symmTensor::XZ));
        const surfaceScalarField E22("E22", EStar.component(symmTensor::YY));
        const surfaceScalarField E23("E23", EStar.component(symmTensor::YZ));
        const surfaceScalarField E33("E33", EStar.component(symmTensor::ZZ));

        // Calculate Q
        const surfaceScalarField Q
        (
            "Q",
            cf_*sqr(E11)
          + ct_*(sqr(E22) + sqr(E33) + 2*sqr(E23))
          + cfs_*(2*sqr(E12) + 2*sqr(E13))
        );

        // Calculate the derivative of Q wrt to EStar
        surfaceSymmTensorField dQdEStar
        (
            IOobject
            (
                "dQdEStar",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            dimensionedSymmTensor("0", dimless, symmTensor::zero)
        );

        dQdEStar.replace(symmTensor::XX, 2*cf_*E11);
        dQdEStar.replace(symmTensor::XY, 2*cfs_*E12);
        dQdEStar.replace(symmTensor::XZ, 2*cfs_*E13);
        dQdEStar.replace(symmTensor::YY, 2*ct_*E22);
        dQdEStar.replace(symmTensor::YZ, 2*ct_*E23);
        dQdEStar.replace(symmTensor::ZZ, 2*ct_*E33);

        expQf_ = exp(Q);
        expQf_.relax();

        // Calculate the local 2nd Piola-Kirchhoff stress (without the
        // hydrostatic term)
        Sf_ = dQdEStar*0.5*k_*expQf_;
        // Sf_ = dQdEStar*0.5*k_*exp(Q);

        // Rotate S from the local fibre coordinate system to the global
        // coordinate system
        Sf_ = symm(Rf_ & Sf_ & RfT);
    }
    else
    {
        // Calculate E.E
        const surfaceSymmTensorField sqrE("sqrE", symm(E & E));

        // Calculate the invariants of E
        const surfaceScalarField I1("I1", tr(E));
        const surfaceScalarField I2("I2", 0.5*(sqr(tr(E)) - tr(sqrE)));
        const surfaceScalarField I4("I4", E && f0f0f_);
        const surfaceScalarField I5("I5", sqrE && f0f0f_);

        // Calculate Q
        const surfaceScalarField Q
        (
            ct_*sqr(I1)
          - 2.0*ct_*I2
          + (cf_ - 2.0*cfs_ + ct_)*sqr(I4)
          + 2.0*(cfs_ - ct_)*I5
        );

        Info << "Q, min: " << gMin(Q.internalField())
             << ", avg: " << gAverage(Q.internalField())
             << ", max: " << gMax(Q.internalField()) << endl;

        // Calculate the derivative of Q wrt to E
        const surfaceSymmTensorField dQdE
        (
            2.0*ct_*E
          + 2.0*(cf_ - 2.0*cfs_ + ct_)*I4*f0f0f_
          + 2.0*(cfs_ - ct_)*symm((E & f0f0f_) + (f0f0f_ & E))
        );

        expQf_ = exp(Q);
        expQf_.relax();

        // Update the 2nd Piola-Kirchhoff stress (without the hydrostatic term)
        Sf_ = dQdE*0.5*k_*expQf_;
        // Sf_ = dQdE*0.5*k_*exp(Q);
    }

    // Convert the second Piola-Kirchhoff stress to the Cauchy stress and take
    // the deviatoric component
    // const surfaceSymmTensorField sf(symm(Ff & Sf_ & FfT));
    const surfaceSymmTensorField sf("sf", dev(symm(Ff & Sf_ & FfT)));
    // const volSymmTensorField s(dev(J*symm(F & S_ & F.T())));

    // Calculate the hydrostatic stress
    // updateSigmaHyd
    // (
    //     0.5*bulkModulus_*(pow(J, 2.0) - 1.0)/J,
    //     (4.0/3.0)*mu_ + bulkModulus_
    // );

    // Lookup pressure field
    const surfaceScalarField& pf =
        mesh().lookupObject<surfaceScalarField>("pf");

    // Convert the second Piola-Kirchhoff deviatoric stress to the Cauchy stress
    // and add hydrostatic stress term
    sigma = sf - pf*I;
    // sigma = s + sigmaHyd()*I;
}


void Foam::pdGuccioneElastic::setRestart()
{
    F().writeOpt() = IOobject::AUTO_WRITE;
    Ff().writeOpt() = IOobject::AUTO_WRITE;
}

//- Calculate Cauchy stress as a function of deformation gradient
void Foam::pdGuccioneElastic::calcDevCauchy
(
    const tensor& F,
    const symmTensor& f0f0,
    const tensor& R,
    symmTensor& devSigma
) const
{
    // Calculate the right Cauchy–Green deformation tensor
    const tensor FT(F.T());
    const symmTensor C(symm(FT & F));

    // Calculate the Green-Lagrange strain
    const symmTensor E(0.5*(C - I));

    const Switch useLocalCoordSys
    (
        dict().lookupOrDefault<Switch>
        (
            "calculateStressInLocalCoordinateSystem",
            Switch(false)
        )
    );

    symmTensor S = symmTensor::zero;

    if (useLocalCoordSys)
    {
        // Calculate the Green strain in the local coordinate system
        const tensor RT(R.T());
        const symmTensor EStar(symm(RT & E & R));

        // Extract the components of EStar
        // Note: EStar is symmetric
        const scalar E11(EStar.xx());
        const scalar E12(EStar.xy());
        const scalar E13(EStar.xz());
        const scalar E22(EStar.yy());
        const scalar E23(EStar.yz());
        const scalar E33(EStar.zz());

        // Calculate Q
        const scalar Q
        (
            cf_*sqr(E11)
          + ct_*(sqr(E22) + sqr(E33) + 2*sqr(E23))
          + cfs_*(2*sqr(E12) + 2*sqr(E13))
        );

        // Calculate the derivative of Q wrt to EStar
        symmTensor dQdEStar = symmTensor::zero;

        dQdEStar.xx() = 2*cf_*E11;
        dQdEStar.xy() = 2*cfs_*E12;
        dQdEStar.xz() = 2*cfs_*E13;
        dQdEStar.yy() = 2*ct_*E22;
        dQdEStar.yz() = 2*ct_*E23;
        dQdEStar.zz() = 2*ct_*E33;

        // Calculate the local 2nd Piola-Kirchhoff stress
        // (without the hydrostatic term)
        S = dQdEStar*0.5*k_.value()*exp(Q);

        // Rotate S from the local fibre coordinate system
        // to the global coordinate system
        S = dev(symm(R & S & RT));
    }
    else
    {
        // Calculate E . E
        const symmTensor sqrE(symm(E & E));

        // Calculate the invariants of E
        const scalar I1(tr(E));
        const scalar I2(0.5*(sqr(tr(E)) - tr(sqrE)));
        const scalar I4(E && f0f0);
        const scalar I5(sqrE && f0f0);

        // Calculate Q
        const scalar Q
        (
            ct_*sqr(I1)
          - 2.0*ct_*I2
          + (cf_ - 2.0*cfs_ + ct_)*sqr(I4)
          + 2.0*(cfs_ - ct_)*I5
        );

        // Calculate the derivative of Q wrt to E
        const symmTensor dQdE
        (
            2.0*ct_*E
          + 2.0*(cf_ - 2.0*cfs_ + ct_)*I4*f0f0
          + 2.0*(cfs_ - ct_)*symm((E & f0f0) + (f0f0 & E))
        );

        // Calculate the 2nd Piola-Kirchhoff stress
        // (without the hydrostatic term)
        S = dQdE*0.5*k_.value()*::exp(Q);
    }

    // Convert the second Piola-Kirchhoff stress to the Cauchy stress and take
    // the deviatoric component
    devSigma = dev(symm(F & S & FT));
    // devSigma = (symm(F & S & F.T()));
    // devSigma = (dev(symm(F & S & F.T())));
}

void Foam::pdGuccioneElastic::calcInitShearModulus()
{
    // max(mu_12, mu_13, mu_23))

    scalarField mu(3, 0);
    {
        scalar gamma = 0.001;
        tensor F = tensor::I;
        F.xy() = gamma;

        symmTensor devSigma = symmTensor::zero;
        calcDevCauchy(F, f0f0_[0], R_[0], devSigma);

        vector dTraction = (vector(0,1,0) & devSigma);

        mu[0] = mag(dTraction.x()/gamma);
    }

    {
        scalar gamma = 0.001;
        tensor F = tensor::I;
        F.xz() = gamma;

        symmTensor devSigma = symmTensor::zero;
        calcDevCauchy(F, f0f0_[0], R_[0], devSigma);

        vector dTraction = (vector(0,0,1) & devSigma);

        mu[1] = mag(dTraction.x()/gamma);
    }

    {
        scalar gamma = 0.001;
        tensor F = tensor::I;
        F.yz() = gamma;

        symmTensor devSigma = symmTensor::zero;
        calcDevCauchy(F, f0f0_[0], R_[0], devSigma);

        vector dTraction = (vector(0,0,1) & devSigma);

        mu[2] = mag(dTraction.y()/gamma);
    }

    Info << "m12 = " << mu[0] << endl;
    Info << "m13 = " << mu[1] << endl;
    Info << "m23 = " << mu[2] << endl;
    Info << "Current mu = " << mu_.value() << endl;

    mu_.value() = max(mu);
}


void Foam::pdGuccioneElastic::calcEffectiveShearModulus()
{
    // This function will update shear modulus muEff_ at the end of time step,
    // using current state of deformation

    const volTensorField& F = this->F();
    const tensorField& FI = F.internalField();
    const tensorField FIT(FI.T()());

    const symmTensorField& f0f0I = f0f0_.internalField();
    const tensorField& RI = R_.internalField();

#ifdef OPENFOAM_NOT_EXTEND
    scalarField& muEffI = muEff_.primitiveFieldRef();
#else
    scalarField& muEffI = muEff_.internalField();
#endif

    forAll(muEffI, cellI)
    {
        scalarField mu(3, 0);

        symmTensor devSigma = symmTensor::zero;
        calcDevCauchy(FI[cellI], f0f0I[cellI], RI[cellI], devSigma);

        tensor pVectors = tensor::zero;
        vector pValues = vector::zero;
        eig3().eigen_decomposition(devSigma, pVectors, pValues);

        tensor pVectorsT = pVectors.T();
        devSigma = symm(pVectors & devSigma & pVectorsT);

        // Info << "cellI: " << cellI << endl;
        // Info << "devSigma: " << devSigma << endl;
        // Info << "pVectors: " << pVectors << endl;
        // Info << "pValues: " << pValues << endl;
        // Info << "transformedDevSigma: " << transformedDevSigma << endl;

        // mu12
        {
            // Disp gradient perturabation in principal coordinate system
            tensor pertGradDD = symmTensor::zero;
            pertGradDD.xy() = 0.001;
            pertGradDD.yx() = 0.001;

            // Transform from local to globa CS

            pertGradDD = (pVectorsT & pertGradDD & pVectors);
            tensor pertGradDDT(pertGradDD.T());

            tensor DF = pertGradDDT;
            tensor newF = FI[cellI] + DF;

            symmTensor newDevSigma = symmTensor::zero;
            calcDevCauchy(newF, f0f0I[cellI], RI[cellI], newDevSigma);

            // Trsnform from global to local CS
            newDevSigma = symm(pVectors & newDevSigma & pVectorsT);

            symmTensor DDevSigma = newDevSigma - devSigma;

            symmTensor DEffStrain =
                symm
                (
                    (FI[cellI] & pertGradDD) +
                    (pertGradDDT & FIT[cellI])
                );

            // Transform from global to local CS
            DEffStrain = symm(pVectors & DEffStrain & pVectorsT);

            scalar mu12 = DDevSigma.xy()/(DEffStrain.xy()+SMALL);

            mu[0] = mu12;
        }

        // mu13
        {
            // In local CS
            tensor pertGradDD = symmTensor::zero;
            pertGradDD.xz() = 0.001;
            pertGradDD.zx() = 0.001;

            // Transform from local to globa CS
            pertGradDD = (pVectors.T() & pertGradDD & pVectors);
            tensor pertGradDDT(pertGradDD.T());

            tensor DF = pertGradDDT;
            tensor newFI = FI[cellI] + DF;

            symmTensor newDevSigma = symmTensor::zero;
            calcDevCauchy(newFI, f0f0I[cellI], RI[cellI], newDevSigma);

            // Trsnform from global to local CS
            newDevSigma = symm(pVectors & newDevSigma & pVectorsT);

            symmTensor DDevSigma = newDevSigma - devSigma;

            symmTensor DEffStrain =
                symm
                (
                    (FI[cellI] & pertGradDD) +
                    (pertGradDDT & FIT[cellI])
                );

            // Transform from global to local CS
            DEffStrain = symm(pVectors & DEffStrain & pVectorsT);

            scalar mu13 = DDevSigma.xz()/(DEffStrain.xz()+SMALL);

            mu[1] = mu13;
        }

        // mu23
        {
            tensor pertGradDD = symmTensor::zero;
            pertGradDD.yz() = 0.001;
            pertGradDD.zy() = 0.001;

            // Transform from local to globa CS
            pertGradDD = (pVectorsT & pertGradDD & pVectors);
            tensor pertGradDDT(pertGradDD.T());

            tensor DF = pertGradDDT;
            tensor newFI = FI[cellI] + DF;

            symmTensor newDevSigma = symmTensor::zero;
            calcDevCauchy(newFI, f0f0I[cellI], RI[cellI], newDevSigma);

            // Trsnform from global to local CS
            newDevSigma = symm(pVectors & newDevSigma & pVectorsT);

            symmTensor DDevSigma = newDevSigma - devSigma;

            symmTensor DEffStrain =
                symm
                (
                    (FI[cellI] & pertGradDD) +
                    (pertGradDDT & FIT[cellI])
                );

            // Transform from global to local CS
            DEffStrain = symm(pVectors & DEffStrain & pVectorsT);

            scalar mu23 = DDevSigma.yz()/(DEffStrain.yz()+SMALL);

            mu[2] = mu23;
        }

        muEffI[cellI] = max(mu);
    }

    muEff_.correctBoundaryConditions();
}

void Foam::pdGuccioneElastic::updateTotalFields()
{
    calcEffectiveShearModulus();
}

// ************************************************************************* //
