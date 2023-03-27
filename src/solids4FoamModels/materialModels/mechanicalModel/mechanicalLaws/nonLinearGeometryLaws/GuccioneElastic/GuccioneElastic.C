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

\*---------------------------------------------------------------------------*/

#include "GuccioneElastic.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(GuccioneElastic, 0);
    addToRunTimeSelectionTable
    (
        mechanicalLaw, GuccioneElastic, nonLinGeomMechLaw
    );
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from dictionary
Foam::GuccioneElastic::GuccioneElastic
(
    const word& name,
    const fvMesh& mesh,
    const dictionary& dict,
    const nonLinearGeometry::nonLinearType& nonLinGeom
)
:
    mechanicalLaw(name, mesh, dict, nonLinGeom),
    bulkModulus_(dict.lookup("bulkModulus")),
    k_(dict.lookup("k")),
    cf_(readScalar(dict.lookup("cf"))),
    ct_(readScalar(dict.lookup("ct"))),
    cfs_(readScalar(dict.lookup("cfs"))),
    mu_(0.75*(cf_ - 2.0*cfs_ + 2.0*cfs_)*k_), // check: is this the equivalent of linear shear modulus?
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
    f0f0_("f0f0", sqr(f0_)),
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
    )
{
    // Check f0 are unit vectors

    if (min(mag(mag(f0_.primitiveField()))) < SMALL)
    {
        FatalErrorIn("GuccioneElastic::GuccioneElastic()")
            << "At least one f0 vector has a length of zero!"
            << abort(FatalError);
    }

    // Normalise f0
    f0_ /= mag(f0_);

    // Re-calculate f0f0
    f0f0_ = sqr(f0_);

    // Store old F
    F().storeOldTime();
    //Ff().storeOldTime();

    // Calculate sheet direction (s0) which is orthogonal to f0. There are
    // infinite vectors which are orthogonal to f0 so we will start with i
    // and remove any component in the f0 direction
    // Remove component in f0 direction
    s0_ = ((I - f0f0_) & s0_);

    // Check for any vectors with zero magnitude; if found, then use the j
    // direction
    const volScalarField magS0(mag(s0_));
    const volScalarField posMagS0(pos(magS0));
    s0_ = posMagS0*s0_ + (1.0 - posMagS0)*((I - f0f0_) & vector(0, 1, 0));

    // Make s0 unit vectors
    s0_ /= mag(s0_);

    // Calculate n0 as orthogonal to f0 and s0
    n0_ = f0_ ^ s0_;
    n0_ /= mag(n0_);

    // Assign the components of R
    R_.replace(tensor::XX, f0_.component(vector::X));
    R_.replace(tensor::YX, f0_.component(vector::Y));
    R_.replace(tensor::ZX, f0_.component(vector::Z));
    R_.replace(tensor::YY, s0_.component(vector::X));
    R_.replace(tensor::YY, s0_.component(vector::Y));
    R_.replace(tensor::ZY, s0_.component(vector::Z));
    R_.replace(tensor::YZ, n0_.component(vector::X));
    R_.replace(tensor::YZ, n0_.component(vector::Y));
    R_.replace(tensor::ZZ, n0_.component(vector::Z));

    if (dict.lookupOrDefault<Switch>("writeS0N0R", Switch(false)))
    {
        Info<< "Writing s0, n0 and R" << endl;
        s0_.write();
        n0_.write();
        R_.write();
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::GuccioneElastic::~GuccioneElastic()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::GuccioneElastic::impK() const
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
            (cf_ - 2.0*cfs_ + 2.0*cfs_)*k_ + bulkModulus_
        )
    );
}


void Foam::GuccioneElastic::correct(volSymmTensorField& sigma)
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

    // Calculate the Jacobian of the deformation gradient
    const volScalarField J(det(F));

    // Calculate the right Cauchy–Green deformation tensor
    const volSymmTensorField C(symm(F.T() & F));
    
    // Calculate the Green-Lagrange strain
    const volSymmTensorField E(0.5*(C - I));

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
        const volSymmTensorField EStar("EStar", symm(R_.T() & E & R_));

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
        S_ = symm(R_ & S_ & R_.T());
    }
    else
    {
        // Calculate E . E
        const volSymmTensorField sqrE(symm(E & E));

        // Calculate the invariants of E
        const volScalarField I1(tr(E));
        const volScalarField I2(0.5*(sqr(tr(E)) - tr(sqrE)));
        const volScalarField I4(E && f0f0_);
        const volScalarField I5(sqrE && f0f0_);

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
    const volSymmTensorField s(dev(J*symm(F & S_ & F.T())));

    // Calculate the hydrostatic stress
    updateSigmaHyd
    (
        0.5*bulkModulus_*(pow(J, 2.0) - 1.0)/J,
        (4.0/3.0)*mu_ + bulkModulus_
    );

    // Convert the second Piola-Kirchhoff deviatoric stress to the Cauchy stress
    // and add hydrostatic stress term
    sigma = s + sigmaHyd()*I;
}


void Foam::GuccioneElastic::correct(surfaceSymmTensorField& sigma)
{
    notImplemented("GuccioneElastic::correct(surfaceSymmTensorField& sigma)");
}


void Foam::GuccioneElastic::setRestart()
{
    F().writeOpt() = IOobject::AUTO_WRITE;
    //Ff().writeOpt() = IOobject::AUTO_WRITE;
}

// ************************************************************************* //
