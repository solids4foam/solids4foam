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

#include "HolzapfelGasserOgdenElastic.H"
#include "addToRunTimeSelectionTable.H"
#include "calculatedFvPatchFields.H"
#include "eig3.H"
#include "zeroGradientFvPatchFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(HolzapfelGasserOgdenElastic, 0);
    addToRunTimeSelectionTable
    (
        mechanicalLaw, HolzapfelGasserOgdenElastic, nonLinGeomMechLaw
    );
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from dictionary
Foam::HolzapfelGasserOgdenElastic::HolzapfelGasserOgdenElastic
(
    const word& name,
    const fvMesh& mesh,
    const dictionary& dict,
    const nonLinearGeometry::nonLinearType& nonLinGeom
)
:
    mechanicalLaw(name, mesh, dict, nonLinGeom),
    rho_(dict.lookup("rho")),
    mu_("mu", dimPressure, 0.0),
    muEff_
    (
        IOobject
        (
            "muEff",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        mu_,
        zeroGradientFvPatchScalarField::typeName
    ),
    K_("K", dimPressure, 0.0),
    fibreAngle_(dict.lookup("fibreAngle")),
    k1_(dict.lookup("k1")),
    k2_(dict.lookup("k2")),
    Ec_
    (
        IOobject
        (
            "Ec",
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),
    Ea_
    (
        IOobject
        (
            "Ea",
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),
    Er_
    (
        IOobject
        (
            "Er",
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),
    Ecf_
    (
        IOobject
        (
            "Ecf",
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),
    Eaf_
    (
        IOobject
        (
            "Eaf",
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),
    Erf_
    (
        IOobject
        (
            "Erf",
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
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
    // Read mechanical properties
    if
    (
        dict.found("E") && dict.found("nu")
     && !dict.found("mu")
    )
    {
        const dimensionedScalar E = dimensionedScalar(dict.lookup("E"));
        const dimensionedScalar nu = dimensionedScalar(dict.lookup("nu"));

        mu_ = (E/(2.0*(1.0 + nu)));

        // Set the bulk modulus
        if (nu.value() < 0.5-SMALL)
        {
            FatalErrorIn(type())
                << "This is incompressible solid model (nu=0.5)"
                << abort(FatalError);
            // K_ = (nu*E/((1.0 + nu)*(1.0 - 2.0*nu))) + (2.0/3.0)*mu_;
        }
        else
        {
            K_.value() = GREAT;
        }
    }
    else if
    (
        dict.found("mu")
     && !dict.found("E") && !dict.found("nu")
    )
    {
        mu_ = dimensionedScalar(dict.lookup("mu"));
        K_.value() = GREAT;
    }
    else
    {
        FatalErrorIn(type())
            << "Either E and nu or mu should be specified"
            << abort(FatalError);
    }

    // Convert degrees to radians
    fibreAngle_.value() = fibreAngle_.value()*M_PI/180;

    // Update shear modulus
    calcInitialShearModulus();

    // Update effective shear modulus
    calcEffectiveShearModulus();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::HolzapfelGasserOgdenElastic::~HolzapfelGasserOgdenElastic()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::HolzapfelGasserOgdenElastic::rho() const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "rhoLaw",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            rho_,
            calculatedFvPatchScalarField::typeName
        )
    );
}


Foam::tmp<Foam::volScalarField> Foam::HolzapfelGasserOgdenElastic::impK() const
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
            // mesh(),
            // impKcoeff_*mu_
            impKcoeff_*muEff_
        )
    );
}


Foam::tmp<Foam::volScalarField> Foam::HolzapfelGasserOgdenElastic::bulkModulus() const
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
            K_
        )
    );
}


void Foam::HolzapfelGasserOgdenElastic::correct(volSymmTensorField& sigma)
{
    // Update the deformation gradient field
    // Note: if true is returned, it means that linearised elasticity was
    // enforced by the solver via the enforceLinear switch
    dimensionedScalar K_("K", mu_.dimensions(), 0);
    if (updateF(sigma, mu_, K_))
    {
        return;
    }

    // NOTE [IMPORTANT]:
    // Do NOT write F.T() & F directly: see the comment in
    // StVenantKirchhoffElastic.C
    const volTensorField& F = this->F();
    const volTensorField FT("FT", F.T());

    // Calculate the Jacobian of the deformation gradient
    // const volScalarField J = det(F()); / J=1

    // Calculate the volume preserving left Cauchy Green strain
    const volSymmTensorField b("b", symm(F & FT));

    // Calculate the deviatoric stress
    const volSymmTensorField s("s", mu_*(b - I));

    // Lookup pressure field
    const volScalarField& p =
        mesh().lookupObject<volScalarField>("p");

    // Calculate the Cauchy stress
    sigma = -p*I + s;


    // Anisotropic (HolzapfelGasserOgdenElastic) part

    //- Right Cauchy-Green tenso
    const volSymmTensorField C("C", symm(FT & F));

    //- Fibre directions
    const volVectorField N4
    (
        "N4",
        cos(fibreAngle_)*Ec_ +
        sin(fibreAngle_)*Ea_ +
        0.0*Er_
    );
    const volVectorField n4("n4", F & N4);

    const volVectorField N6
    (
        "N6",
       -cos(fibreAngle_)*Ec_ +
        sin(fibreAngle_)*Ea_ +
        0.0*Er_
    );
    const volVectorField n6("n6", F & N6);

    //- Direction matrices
    const volSymmTensorField N4N4("N4N4", symm(N4*N4));
    const volSymmTensorField N6N6("N6N6", symm(N6*N6));
    const volSymmTensorField n4n4("n4n4", symm(n4*n4));
    const volSymmTensorField n6n6("n6n6", symm(n6*n6));

    //- Fibre invariants
    const volScalarField I4("I4", C && N4N4);
    const volScalarField I6("I6", C && N6N6);

    //- Anisotropic part of Cauchy stress
    const volSymmTensorField sigmaAniso4
    (
        "sigmaAniso4",
        2*k1_*(I4-1.0)*exp(k2_*pow(I4-1.0, 2))*n4n4
    );
    const volSymmTensorField sigmaAniso6
    (
        "sigmaAniso6",
        2*k1_*(I6-1.0)*exp(k2_*pow(I6-1.0, 2))*n6n6
    );

    //- Add anisotropic part of Cauchy stress
    sigma += sigmaAniso4;
    sigma += sigmaAniso6;
}


void Foam::HolzapfelGasserOgdenElastic::correctF
(
    surfaceSymmTensorField& sigma,
    const surfaceTensorField& Ff
)
{
    // NOTE [IMPORTANT]:
    // Do NOT write F.T() & F directly: see the comment in
    // StVenantKirchhoffElastic.C
    const surfaceTensorField FfT("FfT", Ff.T());

    // Calculate left Cauchy Green strain tensor with volumetric term removed
    const surfaceSymmTensorField b("b", symm(Ff & FfT));

    // Calculate deviatoric stress
    const surfaceSymmTensorField s("s", mu_*(b-I));

    // Lookup pressure field
    const surfaceScalarField& pf =
        mesh().lookupObject<surfaceScalarField>("pf");

    // Calculate the Cauchy stress
    sigma = -pf*I + s;

    // Anisotropic (HolzapfelGasserOgdenElastic) part

    //- Right Cauchy-Green tensor
    const surfaceSymmTensorField Cf("Cf", symm(FfT & Ff));

    //- Fibre directions
    const surfaceVectorField N4
    (
        "N4",
        cos(fibreAngle_)*Ecf_ +
        sin(fibreAngle_)*Eaf_ +
        0.0*Erf_
    );
    const surfaceVectorField n4("n4", Ff & N4);

    const surfaceVectorField N6
    (
        "N6",
       -cos(fibreAngle_)*Ecf_ +
        sin(fibreAngle_)*Eaf_ +
        0.0*Erf_
    );
    const surfaceVectorField n6("n6", Ff & N6);

    //- Direction matrices
    const surfaceSymmTensorField N4N4("N4N4", symm(N4*N4));
    const surfaceSymmTensorField N6N6("N6N6", symm(N6*N6));
    const surfaceSymmTensorField n4n4("n4n4", symm(n4*n4));
    const surfaceSymmTensorField n6n6("n6n6", symm(n6*n6));

    //- Fibre invariants
    const surfaceScalarField I4("I4", Cf && N4N4);
    const surfaceScalarField I6("I6", Cf && N6N6);

    //- Anisotropic part of Cauchy stress
    const surfaceSymmTensorField sigmaAniso4
    (
        "sigmaAniso4",
        2*k1_*(I4-1.0)*exp(k2_*pow(I4-1.0, 2))*n4n4
    );
    const surfaceSymmTensorField sigmaAniso6
    (
        "sigmaAniso6",
        2*k1_*(I6-1.0)*exp(k2_*pow(I6-1.0, 2))*n6n6
    );

    //- Add anisotropic part of Cauchy stress
    sigma += sigmaAniso4;
    sigma += sigmaAniso6;
}


void Foam::HolzapfelGasserOgdenElastic::correct(surfaceSymmTensorField& sigma)
{
    // Update the deformation gradient field
    // Note: if true is returned, it means that linearised elasticity was
    // enforced by the solver via the enforceLinear switch
    dimensionedScalar K_("K", mu_.dimensions(), 0);
    if (updateF(sigma, mu_, K_))
    {
        return;
    }

    correctF(sigma, this->Ff());
}

void Foam::HolzapfelGasserOgdenElastic::setRestart()
{
    F().writeOpt() = IOobject::AUTO_WRITE;
    Ff().writeOpt() = IOobject::AUTO_WRITE;
}

//- Calculate Cauchy stress as a function of deformation gradient
void Foam::HolzapfelGasserOgdenElastic::calcDevCauchy
(
    const tensor& F,
    const vector& Ec,
    const vector& Ea,
    const vector& Er,
    symmTensor& devSigma
) const
{
    // Calculate the volume preserving left Cauchy Green strain
    const tensor FT = F.T();
    const symmTensor b = symm(F & FT);

    // Calculate the deviatoric stress
    const symmTensor s = mu_.value()*(b - I);

    // Calculate the Cauchy stress
    devSigma = s;


    // Anisotropic (HolzapfelGasserOgdenElastic) part

    //- Right Cauchy-Green tenso
    const symmTensor C = symm(FT & F);

    //- Fibre directions
    const vector N4 =
        cos(fibreAngle_.value())*Ec +
        sin(fibreAngle_.value())*Ea +
        0.0*Er;
    const vector n4 = (F & N4);

    const vector N6 =
       -cos(fibreAngle_.value())*Ec +
        sin(fibreAngle_.value())*Ea +
        0.0*Er;
    const vector n6 = (F & N6);

    //- Direction matrices
    const symmTensor N4N4 = symm(N4*N4);
    const symmTensor N6N6 = symm(N6*N6);
    const symmTensor n4n4 = symm(n4*n4);
    const symmTensor n6n6 = symm(n6*n6);

    //- Fibre invariants
    const scalar I4 = (C && N4N4);
    const scalar I6 = (C && N6N6);

    //- Anisotropic part of Cauchy stress
    const symmTensor sigmaAniso4 =
        2*k1_.value()*(I4-1.0)*exp(k2_.value()*pow(I4-1.0, 2))*n4n4;
    const symmTensor sigmaAniso6 =
        2*k1_.value()*(I6-1.0)*exp(k2_.value()*pow(I6-1.0, 2))*n6n6;

    //- Add anisotropic part of Cauchy stress
    devSigma += sigmaAniso4;
    devSigma += sigmaAniso6;
}

void Foam::HolzapfelGasserOgdenElastic::calcInitialShearModulus()
{
    // max(mu_12, mu_13, mu_23))

    scalarField mu(3, 0);
    {
        scalar gamma = 0.001;
        tensor F = tensor::I;
        F.xy() = gamma;

        symmTensor devSigma = symmTensor::zero;
        calcDevCauchy(F, Ec_[0], Ea_[0], Er_[0], devSigma);

        vector dTraction = (vector(0,1,0) & devSigma);

        mu[0] = mag(dTraction.x()/gamma);
    }

    {
        scalar gamma = 0.001;
        tensor F = tensor::I;
        F.xz() = gamma;

        symmTensor devSigma = symmTensor::zero;
        calcDevCauchy(F, Ec_[0], Ea_[0], Er_[0], devSigma);

        vector dTraction = (vector(0,0,1) & devSigma);

        mu[1] = mag(dTraction.x()/gamma);
    }

    {
        scalar gamma = 0.001;
        tensor F = tensor::I;
        F.yz() = gamma;

        symmTensor devSigma = symmTensor::zero;
        calcDevCauchy(F, Ec_[0], Ea_[0], Er_[0], devSigma);

        vector dTraction = (vector(0,0,1) & devSigma);

        mu[2] = mag(dTraction.y()/gamma);
    }

    // Info << "m12 = " << mu[0] << endl;
    // Info << "m13 = " << mu[1] << endl;
    // Info << "m23 = " << mu[2] << endl;

    Info << "Specified mu = " << mu_.value() << endl;
    Info << "Calculated initial mu[]: " << mu << endl;

    // mu_.value() = max(mu);
}


void Foam::HolzapfelGasserOgdenElastic::calcEffectiveShearModulus()
{
    // This function will update shear modulus muEff_ at the end of time step,
    // using current state of deformation

    const volTensorField& F = this->F();
    const tensorField& FI = F.internalField();
    const tensorField FIT(FI.T()());

    const vectorField& EcI = Ec_.internalField();
    const vectorField& EaI = Ea_.internalField();
    const vectorField& ErI = Er_.internalField();

#ifdef OPENFOAM_NOT_EXTEND
    scalarField& muEffI = muEff_.primitiveFieldRef();
#else
    scalarField& muEffI = muEff_.internalField();
#endif

    forAll(muEffI, cellI)
    {
        scalarField mu(3, 0);

        symmTensor devSigma = symmTensor::zero;
        calcDevCauchy(FI[cellI], EcI[cellI], EaI[cellI], ErI[cellI], devSigma);

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
            const tensor pertGradDDT = pertGradDD.T();

            tensor DF = pertGradDDT;
            tensor newF = FI[cellI] + DF;

            symmTensor newDevSigma = symmTensor::zero;
            calcDevCauchy(newF, EcI[cellI], EaI[cellI], ErI[cellI], newDevSigma);

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
            pertGradDD = (pVectorsT & pertGradDD & pVectors);
            const tensor pertGradDDT = pertGradDD.T();

            tensor DF = pertGradDDT;
            tensor newFI = FI[cellI] + DF;

            symmTensor newDevSigma = symmTensor::zero;
            calcDevCauchy(newFI, EcI[cellI], EaI[cellI], ErI[cellI], newDevSigma);

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
            const tensor pertGradDDT = pertGradDD.T();

            tensor DF = pertGradDDT;
            tensor newFI = FI[cellI] + DF;

            symmTensor newDevSigma = symmTensor::zero;
            calcDevCauchy(newFI, EcI[cellI], EaI[cellI], ErI[cellI], newDevSigma);

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


void Foam::HolzapfelGasserOgdenElastic::updateTotalFields()
{
    calcEffectiveShearModulus();
}

// ************************************************************************* //
