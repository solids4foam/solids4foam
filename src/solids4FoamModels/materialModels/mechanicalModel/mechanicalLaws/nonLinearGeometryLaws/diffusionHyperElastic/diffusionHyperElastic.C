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

#ifdef OPENFOAM_NOT_EXTEND

#include "diffusionHyperElastic.H"
#include "addToRunTimeSelectionTable.H"
#include "fvc.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(diffusionHyperElastic, 0);
    addToRunTimeSelectionTable
    (
        mechanicalLaw, diffusionHyperElastic, nonLinGeomMechLaw
    );
}


// * * * * * * * * * * * * * Private Member Functions * * * * * * * * * * * * //


void Foam::diffusionHyperElastic::updateD() const
{
    // Recalculate the distances
    motionDiffPtr_->correct();

    // Limit the range and normalise wrt to the average
    df_ = motionDiffPtr_()();
    const dimensionedScalar avDf = average(df_);
    df_ = min(maxFactor_*avDf, max(minFactor_*avDf, df_));
    df_ /= avDf;
    d_ = fvc::average(df_);

    // Enforce zero gradient boundaries
    d_.correctBoundaryConditions();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from dictionary
Foam::diffusionHyperElastic::diffusionHyperElastic
(
    const word& name,
    const fvMesh& mesh,
    const dictionary& dict,
    const nonLinearGeometry::nonLinearType& nonLinGeom
)
:
    mechanicalLaw(name, mesh, dict, nonLinGeom),
    mu_(dict.lookup("mu")),
    K_(dict.lookup("K")),
    maxFactor_(dict.lookupOrDefault<scalar>("maxFactor", 100)),
    minFactor_(dict.lookupOrDefault<scalar>("minFactor", 0.1)),
    motionDiffPtr_(motionDiffusivity::New(mesh, dict.lookup("diffusivity"))),
    df_("distf", motionDiffPtr_()()),
    d_
    (
        IOobject
        (
            "stiffnessScaleFactor",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        fvc::average(df_),
        "zeroGradient"
    )
{
    // Store old F
    F().storeOldTime();
    // Ff().storeOldTime();

    if (maxFactor_ < minFactor_)
    {
        FatalErrorInFunction
            << "maxFactor cannot be less than minFactor!" << exit(FatalError);
    }

    // Enforce zero gradient boundaries
    d_.correctBoundaryConditions();

    // Update D
    updateD();

    if (Switch(dict.lookup("writeStiffScaleFactor")))
    {
        d_.writeOpt() = IOobject::AUTO_WRITE;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


Foam::tmp<Foam::volScalarField> Foam::diffusionHyperElastic::bulkModulus() const
{
    // updateD();

    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "bulkModulus",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            d_*K_
        )
    );
}


Foam::tmp<Foam::volScalarField> Foam::diffusionHyperElastic::shearModulus() const
{
    // updateD();

    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "shearModulus",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            d_*mu_
        )
    );
}


Foam::tmp<Foam::volScalarField> Foam::diffusionHyperElastic::impK() const
{
    // updateD();

    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "k",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            d_*((4.0/3.0)*mu_ + K_)
        )
    );
}


void Foam::diffusionHyperElastic::correct(volSymmTensorField& sigma)
{
    // Update the deformation gradient field
    // Note: if true is returned, it means that linearised elasticity was
    // enforced by the solver via the enforceLinear switch
    if (updateF(sigma, mu_, K_))
    {
        return;
    }

    // Calculate the Jacobian of the deformation gradient
    const volScalarField J(det(F()));

    // Calculate the volume preserving left Cauchy Green strain
    const volSymmTensorField bEbar(pow(J, -2.0/3.0)*symm(F() & F().T()));

    // Calculate the deviatoric stress
    const volSymmTensorField s(mu_*dev(bEbar));

    // Calculate the Cauchy stress
    sigma = d_*(s + 0.5*K()*(pow(J, 2.0) - 1.0)*I)/J;
}


void Foam::diffusionHyperElastic::correct(surfaceSymmTensorField& sigma)
{
    NotImplemented;
}

#endif // OPENFOAM_NOT_EXTEND

// ************************************************************************* //
