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

#include "diffusionElastic.H"
#include "addToRunTimeSelectionTable.H"
#include "fvc.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(diffusionElastic, 0);
    addToRunTimeSelectionTable
    (
        mechanicalLaw, diffusionElastic, linGeomMechLaw
    );
}


// * * * * * * * * * * * * * Private Member Functions * * * * * * * * * * * * //


void Foam::diffusionElastic::updateD() const
{
    // Recalculate the distances
    motionDiffPtr_->correct();

    // Limit the range and normalise wrt to the average
    df_ = motionDiffPtr_()();
    const dimensionedScalar avDf = average(df_);
    df_ = min(maxFactor_*avDf, max(minFactor_*avDf, df_));
    df_ /= avDf;
    d_ = fvc::average(df_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from dictionary
Foam::diffusionElastic::diffusionElastic
(
    const word& name,
    const fvMesh& mesh,
    const dictionary& dict,
    const nonLinearGeometry::nonLinearType& nonLinGeom
)
:
    mechanicalLaw(name, mesh, dict, nonLinGeom),
    mu_(dict.lookup("mu")),
    kappa_(dict.lookup("kappa")),
    maxFactor_(dict.lookupOrDefault<scalar>("maxFactor", 10)),
    minFactor_(dict.lookupOrDefault<scalar>("minFactor", 0.1)),
    motionDiffPtr_(motionDiffusivity::New(mesh, dict.lookup("diffusivity"))),
    df_("distf", motionDiffPtr_()()),
    d_("dist", fvc::average(df_))
{
    // Store old times
    epsilon().storeOldTime();
    epsilonf().storeOldTime();

    if (maxFactor_ < minFactor_)
    {
        FatalErrorInFunction
            << "maxFactor cannot be less than minFactor!" << exit(FatalError);
    }

    updateD();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


Foam::tmp<Foam::volScalarField> Foam::diffusionElastic::bulkModulus() const
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
            kappa_*d_
        )
    );
}


Foam::tmp<Foam::volScalarField> Foam::diffusionElastic::shearModulus() const
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


Foam::tmp<Foam::volScalarField> Foam::diffusionElastic::impK() const
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
            d_*((4.0/3.0)*mu_ + kappa_)
        )
    );
}


void Foam::diffusionElastic::correct(volSymmTensorField& sigma)
{
    // Update epsilon
    updateEpsilon();

    sigma = d_*(2*mu_*dev(epsilon()) + kappa_*tr(epsilon())*I);
}


void Foam::diffusionElastic::correct(surfaceSymmTensorField& sigma)
{
    // Update epsilon
    updateEpsilonf();

    sigma = df_*(2*mu_*dev(epsilonf()) + kappa_*tr(epsilonf())*I);
}


// ************************************************************************* //
