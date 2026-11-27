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

#include "wallNormalMomentum.H"
#include "surfaceFields.H"
#include "fvc.H"
#include "compatibilityFunctions.H"
#include "linear.H"
#ifdef OPENFOAM_COM
    #include "fvOptions.H"
    #include "fvMatrices.H"
#endif
#ifdef OPENFOAM_ORG
    #include "momentumTransportModel.H"
#else
    #include "turbulenceModel.H"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace wallNormalMomentum
{

// * * * * * * * * * * * * * * * * Functions * * * * * * * * * * * * * * * * //

tmp<volVectorField> field
(
    const objectRegistry& db,
    const word& UName,
    const word& patchName
)
{
#ifdef OPENFOAM_ORG
    typedef momentumTransportModel viscosityModelType;
#elif defined(OPENFOAM_COM)
    typedef turbulenceModel viscosityModelType;
#else
    typedef incompressible::turbulenceModel viscosityModelType;
#endif

    const wordList modelNames(db.lookupClass<viscosityModelType>().toc());

    if (modelNames.size() != 1)
    {
        FatalErrorInFunction
            << "Patch " << patchName << ": the viscous term needs one "
            << "turbulence (momentum transport) model in the fluid region, "
            << "found " << modelNames.size() << nl
            << "Set viscousTerm no to omit the viscous term"
            << abort(FatalError);
    }

    const viscosityModelType& model =
        db.lookupObject<viscosityModelType>(modelNames[0]);

    const volVectorField& U = db.lookupObject<volVectorField>(UName);

    // Named intermediates: a reused tmp would alias an input of the products
    const volScalarField nuEff(model.nuEff());
    const volTensorField gradU(fvc::grad(U));
    const volTensorField gradUT(gradU.T());
    const volTensorField devGradUT(dev2(gradUT));

    // Variable density (a registered rho field, e.g. interFluid): the
    // momentum equation has the stress divergence of rho*nuEff, which is
    // divided by rho to give an acceleration. Such solvers do not apply the
    // fvOptions to the momentum equation, so the sources are not included
    const bool variableDensity = db.foundObject<volScalarField>("rho");

#ifndef OPENFOAM_COM
    // Only the OpenFOAM.com interFluid adds the corresponding interface flux
    // and uses this form of the viscous term: elsewhere the term is omitted
    // for variable density
    if (variableDensity)
    {
        return tmp<volVectorField>
        (
            new volVectorField
            (
                IOobject
                (
                    "wallNormalMomentum",
                    U.time().timeName(),
                    U.mesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    false
                ),
                U.mesh(),
                dimensionedVector("zero", dimAcceleration, vector::zero)
            )
        );
    }
#endif

    const volScalarField muEff
    (
        variableDensity
      ? volScalarField(db.lookupObject<volScalarField>("rho")*nuEff)
      : nuEff
    );

    const volTensorField muDevGradUT(muEff*devGradUT);

    // The transpose part is the Gauss linear divergence, as usually in the
    // fluid momentum equation, written with linearly interpolated face values
    // so that it needs no divergence or interpolation scheme of its own
    const surfaceVectorField muDevGradUTf
    (
        U.mesh().Sf() & linearInterpolate(muDevGradUT)
    );

    tmp<volVectorField> tterm
    (
        fvc::laplacian(muEff, U) + fvc::div(muDevGradUTf)
    );

    if (variableDensity)
    {
        tmpRef(tterm) /= db.lookupObject<volScalarField>("rho");
    }

#ifdef OPENFOAM_COM
    // Momentum sources (fvOptions), e.g. a body force: "M & U" evaluates the
    // source matrix (Su + Sp*U) per unit volume
    fv::options& fvOptions = fv::options::New(U.mesh());

    if (!variableDensity && fvOptions.optionList::size())
    {
        volVectorField& Uref = const_cast<volVectorField&>(U);
        const tmp<fvVectorMatrix> tM(fvOptions(Uref));
        const volVectorField source(tM() & Uref);

        tmpRef(tterm) += source;
    }
#endif

    return tterm;
}


tmp<scalarField> patchNormal
(
    const fvPatch& patch,
    const volVectorField& term
)
{
    const vectorField n(patch.nf());
    const vectorField viscousP
    (
        term.boundaryField()[patch.index()].patchInternalField()
    );

    return tmp<scalarField>(new scalarField(n & viscousP));
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace wallNormalMomentum

} // End namespace Foam

// ************************************************************************* //
