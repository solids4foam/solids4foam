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

#include "wallNormalViscousTerm.H"
#include "surfaceFields.H"
#include "fvc.H"
#include "linear.H"
#ifdef OPENFOAM_ORG
    #include "momentumTransportModel.H"
#else
    #include "turbulenceModel.H"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace wallNormalViscousTerm
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
    const volTensorField nuDevGradUT(nuEff*dev2(gradUT));

    // The transpose part is the Gauss linear divergence, as usually in the
    // fluid momentum equation, written with linearly interpolated face values
    // so that it needs no divergence or interpolation scheme of its own
    const surfaceVectorField nuDevGradUTf
    (
        U.mesh().Sf() & linearInterpolate(nuDevGradUT)
    );

    return fvc::laplacian(nuEff, U) + fvc::div(nuDevGradUTf);
}


tmp<scalarField> patchNormal
(
    const fvPatch& patch,
    const volVectorField& viscousTerm
)
{
    const vectorField n(patch.nf());
    const vectorField viscousP
    (
        viscousTerm.boundaryField()[patch.index()].patchInternalField()
    );

    return tmp<scalarField>(new scalarField(n & viscousP));
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace wallNormalViscousTerm

} // End namespace Foam

// ************************************************************************* //
