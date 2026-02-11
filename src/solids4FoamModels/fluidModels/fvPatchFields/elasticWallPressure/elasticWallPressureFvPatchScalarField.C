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

#include "elasticWallPressureFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "fluidSolidInterface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Functions * * * * * * * * * * * * * * * //

const scalarField& elasticWallPressureFvPatchScalarField::rhoSolidHs() const
{
    if (rhoSolidHsPtr_.empty())
    {
    #ifdef OPENFOAM_NOT_EXTEND
        const fvMesh& mesh = internalField().mesh();
    #else
        const fvMesh& mesh = dimensionedInternalField().mesh();
    #endif

        // Looking up the FSI solver
        const fluidSolidInterface& fsi =
            mesh.objectRegistry::parent().lookupObject<fluidSolidInterface>
            (
                "fsiProperties"
            );

        // Find the solid patch ID corresponding to the current fluid patch
        label solidPatchID = -1;
        label interfaceID = -1;
        forAll(fsi.fluidPatchIndices(), i)
        {
            if (fsi.fluidPatchIndices()[i] == patch().index())
            {
                // Take the corresponding solid patch ID
                solidPatchID = fsi.solidPatchIndices()[i];
                interfaceID = i;
                break;
            }
        }

        if (solidPatchID == -1)
        {
            FatalErrorInFunction
                << "Are you sure this patch is an FSI interface?"
                << abort(FatalError);
        }

        // Get solid density
        const scalarField& rho =
            fsi.solidMesh().lookupObject<volScalarField>
            (
                "rho"
            ).boundaryField()[solidPatchID];

        // Get solid stiffness (impK for generality)
        const scalarField& impK =
            fsi.solidMesh().lookupObject<volScalarField>
            (
                "impK"
            ).boundaryField()[solidPatchID];

        // p-wave propagation speed, ap, on the solid patch
        const scalarField ap(sqrt(impK/rho));

        // Solid "virtual thickness"
        scalarField hs(rho.size(), constantHs_);
        if (constantHs_ < SMALL)
        {
            // Calculate a virtual thickness based on the speed of sound and time
            // step
            hs = ap*mesh.time().deltaT().value();

            if (debug)
            {
                Info<< "hs: min = " << min(hs) << ", max = " << max(hs)
                    << ", mean = " << average(hs) << endl;
            }
        }

        // Calculate rhoHs at the solid
        const scalarField rhoHs(rho*hs);

        // Take references to zones
        const standAlonePatch& fluidZone =
            fsi.fluid().globalPatches()[interfaceID].globalPatch();
        const standAlonePatch& solidZone =
            fsi.solid().globalPatches()[interfaceID].globalPatch();

        // Map the solid patch field to the global zone
        scalarField rhoHsZone
        (
            fsi.solid().globalPatches()[interfaceID].patchFaceToGlobal(rhoHs)
        );
        scalarField rhoHsZoneAtFluid(fluidZone.size(), 0.0);

        // Map rhoHs from the solid patch to the fluid patch
        // Transfer the field from the solid interface to the fluid interface
        fsi.interfaceToInterfaceList()[interfaceID].transferFacesZoneToZone
        (
            solidZone,          // from zone
            fluidZone,          // to zone
            rhoHsZone,          // from field
            rhoHsZoneAtFluid    // to field
        );

        // Initialise the rhoHs field for the fluid patch and map zone to patch
        rhoSolidHsPtr_.set
        (
            new scalarField
            (
                fsi.fluid().globalPatches()
                [
                    interfaceID
                ].globalFaceToPatch(rhoHsZoneAtFluid)
            )
        );
    }

    return rhoSolidHsPtr_();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

elasticWallPressureFvPatchScalarField::elasticWallPressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    robinFvPatchScalarField(p, iF),
    prevPressure_(p.patch().size(), 0),
    prevAcceleration_(p.patch().size(), vector::zero),
    rhoSolidHsPtr_(),
    constantHs_(-1.0)
{}


elasticWallPressureFvPatchScalarField::elasticWallPressureFvPatchScalarField
(
    const elasticWallPressureFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    robinFvPatchScalarField(ptf, p, iF, mapper),
    prevPressure_(p.patch().size(), 0),
    prevAcceleration_(p.patch().size(), vector::zero),
    rhoSolidHsPtr_(),
    constantHs_(ptf.constantHs_)
{}


elasticWallPressureFvPatchScalarField::elasticWallPressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    robinFvPatchScalarField(p, iF),
    prevPressure_(p.patch().size(), 0),
    prevAcceleration_(p.patch().size(), vector::zero),
    rhoSolidHsPtr_(),
    constantHs_(dict.lookupOrDefault<scalar>("constantHs", -1.0))
{
    if (dict.found("value"))
    {
        Field<scalar>::operator=(scalarField("value", dict, p.size()));
    }

    if (dict.found("prevPressure"))
    {
        Field<scalar>::operator=(scalarField("prevPressure", dict, p.size()));
    }

    if (constantHs_ < SMALL)
    {
        Info<< type() << " " << patch().name() << ": constantHs unused" << endl;
    }
    else
    {
        Info<< type() << " " << patch().name() << ": constantHs = "
            << constantHs_ << endl;
    }

    this->coeff0() = 1.0;
    this->coeff1() = 1.0;
}


#ifndef OPENFOAM_ORG
elasticWallPressureFvPatchScalarField::elasticWallPressureFvPatchScalarField
(
    const elasticWallPressureFvPatchScalarField& pivpvf
)
:
    robinFvPatchScalarField(pivpvf),
    prevPressure_(pivpvf.prevPressure_),
    prevAcceleration_(pivpvf.prevAcceleration_),
    constantHs_(pivpvf.constantHs_)
{}
#endif


elasticWallPressureFvPatchScalarField::elasticWallPressureFvPatchScalarField
(
    const elasticWallPressureFvPatchScalarField& pivpvf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    robinFvPatchScalarField(pivpvf, iF),
    prevPressure_(pivpvf.prevPressure_),
    prevAcceleration_(pivpvf.prevAcceleration_),
    rhoSolidHsPtr_(),
    constantHs_(pivpvf.constantHs_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void elasticWallPressureFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fvPatchField<scalar>::autoMap(m);
}


void elasticWallPressureFvPatchScalarField::rmap
(
    const fvPatchField<scalar>& ptf,
    const labelList& addr
)
{
    fvPatchField<scalar>::rmap(ptf, addr);
}

void elasticWallPressureFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

#ifdef OPENFOAM_NOT_EXTEND
    const fvMesh& mesh = internalField().mesh();
#else
    const fvMesh& mesh = dimensionedInternalField().mesh();
#endif

    // Looking up the FSI solver
    const fluidSolidInterface& fsi =
        mesh.objectRegistry::parent().lookupObject<fluidSolidInterface>
        (
            "fsiProperties"
        );

    // Map the solid density times the solid virtual thickness mapped to the
    // current (fluid) patch
    const scalarField& rhoSolidHs = this->rhoSolidHs();

    // Lookup the pressure dimensions
#ifdef OPENFOAM_NOT_EXTEND
    const word fieldName = internalField().name();
#else
    const word fieldName = dimensionedInternalField().name();
#endif

    const dimensionSet& pDims =
        mesh.lookupObject<volScalarField>(fieldName).dimensions();;

    // The previous acceleration is updated at the end of each
    // time step in the fluidSolidInterface
    // const vectorField n(p.nf());
    const scalarField prevDdtUn(patch().nf() & prevAcceleration_);

    // Check if a density field is present
    if (fsi.fluid().mesh().foundObject<volScalarField>("rho"))
    {
        const scalarField& rhoFluid =
            patch().lookupPatchField<volScalarField, scalar>("rho");
        const scalarField& phig =
            patch().lookupPatchField<surfaceScalarField, scalar>("phig");

        const scalarField c1(rhoSolidHs/rhoFluid);

        if (pDims == dimPressure/dimDensity)
        {
            // p/rho
            // Divide RHS by rhoFluid
            this->coeff0() = 1.0;
            this->coeff1() = c1;
            this->rhs() =
                (prevPressure_ - rhoSolidHs*prevDdtUn + c1*phig)/rhoFluid;
        }
        else
        {
            // p
            this->coeff0() = 1.0;
            this->coeff1() = c1;
            this->rhs() = prevPressure_ - rhoSolidHs*prevDdtUn + c1*phig;
        }
    }
    else
    {
        if (debug)
        {
            Info<< "Did not find rho: looking up from transportProperties"
                << endl;
        }

        // Fluid properties
        const dictionary& transportProperties =
            db().lookupObject<IOdictionary>("transportProperties");

        // Lookup the density from the transport properties
        const dimensionedScalar rhoFluid
        (
            transportProperties.lookup("rho")
        );

        if (debug)
        {
            Info<< "rhoSolidHs = " << max(rhoSolidHs)
                << ", rhoFluid = " << rhoFluid.value()
                << endl;
        }

        if (pDims == dimPressure/dimDensity)
        {
            // p/rho
            this->coeff0() = 1.0;
            this->coeff1() = rhoSolidHs/rhoFluid.value();
            this->rhs() =
                prevPressure_/rhoFluid.value()
              - rhoSolidHs*prevDdtUn/rhoFluid.value();
        }
        else
        {
            // p
            this->coeff0() = 1.0;
            this->coeff1() = rhoSolidHs/rhoFluid.value();
            this->rhs() = prevPressure_ - rhoSolidHs*prevDdtUn;
        }
    }

    robinFvPatchField<scalar>::updateCoeffs();
}


void elasticWallPressureFvPatchScalarField::patchFlux
(
    GeometricField<scalar, fvsPatchField, surfaceMesh>& flux,
    const fvMatrix<scalar>& matrix
) const
{
    scalarField rAU(patch().size(), 0.0);
    if (db().foundObject<volScalarField>("rAU"))
    {
        rAU = patch().lookupPatchField<volScalarField, scalar>("rAU");
    }
    else
    {
        rAU = patch().lookupPatchField<surfaceScalarField, scalar>("rAUf");
    }

#ifdef OPENFOAM_NOT_EXTEND
    flux.boundaryFieldRef()[patch().index()] = rAU*snGrad()*patch().magSf();
#else
    flux.boundaryField()[patch().index()] = rAU*snGrad()*patch().magSf();
#endif
}


void elasticWallPressureFvPatchScalarField::write(Ostream& os) const
{
    robinFvPatchScalarField::write(os);
#ifdef OPENFOAM_ORG
    writeEntry(os, "prevPressure", prevPressure_);
#else
    prevPressure_.writeEntry("prevPressure", os);
#endif
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    elasticWallPressureFvPatchScalarField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
