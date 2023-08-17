/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "thermalRobinFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "fluidSolidInterface.H"
#include "fvcMeshPhi.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

thermalRobinFvPatchScalarField::
thermalRobinFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    robinFvPatchScalarField(p, iF),
    neumann_(true),
    dirichlet_(true),
    neiTemperature_(p.patch().size(), 0),
    neiHeatFlux_(p.patch().size(), 0),
    eqInterHeatTransferCoeff_(p.patch().size(), 0),
    lambdaName_("lambda")
{}


thermalRobinFvPatchScalarField::
thermalRobinFvPatchScalarField
(
    const thermalRobinFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    robinFvPatchScalarField(ptf, p, iF, mapper),
    neumann_(ptf.neumann_),
    dirichlet_(ptf.dirichlet_),
#ifdef OPENFOAMFOUNDATION
    neiTemperature_(mapper(ptf.neiTemperature_)),
    neiHeatFlux_(mapper(ptf.neiHeatFlux_)),
    eqInterHeatTransferCoeff_(mapper(ptf.eqInterHeatTransferCoeff_)),
#else
    neiTemperature_(ptf.neiTemperature_, mapper),
    neiHeatFlux_(ptf.neiHeatFlux_, mapper),
    eqInterHeatTransferCoeff_(ptf.eqInterHeatTransferCoeff_, mapper),
#endif
    lambdaName_(ptf.lambdaName_)
{}


thermalRobinFvPatchScalarField::
thermalRobinFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    robinFvPatchScalarField(p, iF),
    neumann_(dict.lookupOrDefault<Switch>("neumann", true)),
    dirichlet_(dict.lookupOrDefault<Switch>("dirichlet", true)),
    neiTemperature_(p.patch().size(), 0),
    neiHeatFlux_(p.patch().size(), 0),
    eqInterHeatTransferCoeff_(p.patch().size(), 1),
    lambdaName_(dict.lookupOrDefault<word>("lambda", "lambda"))
{
    if (dict.found("neiTemperature"))
    {
        neiTemperature_ =
            scalarField("neiTemperature", dict, p.size());
    }
    else
    {
        neiTemperature_ =
            patchInternalField();
    }
    
    if (dict.found("neiHeatFlux"))
    {
        neiHeatFlux_ =
            scalarField("neiHeatFlux", dict, p.size());
    }
    
    if (dict.found("eqInterHeatTransferCoeff"))
    {
        eqInterHeatTransferCoeff_ =
            scalarField("eqInterHeatTransferCoeff", dict, p.size());
    }

    Field<scalar>::operator=(patchInternalField());
}


#ifndef OPENFOAMFOUNDATION
thermalRobinFvPatchScalarField::
thermalRobinFvPatchScalarField
(
    const thermalRobinFvPatchScalarField& pivpvf
)
:
    robinFvPatchScalarField(pivpvf),
    neumann_(pivpvf.neumann_),
    dirichlet_(pivpvf.dirichlet_),
    neiTemperature_(pivpvf.neiTemperature_),
    neiHeatFlux_(pivpvf.neiHeatFlux_),
    eqInterHeatTransferCoeff_(pivpvf.eqInterHeatTransferCoeff_),
    lambdaName_(pivpvf.lambdaName_)
{}
#endif


thermalRobinFvPatchScalarField::
thermalRobinFvPatchScalarField
(
    const thermalRobinFvPatchScalarField& pivpvf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    robinFvPatchScalarField(pivpvf, iF),
    neumann_(pivpvf.neumann_),
    dirichlet_(pivpvf.dirichlet_),
    neiTemperature_(pivpvf.neiTemperature_),
    neiHeatFlux_(pivpvf.neiHeatFlux_),
    eqInterHeatTransferCoeff_(pivpvf.eqInterHeatTransferCoeff_),
    lambdaName_(pivpvf.lambdaName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void thermalRobinFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fvPatchField<scalar>::autoMap(m);
}


void thermalRobinFvPatchScalarField::rmap
(
    const fvPatchField<scalar>& ptf,
    const labelList& addr
)
{
    fvPatchField<scalar>::rmap(ptf, addr);
}


void thermalRobinFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const scalarField lambda =
        patch().lookupPatchField<volScalarField, scalar>(lambdaName_);

    if (dirichlet_ && neumann_)
    {
        this->coeff0() = 1.0;
        this->coeff1() = eqInterHeatTransferCoeff_*lambda;
        this->rhs() = neiTemperature_ + eqInterHeatTransferCoeff_*neiHeatFlux_;
    }
    else if (dirichlet_ && !neumann_)
    {
        this->coeff0() = 1.0;
        this->coeff1() = 0.0;
        this->rhs() = neiTemperature_;
    }
    else if (!dirichlet_ && neumann_)
    {
        this->coeff0() = 0.0;
        this->coeff1() = lambda;
        this->rhs() = neiHeatFlux_;
    }

    robinFvPatchScalarField::updateCoeffs();
}


void thermalRobinFvPatchScalarField::write(Ostream& os) const
{
    robinFvPatchScalarField::write(os);

#ifdef OPENFOAMESI
    os.writeEntryIfDifferent<Switch>("neumann", true, neumann_);
    os.writeEntryIfDifferent<Switch>("dirichlet", true, dirichlet_);
#else
    writeEntryIfDifferent<Switch>(os, "neumann", true, neumann_);
    writeEntryIfDifferent<Switch>(os, "dirichlet", true, dirichlet_);
#endif

    
#ifdef OPENFOAMFOUNDATION
    writeEntry(os, "neiTemperature", neiTemperature_);
    writeEntry(os, "neiHeatFlux", neiHeatFlux_);
    writeEntry(os, "eqInterHeatTransferCoeff", eqInterHeatTransferCoeff_);
#else
    neiTemperature_.writeEntry("neiTemperature", os);
    neiHeatFlux_.writeEntry("neiHeatFlux", os);
    eqInterHeatTransferCoeff_.writeEntry("eqInterHeatTransferCoeff", os);
#endif

    if (lambdaName_ != "lambda")
    {
        os.writeKeyword("lambda") << lambdaName_ << token::END_STATEMENT << nl;
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    thermalRobinFvPatchScalarField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
