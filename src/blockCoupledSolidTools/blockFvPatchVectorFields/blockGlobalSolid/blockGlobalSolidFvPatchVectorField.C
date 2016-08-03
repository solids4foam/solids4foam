/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "blockGlobalSolidFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "fvcMeshPhi.H"
#include "solidPolyMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

blockGlobalSolidFvPatchVectorField::
blockGlobalSolidFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    coupledFvPatchField<vector>(p, iF),
    blockFvPatchVectorField()
{}


blockGlobalSolidFvPatchVectorField::
blockGlobalSolidFvPatchVectorField
(
    const blockGlobalSolidFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    coupledFvPatchField<vector>(ptf, p, iF, mapper),
    blockFvPatchVectorField()
{}


blockGlobalSolidFvPatchVectorField::
blockGlobalSolidFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    //coupledFvPatchField<vector>(p, iF, dict),
    coupledFvPatchField<vector>(p, iF),
    blockFvPatchVectorField()
{
    fvPatchField<vector>::operator=
    (
        vectorField(p.size(), vector::zero)
    );
}


blockGlobalSolidFvPatchVectorField::
blockGlobalSolidFvPatchVectorField
(
    const blockGlobalSolidFvPatchVectorField& ptf
)
:
    coupledFvPatchField<vector>(ptf),
    blockFvPatchVectorField()
{}


blockGlobalSolidFvPatchVectorField::
blockGlobalSolidFvPatchVectorField
(
    const blockGlobalSolidFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    coupledFvPatchField<vector>(ptf, iF),
    blockFvPatchVectorField()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp< Foam::Field<vector> >
blockGlobalSolidFvPatchVectorField::patchNeighbourField() const
{
    if (patch().size())
    {
        FatalErrorIn("patchNeighbourField()")
            << "The global patch should not have any faces!"
            << abort(FatalError);
    }

    // keep compiler happy
    return patchInternalField();
}


void blockGlobalSolidFvPatchVectorField::transformCoupleField
(
    scalarField &f, const direction cmpt
) const
{
    FatalErrorIn("transformCoupleField()")
        << "This function should not be called"
        << abort(FatalError);
}


Foam::tmp<Foam::Field<vector> >
blockGlobalSolidFvPatchVectorField::snGrad() const
{
    if (patch().size())
    {
        FatalErrorIn("snGrad()")
            << "The global patch should not have any faces!"
            << abort(FatalError);
    }

    // keep compiler happy
    return patchInternalField();
}


tmp<Field<vector> > blockGlobalSolidFvPatchVectorField::
gradientBoundaryCoeffs() const
{
    FatalErrorIn("gradientBoundaryCoeffs()")
        << "This function should not be called!" << nl
        << "This boundary condition is only for use with the block coupled"
        << " solid solver"
        << abort(FatalError);

    // Keep the compiler happy
    return *this;
}


void blockGlobalSolidFvPatchVectorField::updateInterfaceMatrix
(
    const scalarField &psiInternal,
    scalarField &result,
    const lduMatrix &,
    const scalarField &coeffs,
    const direction,
    const Pstream::commsTypes commsType,
    const bool switchToLhs
) const
{
    FatalErrorIn("updateInterfaceMatrix()")
        << "This function should not be called!"
        << abort(FatalError);
}


void blockGlobalSolidFvPatchVectorField::updateInterfaceMatrix
(
    const Field< vector >& x,
    Field< vector >& Ax,
    const BlockLduMatrix< vector >& matrix,
    const CoeffField< vector >& coeffs,
    const Pstream::commsTypes commsType,
    const bool switchToLhs
) const
{
    // disabled
    // if (Pstream::parRun())
    // {
    //     patch().boundaryMesh().mesh().thisDb().lookupObject<solidPolyMesh>
    //     (
    //         "solidPolyMesh"
    //     ).updateGlobalFields(Ax, x);
    // }
}


void blockGlobalSolidFvPatchVectorField::insertBlockCoeffs
(
    const solidPolyMesh& solidMesh,
    const surfaceScalarField& muf,
    const surfaceScalarField& lambdaf,
    const GeometricField<vector, fvPatchField, volMesh>& U,
    Field<vector>& blockB,
    BlockLduMatrix<vector>& blockM
) const
{
    if (patch().size())
    {
        FatalErrorIn("insertBlockCoeffs()")
            << "The global patch should not have any faces!"
            << abort(FatalError);
    }
}

void blockGlobalSolidFvPatchVectorField::write(Ostream& os) const
{
    coupledFvPatchField<vector>::write(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchVectorField,
    blockGlobalSolidFvPatchVectorField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
