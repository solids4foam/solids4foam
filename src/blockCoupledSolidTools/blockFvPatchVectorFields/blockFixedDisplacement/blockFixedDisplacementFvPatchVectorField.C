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

#include "blockFixedDisplacementFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "fvcMeshPhi.H"
#include "fixedValuePointPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

blockFixedDisplacementFvPatchVectorField::
blockFixedDisplacementFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF),
    blockFvPatchVectorField(),
    totalDisp_(vector::zero),
    dispSeries_()
{}


blockFixedDisplacementFvPatchVectorField::
blockFixedDisplacementFvPatchVectorField
(
    const blockFixedDisplacementFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper),
    blockFvPatchVectorField(),
    totalDisp_(ptf.totalDisp_),
    dispSeries_(ptf.dispSeries_)
{}


blockFixedDisplacementFvPatchVectorField::
blockFixedDisplacementFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF, dict),
    blockFvPatchVectorField(),
    totalDisp_(vector::zero),
    dispSeries_()
{
    // Check if displacement is time-varying
    if (dict.found("displacementSeries"))
    {
        Info<< "    displacement is time-varying" << endl;
        dispSeries_ =
            interpolationTable<vector>(dict.subDict("displacementSeries"));

        fvPatchField<vector>::operator==
        (
            dispSeries_(this->db().time().timeOutputValue())
        );
    }
    else
    {
        if (dict.found("value"))
        {
            vectorField val = vectorField("value", dict, p.size());

            totalDisp_ = gAverage(val);
        }
    }
}

blockFixedDisplacementFvPatchVectorField::
blockFixedDisplacementFvPatchVectorField
(
    const blockFixedDisplacementFvPatchVectorField& ptf
)
:
    fixedValueFvPatchVectorField(ptf),
    blockFvPatchVectorField(),
    totalDisp_(ptf.totalDisp_),
    dispSeries_(ptf.dispSeries_)
{}


blockFixedDisplacementFvPatchVectorField::
blockFixedDisplacementFvPatchVectorField
(
    const blockFixedDisplacementFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(ptf, iF),
    blockFvPatchVectorField(),
    totalDisp_(ptf.totalDisp_),
    dispSeries_(ptf.dispSeries_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void blockFixedDisplacementFvPatchVectorField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    vectorField disp(patch().size(), totalDisp_);
    vectorField pointDisp(patch().patch().nPoints(), totalDisp_);

    if (dispSeries_.size())
    {
        disp = dispSeries_(this->db().time().timeOutputValue());
        pointDisp = dispSeries_(this->db().time().timeOutputValue());
    }

    if (dimensionedInternalField().name() == "DD")
    {
        // Incremental approach, so we wil set the increment of displacement
        // Lookup the old displacement field and subtract it from the total
        // displacement
        const volVectorField& Dold =
            db().lookupObject<volVectorField>("D").oldTime();

        disp -= Dold.boundaryField()[patch().index()];
    }

    fvPatchField<vector>::operator==(disp);


    // Update point displacement field

    const fvMesh& mesh = patch().boundaryMesh().mesh();

    if
    (
        mesh.lookupObject<pointVectorField>
        (
            "point" + dimensionedInternalField().name()
        ).boundaryField()[patch().index()].type() == "fixedValue"
    )
    {
        if (dimensionedInternalField().name() == "DD")
        {
            const pointVectorField& pointD =
                mesh.objectRegistry::lookupObject
                <
                pointVectorField
                >("pointD");

            const pointVectorField& pointDOld = pointD.oldTime();

            const labelList& meshPoints = patch().patch().meshPoints();

            forAll(meshPoints, pI)
            {
                const label pointID = meshPoints[pI];

                pointDisp[pI] -= pointDOld[pointID];
            }
        }

        fixedValuePointPatchVectorField& patchPointDD =
            refCast<fixedValuePointPatchVectorField>
            (
                const_cast<pointVectorField&>
                (
                    mesh.objectRegistry::lookupObject<pointVectorField>
                    ("point" + dimensionedInternalField().name())
                ).boundaryField()[patch().index()]
            );

        patchPointDD == pointDisp;
    }

    fixedValueFvPatchVectorField::updateCoeffs();
}


Foam::tmp<Foam::Field<vector> >
blockFixedDisplacementFvPatchVectorField::snGrad() const
{
    // snGrad without non-orthogonal correction
    // return (*this - patchInternalField())*this->patch().deltaCoeffs();

    // Lookup previous boundary gradient
    const fvPatchField<tensor>& gradField =
        patch().lookupPatchField<volTensorField, tensor>
        (
            "grad(" + dimensionedInternalField().name() + ")"
        );

    // Calculate correction vector
    vectorField n = this->patch().nf();
    vectorField delta = this->patch().delta();
    vectorField k = ((I - sqr(n)) & delta);

    return
    (
        //*this - patchInternalField()
        //*this - (patchInternalField() + (k & gradField.patchInternalField()))
        //*this - (patchInternalField() + (k & gradField))
        (*this - (k & gradField)) - patchInternalField()
    )*this->patch().deltaCoeffs();
}


tmp<Field<vector> > blockFixedDisplacementFvPatchVectorField::
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


void blockFixedDisplacementFvPatchVectorField::insertBlockCoeffs
(
    const solidPolyMesh& solidMesh,
    const surfaceScalarField& muf,
    const surfaceScalarField& lambdaf,
    const GeometricField<vector, fvPatchField, volMesh>& U,
    Field<vector>& blockB,
    BlockLduMatrix<vector>& blockM
) const
{
    // Const reference to polyPatch and the fvMesh
    const polyPatch& ppatch = patch().patch();

    // Update the displacement
    // We shouldn't have to use const_cast ...
    const_cast<blockFixedDisplacementFvPatchVectorField&>(*this).updateCoeffs();

    // Const reference to the patch field
    const vectorField& pU = *this;

    // Grab block diagonal
    Field<tensor>& d = blockM.diag().asSquare();

    // Index offset for addressing the diagonal of the boundary faces
    const label start = ppatch.start();

    // We currently assume that the tangential gradient along fixedValue
    // boundaries is zero i.e. the values are uniform... do we?
    // This BC just forces the value so we don't make any simplifying
    // assumptions.
    forAll(ppatch, faceI)
    {
        // For displacement, the value is fixed so we will set the
        // diag and source so as to force this fixed value

        // The preconditioner will take care of scaling this coefficients;
        // however, scaling the coefficients here does affect the number of
        // iterations

        //const scalar scaleFac = (d[0].xx() + d[0].yy() + d[0].zz())/3.0;
        // if (mag(scaleFac) < SMALL)
        // {
        //     FatalErrorIn
        //     (
        //         "void blockFixedDisplacementFvPatchVectorField::"
        //         "insertBlockCoeffs\n"
        //         "(\n"
        //         "    const solidPolyMesh& solidMesh,\n"
        //         "    const surfaceScalarField& muf,\n"
        //         "    const surfaceScalarField& lambdaf,\n"
        //         "    const GeometricField<vector, fvPatchField, volMesh>&,\n"
        //         "    Field<vector>& blockB,\n"
        //         "    BlockLduMatrix<vector>& blockM\n"
        //         ") const"
        //     )   << "displacement scaleFac is zero" << abort(FatalError);
        // }

        const label varI = solidMesh.findOldVariableID(start + faceI);

        // Diagonal contribution for the boundary face
        //d[varI] += tensor(scaleFac*I);
        d[varI] += tensor(I);

        // Source contribution
        //blockB[varI] += scaleFac*pU[faceI];
        blockB[varI] += pU[faceI];
    }
}

void blockFixedDisplacementFvPatchVectorField::write(Ostream& os) const
{
    if (dispSeries_.size())
    {
        os.writeKeyword("displacementSeries") << nl;
        os << token::BEGIN_BLOCK << nl;
        dispSeries_.write(os);
        os << token::END_BLOCK << nl;
    }

    fixedValueFvPatchVectorField::write(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchVectorField,
    blockFixedDisplacementFvPatchVectorField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
