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

#include "blockSolidVelocityFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "fvcMeshPhi.H"
#include "fixedValuePointPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * * * //

void blockSolidVelocityFvPatchVectorField::makeInterp() const
{
    if (interpPtr_.valid())
    {
        FatalErrorIn
        (
            "void blockSolidVelocityFvPatchVectorField::makeInterp() const"
        ) << "pointer already set" << abort(FatalError);
    }

    interpPtr_.set(new primitivePatchInterpolation(patch().patch()));
}


primitivePatchInterpolation& blockSolidVelocityFvPatchVectorField::interp()
{
    if (interpPtr_.empty())
    {
        makeInterp();
    }

    return interpPtr_();
}


void blockSolidVelocityFvPatchVectorField::setPointDisplacement
(
    const vectorField& faceDisp
)
{
    const fvMesh& mesh = patch().boundaryMesh().mesh();

    if
    (
        mesh.foundObject<pointVectorField>
        (
            "point" + dimensionedInternalField().name()
        )
    )
    {
        const pointVectorField& pointD =
            mesh.lookupObject<pointVectorField>
            (
                "point" + dimensionedInternalField().name()
            );

        // Check if the boundary is fixedValue
        if
        (
            pointD.boundaryField()[patch().index()].type()
         == fixedValuePointPatchVectorField::typeName
        )
        {
            // Use const_cast to set boundary condition
            fixedValuePointPatchVectorField& patchPointD =
                refCast<fixedValuePointPatchVectorField>
                (
                    const_cast<pointVectorField&>
                    (
                        pointD
                    ).boundaryField()[patch().index()]
                );

            // Interpolate face values to the points
            patchPointD == interp().faceToPointInterpolate(faceDisp);
        }
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

blockSolidVelocityFvPatchVectorField::
blockSolidVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF),
    blockFvPatchVectorField(),
    velocity_(p.size(), vector::zero),
    velocitySeries_(),
    interpPtr_(NULL)
{}


blockSolidVelocityFvPatchVectorField::
blockSolidVelocityFvPatchVectorField
(
    const blockSolidVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper),
    blockFvPatchVectorField(),
    velocity_(ptf.velocity_),
    velocitySeries_(ptf.velocitySeries_),
    interpPtr_(NULL)
{}


blockSolidVelocityFvPatchVectorField::
blockSolidVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF),
    blockFvPatchVectorField(),
    velocity_(p.size(), vector::zero),
    velocitySeries_(),
    interpPtr_(NULL)
{
    Info<< "Creating " << type() << " boundary condition" << endl;

    // Read velocity
    if (dict.found("velocity"))
    {
        velocity_ = vectorField("velocity", dict, p.size());
    }
    else if (dict.found("velocitySeries"))
    {
        Info<< "    velocity is time-varying" << endl;
        velocitySeries_ =
            interpolationTable<vector>(dict.subDict("velocitySeries"));

        fvPatchField<vector>::operator==
        (
            velocitySeries_(this->db().time().timeOutputValue())
        );
    }
    else
    {
        FatalErrorIn(type() + "::solidVelocityFvPatchVectorField(...)")
            << "Either 'velocity' or 'velocitySeries' should be specified!"
            << abort(FatalError);
    }

    if (dict.found("value"))
    {
        Field<vector>::operator=(vectorField("value", dict, p.size()));
    }
    else
    {
        fvPatchVectorField::operator=(patchInternalField());
    }
}

blockSolidVelocityFvPatchVectorField::
blockSolidVelocityFvPatchVectorField
(
    const blockSolidVelocityFvPatchVectorField& ptf
)
:
    fixedValueFvPatchVectorField(ptf),
    blockFvPatchVectorField(),
    velocity_(ptf.velocity_),
    velocitySeries_(ptf.velocitySeries_),
    interpPtr_(NULL)
{}


blockSolidVelocityFvPatchVectorField::
blockSolidVelocityFvPatchVectorField
(
    const blockSolidVelocityFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(ptf, iF),
    blockFvPatchVectorField(),
    velocity_(ptf.velocity_),
    velocitySeries_(ptf.velocitySeries_),
    interpPtr_(NULL)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Map from self
void blockSolidVelocityFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedValueFvPatchVectorField::autoMap(m);

    velocity_.autoMap(m);
}


// Reverse-map the given fvPatchField onto this fvPatchField
void blockSolidVelocityFvPatchVectorField::rmap
(
    const fvPatchField<vector>& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchVectorField::rmap(ptf, addr);

    const blockSolidVelocityFvPatchVectorField& dmptf =
       refCast<const blockSolidVelocityFvPatchVectorField>(ptf);

    velocity_.rmap(dmptf.velocity_, addr);
}


void blockSolidVelocityFvPatchVectorField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    // Check if the velocity is time-varying
    if (velocitySeries_.size())
    {
        velocity_ = velocitySeries_(this->db().time().timeOutputValue());
    }

    vectorField disp = vectorField(patch().size(), vector::zero);

    if (dimensionedInternalField().name() == "DD")
    {
        // Incremental approach, so we wil set the increment of displacement for
        // this time-step
        disp = velocity_*db().time().deltaTValue();
    }
    else
    {
        // Lookup the old time total displacement
        const volVectorField& Dold =
            db().lookupObject<volVectorField>("D").oldTime();

        // The new total displacement is equal to Dold plus the increment of
        // displacement based on the current velocity and time-step
        disp =
            Dold.boundaryField()[patch().index()]
          + velocity_*db().time().deltaTValue();
    }

    // Set the displacement (or displacement increment) on the patch
    fvPatchField<vector>::operator==(disp);
    fixedValueFvPatchVectorField::updateCoeffs();

    // If the corresponding point displacement field has a fixedValue type
    // boundary condition, then we wil update it
    setPointDisplacement(disp);
}


Foam::tmp<Foam::Field<vector> >
blockSolidVelocityFvPatchVectorField::snGrad() const
{
    // Lookup previous boundary gradient
    const fvPatchField<tensor>& gradField =
        patch().lookupPatchField<volTensorField, tensor>
        (
            "grad(" + dimensionedInternalField().name() + ")"
        );

    // Unit normals
    const vectorField n = patch().nf();

    // Delta vectors
    const vectorField delta = patch().delta();

    // Correction vectors
    const vectorField k = delta - n*(n&delta);

    return
    (
        //*this - patchInternalField()
        //*this - (patchInternalField() + (k & gradField.patchInternalField()))
        //*this - (patchInternalField() + (k & gradField))
        (*this - (k & gradField)) - patchInternalField()
    )*this->patch().deltaCoeffs();
}


tmp<Field<vector> > blockSolidVelocityFvPatchVectorField::
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


void blockSolidVelocityFvPatchVectorField::insertBlockCoeffs
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
    const_cast<blockSolidVelocityFvPatchVectorField&>(*this).updateCoeffs();

    // Const reference to the patch field
    const vectorField& pU = *this;

    // Grab block diagonal
    Field<tensor>& d = blockM.diag().asSquare();

    // Index offset for addressing the diagonal of the boundary faces
    const label start = ppatch.start();

    // We will calculate the average of the current diagonal to scale the
    // coefficients for the fixedValue boundary conditions
    // If we don't do this then the convergence can be worse and also the
    // results can be strange
    // Note that this scale factor is just approximate and its exact value
    // does not matter
    const tensor averageDiag = gAverage(d);
    scalar diagSign =
        (averageDiag.xx() + averageDiag.yy() + averageDiag.zz());
    diagSign /= mag(diagSign);
    const scalar scaleFac = diagSign*(1.0/sqrt(3.0))*mag(averageDiag);

    if (mag(scaleFac) < SMALL)
    {
        FatalErrorIn
        (
            "void blockSolidVelocityFvPatchVectorField::insertBlockCoeffs\n"
            "(\n"
            "    const solidPolyMesh& solidMesh,\n"
            "    const surfaceScalarField& muf,\n"
            "    const surfaceScalarField& lambdaf,\n"
            "    const GeometricField<vector, fvPatchField, volMesh>& U,\n"
            "    Field<vector>& blockB,\n"
            "    BlockLduMatrix<vector>& blockM\n"
            ") const"
        )   << "The average diagonal coefficient is zero! The internal faces "
            << "should be discretised before inserting the boundary condition "
            << "equations" << abort(FatalError);
    }

    forAll(ppatch, faceI)
    {
        // Find the face index in the linear system
        const label varI = solidMesh.findOldVariableID(start + faceI);

        // For displacement, the value is fixed so we will set the
        // diag and source so as to force this fixed value

        // Diagonal contribution for the boundary face
        d[varI] += tensor(scaleFac*I);

        // Source contribution
        blockB[varI] += scaleFac*pU[faceI];
    }
}

void blockSolidVelocityFvPatchVectorField::write(Ostream& os) const
{
    if (velocitySeries_.size())
    {
        os.writeKeyword("velocitySeries") << nl;
        os << token::BEGIN_BLOCK << nl;
        velocitySeries_.write(os);
        os << token::END_BLOCK << nl;
    }
    else
    {
        velocity_.writeEntry("velocity", os);
    }

    fixedValueFvPatchVectorField::write(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchVectorField,
    blockSolidVelocityFvPatchVectorField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
