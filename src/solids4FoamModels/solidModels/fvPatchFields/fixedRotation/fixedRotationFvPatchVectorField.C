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

#include "fixedRotationFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "fvcMeshPhi.H"
#include "RodriguesRotation.H"
#include "fixedValuePointPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

fixedRotationFvPatchVectorField::fixedRotationFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF),
    rotationAngle_(0.0),
    rotationAxis_(vector::zero),
    rotationOrigin_(vector::zero),
    origFaceCentres_(0, vector::zero),
    origPatchPoints_(0, vector::zero),
    angleSeries_(),
    dispSeries_(),
    originSeries_()
{}


fixedRotationFvPatchVectorField::fixedRotationFvPatchVectorField
(
    const fixedRotationFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper),
    rotationAngle_(ptf.rotationAngle_),
    rotationAxis_(ptf.rotationAxis_),
    rotationOrigin_(ptf.rotationOrigin_),
    origFaceCentres_(ptf.origFaceCentres_),
    origPatchPoints_(ptf.origPatchPoints_),
    angleSeries_(ptf.angleSeries_),
    dispSeries_(ptf.dispSeries_),
    originSeries_(ptf.originSeries_)
{}


fixedRotationFvPatchVectorField::fixedRotationFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF, dict),
    rotationAngle_(0.0),
    rotationAxis_(dict.lookup("rotationAxis")),
    rotationOrigin_(vector::zero),
    origFaceCentres_(patch().patch().faceCentres()),
    origPatchPoints_(patch().patch().localPoints()),
    angleSeries_(),
    dispSeries_(),
    originSeries_()
{
    // Check if angle is time-varying
    if (dict.found("rotationAngleSeries"))
    {
        Info<< "    angle is time-varying" << endl;
        angleSeries_ =
            interpolationTable<scalar>(dict.subDict("rotationAngleSeries"));
    }
    else
    {
        rotationAngle_ = readScalar(dict.lookup("rotationAngle"));
    }

    // Check if there is a time-varying translation
    if (dict.found("displacementSeries"))
    {
        Info<< "    translation is time-varying" << endl;
        dispSeries_ =
            interpolationTable<vector>(dict.subDict("displacementSeries"));
    }

    // Check if origin is time-varying
    if (dict.found("rotationOriginSeries"))
    {
        Info<< "    origin is time-varying" << endl;
        originSeries_ =
            interpolationTable<vector>(dict.subDict("rotationOriginSeries"));
    }
    else
    {
        rotationOrigin_ = vector(dict.lookup("rotationOrigin"));
    }

    if (dict.found("value"))
    {
        fvPatchField<vector>::operator==
        (
            vectorField("value", dict, p.size())
        );
    }
    else
    {
        fvPatchField<vector>::operator==
        (
            vectorField(p.size(), vector::zero)
        );
    }

    if
    (
#ifdef OPENFOAMESIORFOUNDATION
        internalField().name() != "D"
     && internalField().name() != "DD"
#else
        dimensionedInternalField().name() != "D"
     && dimensionedInternalField().name() != "DD"
#endif
    )
    {
        FatalErrorIn
        (
            "fixedRotationFvPatchVectorField::"
            "fixedRotationFvPatchVectorField(...)"
        )   << "The displacement field should be D or DD"
            << abort(FatalError);
    }
}

#ifndef OPENFOAMFOUNDATION
fixedRotationFvPatchVectorField::fixedRotationFvPatchVectorField
(
    const fixedRotationFvPatchVectorField& pivpvf
)
:
    fixedValueFvPatchVectorField(pivpvf),
    rotationAngle_(pivpvf.rotationAngle_),
    rotationAxis_(pivpvf.rotationAxis_),
    rotationOrigin_(pivpvf.rotationOrigin_),
    origFaceCentres_(pivpvf.origFaceCentres_),
    origPatchPoints_(pivpvf.origPatchPoints_),
    angleSeries_(pivpvf.angleSeries_),
    dispSeries_(pivpvf.dispSeries_),
    originSeries_(pivpvf.originSeries_)
{}
#endif

fixedRotationFvPatchVectorField::fixedRotationFvPatchVectorField
(
    const fixedRotationFvPatchVectorField& pivpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(pivpvf, iF),
    rotationAngle_(pivpvf.rotationAngle_),
    rotationAxis_(pivpvf.rotationAxis_),
    rotationOrigin_(pivpvf.rotationOrigin_),
    origFaceCentres_(pivpvf.origFaceCentres_),
    origPatchPoints_(pivpvf.origPatchPoints_),
    angleSeries_(pivpvf.angleSeries_),
    dispSeries_(pivpvf.dispSeries_),
    originSeries_(pivpvf.originSeries_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::Field<vector> > fixedRotationFvPatchVectorField::
snGrad() const
{
    // fixedValue snGrad with no correction
    // return (*this - patchInternalField())*this->patch().deltaCoeffs();

    const fvPatchField<tensor>& gradField =
        patch().lookupPatchField<volTensorField, tensor>
        (
#ifdef OPENFOAMESIORFOUNDATION
            "grad(" + internalField().name() + ")"
#else
            "grad(" + dimensionedInternalField().name() + ")"
#endif
        );

    // Face unit normals
    const vectorField n(patch().nf());

    // Delta vectors
    const vectorField delta(patch().delta());

    // Non-orthogonal correction vectors
    const vectorField k((I - sqr(n)) & delta);

    return
    (
        *this
      - (patchInternalField() + (k & gradField.patchInternalField()))
    )*patch().deltaCoeffs();
}


void fixedRotationFvPatchVectorField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    if (angleSeries_.size())
    {
        rotationAngle_ = angleSeries_(this->db().time().timeOutputValue());
    }

    if (originSeries_.size())
    {
        rotationOrigin_ = originSeries_(this->db().time().timeOutputValue());
    }

    // Rotation tensor
    const tensor rotMat = RodriguesRotation(rotationAxis_, rotationAngle_);

    vectorField newFaceCentres
    (
        (rotMat & (origFaceCentres_ - rotationOrigin_)) + rotationOrigin_
    );

    vectorField disp(newFaceCentres - origFaceCentres_);

    // Superimposed translation
    if (dispSeries_.size())
    {
        disp += dispSeries_(this->db().time().timeOutputValue());
    }

    const fvMesh& mesh = patch().boundaryMesh().mesh();

#ifdef OPENFOAMESIORFOUNDATION
    if (internalField().name() == "DD")
#else
    if (dimensionedInternalField().name() == "DD")
#endif
    {
        const volVectorField& D = mesh.lookupObject<volVectorField>("D");

        disp -= D.oldTime().boundaryField()[patch().index()];
    }

    // Set the face displacement
    fvPatchField<vector>::operator==
    (
        disp
    );

    // If the point displacement field is found and is fixedValue then we will
    // update it using const_cast
    if
    (
        mesh.foundObject<pointVectorField>
        (
#ifdef OPENFOAMESIORFOUNDATION
            "point" + internalField().name()
#else
            "point" + dimensionedInternalField().name()
#endif
        )
    )
    {
        pointVectorField& pointDField =
            const_cast<pointVectorField&>
            (
                mesh.lookupObject<pointVectorField>
                (
#ifdef OPENFOAMESIORFOUNDATION
                    "point" + internalField().name()
#else
                    "point" + dimensionedInternalField().name()
#endif
                )
            );

        if
        (
            pointDField.boundaryField()[patch().index()].type()
         == fixedValuePointPatchVectorField::typeName
        )
        {
            fixedValuePointPatchVectorField& pointD =
                refCast<fixedValuePointPatchVectorField>
                (
#ifdef OPENFOAMESIORFOUNDATION
                    pointDField.boundaryFieldRef()[patch().index()]
#else
                    pointDField.boundaryField()[patch().index()]
#endif
                );

            const vectorField newPatchPoints
            (
                (rotMat & (origPatchPoints_ - rotationOrigin_))
              + rotationOrigin_
            );

            vectorField pointDisp(newPatchPoints - origPatchPoints_);

            // Superimposed translation
            if (dispSeries_.size())
            {
                pointDisp += dispSeries_(this->db().time().timeOutputValue());
            }

            const labelList& meshPoints =
                mesh.boundaryMesh()[patch().index()].meshPoints();

#ifdef OPENFOAMESIORFOUNDATION
            if (internalField().name() == "DD")
#else
            if (dimensionedInternalField().name() == "DD")
#endif
            {
                // Lookup the accumulated total displacement
                const pointVectorField& pointDFieldOld =
                    mesh.lookupObject<pointVectorField>("pointD").oldTime();

                forAll(meshPoints, pointI)
                {
                    pointDisp[pointI] -= pointDFieldOld[meshPoints[pointI]];
                }
            }

            pointD == pointDisp;
        }
    }

    fixedValueFvPatchVectorField::updateCoeffs();
}


tmp<Field<vector> > fixedRotationFvPatchVectorField::
gradientBoundaryCoeffs() const
{
    const fvPatchField<tensor>& gradField =
        patch().lookupPatchField<volTensorField, tensor>
        (
#ifdef OPENFOAMESIORFOUNDATION
            "grad(" + internalField().name() + ")"
#else
            "grad(" + dimensionedInternalField().name() + ")"
#endif
        );

    // Face unit normals
    const vectorField n(patch().nf());

    // Delta vectors
    const vectorField delta(patch().delta());

    // Non-orthogonal correction vectors
    const vectorField k((I - sqr(n)) & delta);

    return patch().deltaCoeffs()*
    (
       *this - (k & gradField.patchInternalField())
    );
}


void fixedRotationFvPatchVectorField::write(Ostream& os) const
{
    fixedValueFvPatchVectorField::write(os);

    if (angleSeries_.size())
    {
        os.writeKeyword("rotationAngleSeries") << nl;
        os << token::BEGIN_BLOCK << nl;
        angleSeries_.write(os);
        os << token::END_BLOCK << nl;
    }
    else
    {
        os.writeKeyword("rotationAngle")
            << rotationAngle_
            << token::END_STATEMENT << nl;
    }
    if (dispSeries_.size())
    {
        os.writeKeyword("displacementSeries") << nl;
        os << token::BEGIN_BLOCK << nl;
        dispSeries_.write(os);
        os << token::END_BLOCK << nl;
    }
    if (originSeries_.size())
    {
        os.writeKeyword("rotationOriginSeries") << nl;
        os << token::BEGIN_BLOCK << nl;
        originSeries_.write(os);
        os << token::END_BLOCK << nl;
    }
    else
    {
        os.writeKeyword("rotationOrigin")
            << rotationOrigin_
            << token::END_STATEMENT << nl;
    }
    os.writeKeyword("rotationAxis")
        << rotationAxis_
        << token::END_STATEMENT << nl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchVectorField,
    fixedRotationFvPatchVectorField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
