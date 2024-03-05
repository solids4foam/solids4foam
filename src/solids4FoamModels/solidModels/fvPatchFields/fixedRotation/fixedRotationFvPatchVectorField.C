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

#include "fixedRotationFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "fvcMeshPhi.H"
#include "RodriguesRotation.H"
#include "fixedValuePointPatchFields.H"
#include "patchCorrectionVectors.H"

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
    nonOrthogonalCorrections_(true),
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
    const fixedRotationFvPatchVectorField& pvf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(pvf, p, iF, mapper),
    nonOrthogonalCorrections_(pvf.nonOrthogonalCorrections_),
    rotationAngle_(pvf.rotationAngle_),
    rotationAxis_(pvf.rotationAxis_),
    rotationOrigin_(pvf.rotationOrigin_),
    origFaceCentres_(pvf.origFaceCentres_),
    origPatchPoints_(pvf.origPatchPoints_),
    angleSeries_(pvf.angleSeries_),
    dispSeries_(pvf.dispSeries_),
    originSeries_(pvf.originSeries_)
{}


fixedRotationFvPatchVectorField::fixedRotationFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF, dict),
    nonOrthogonalCorrections_
    (
        dict.lookupOrDefault<Switch>("nonOrthogonalCorrections", true)
    ),
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
#ifdef OPENFOAM_NOT_EXTEND
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

#ifndef OPENFOAM_ORG
fixedRotationFvPatchVectorField::fixedRotationFvPatchVectorField
(
    const fixedRotationFvPatchVectorField& pvf
)
:
    fixedValueFvPatchVectorField(pvf),
    nonOrthogonalCorrections_(pvf.nonOrthogonalCorrections_),
    rotationAngle_(pvf.rotationAngle_),
    rotationAxis_(pvf.rotationAxis_),
    rotationOrigin_(pvf.rotationOrigin_),
    origFaceCentres_(pvf.origFaceCentres_),
    origPatchPoints_(pvf.origPatchPoints_),
    angleSeries_(pvf.angleSeries_),
    dispSeries_(pvf.dispSeries_),
    originSeries_(pvf.originSeries_)
{}
#endif

fixedRotationFvPatchVectorField::fixedRotationFvPatchVectorField
(
    const fixedRotationFvPatchVectorField& pvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(pvf, iF),
    nonOrthogonalCorrections_(pvf.nonOrthogonalCorrections_),
    rotationAngle_(pvf.rotationAngle_),
    rotationAxis_(pvf.rotationAxis_),
    rotationOrigin_(pvf.rotationOrigin_),
    origFaceCentres_(pvf.origFaceCentres_),
    origPatchPoints_(pvf.origPatchPoints_),
    angleSeries_(pvf.angleSeries_),
    dispSeries_(pvf.dispSeries_),
    originSeries_(pvf.originSeries_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::Field<vector> > fixedRotationFvPatchVectorField::
snGrad() const
{
    if (nonOrthogonalCorrections_)
    {
        const fvPatchField<tensor>& gradField =
            patch().lookupPatchField<volTensorField, tensor>
            (
            #ifdef OPENFOAM_NOT_EXTEND
                "grad(" + internalField().name() + ")"
            #else
                "grad(" + dimensionedInternalField().name() + ")"
            #endif
            );

        // Non-orthogonal correction vectors
        const vectorField k(patchCorrectionVectors(patch()));

        return
        (
            *this
          - (patchInternalField() + (k & gradField.patchInternalField()))
        )*patch().deltaCoeffs();
    }
    else
    {
        // fixedValue snGrad with no correction
        return (*this - patchInternalField())*patch().deltaCoeffs();
    }
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

#ifdef OPENFOAM_NOT_EXTEND
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
#ifdef OPENFOAM_NOT_EXTEND
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
#ifdef OPENFOAM_NOT_EXTEND
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
#ifdef OPENFOAM_NOT_EXTEND
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

#ifdef OPENFOAM_NOT_EXTEND
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
    if (nonOrthogonalCorrections_)
    {
        const fvPatchField<tensor>& gradField =
            patch().lookupPatchField<volTensorField, tensor>
            (
            #ifdef OPENFOAM_NOT_EXTEND
                "grad(" + internalField().name() + ")"
            #else
                "grad(" + dimensionedInternalField().name() + ")"
            #endif
            );

        // Non-orthogonal correction vectors
        const vectorField k(patchCorrectionVectors(patch()));

        return
        (
            patch().deltaCoeffs()
           *(*this - (k & gradField.patchInternalField()))
        );
    }
    else
    {
        return (patch().deltaCoeffs()*(*this));
    }
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
