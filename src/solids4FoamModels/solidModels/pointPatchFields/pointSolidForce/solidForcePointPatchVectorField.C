/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
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

#include "solidForcePointPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "transformField.H"
#include "pointPatchFields.H"
#include "pointBoundaryMesh.H"
#include "pointMesh.H"
#ifdef OPENFOAMESIORFOUNDATION
    #include "Time.H"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

solidForcePointPatchVectorField::solidForcePointPatchVectorField
(
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF
)
:
    solidTractionPointPatchVectorField(p, iF),
    force_(p.size(), vector::zero),
    forceFieldPtr_(),
    curTimeIndex_(-1)
{}


solidForcePointPatchVectorField::solidForcePointPatchVectorField
(
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF,
    const dictionary& dict
)
:
    solidTractionPointPatchVectorField(p, iF),
    force_(p.size(), vector::zero),
    forceFieldPtr_(),
    curTimeIndex_(-1)
{
    // Check how force is defined
    if (dict.found("force") && dict.found("forceField"))
    {
        FatalErrorIn
        (
            "solidForcePointPatchVectorField::solidForcePointPatchVectorField"
        )   << "Only force or forceField can be "
            << "specified, not both!"
            << abort(FatalError);
    }
    else if (dict.found("forceField"))
    {
        Info<< "    force is specified as a field" << endl;

        // Lookup region name
        const word regionName
        (
            dict.lookupOrDefault<word>("region", fvMesh::defaultRegion)
        );

        // Lookup the fvMesh associated with the
        const fvMesh& fMesh =
            patch().boundaryMesh().mesh().time().lookupObject<fvMesh>
            (
                regionName
            );

        forceFieldPtr_.set
        (
            new volVectorField
            (
                IOobject
                (
                    word(dict.lookup("forceField")),
                    fMesh.time().timeName(),
                    fMesh,
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                fMesh
            )
        );
    }
    else
    {
        force_ = vectorField("force", dict, p.size());
    }
}


solidForcePointPatchVectorField::solidForcePointPatchVectorField
(
    const solidForcePointPatchVectorField& ptf,
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF,
    const PointPatchFieldMapper& mapper
)
:
    solidTractionPointPatchVectorField(p, iF),
#ifdef OPENFOAMFOUNDATION
    force_(mapper(ptf.force_)),
#else
    force_(ptf.force_, mapper),
#endif
    forceFieldPtr_(),
    curTimeIndex_(ptf.curTimeIndex_)
{}


#ifndef OPENFOAMFOUNDATION
solidForcePointPatchVectorField::solidForcePointPatchVectorField
(
    const solidForcePointPatchVectorField& ptf
)
:
    solidTractionPointPatchVectorField(ptf),
    force_(ptf.force_),
    forceFieldPtr_(),
    curTimeIndex_(ptf.curTimeIndex_)
{}
#endif


solidForcePointPatchVectorField::solidForcePointPatchVectorField
(
    const solidForcePointPatchVectorField& ptf,
    const DimensionedField<vector, pointMesh>& iF
)
:
    solidTractionPointPatchVectorField(ptf, iF),
    force_(ptf.force_),
    forceFieldPtr_(),
    curTimeIndex_(ptf.curTimeIndex_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Map and resize from self given a mapper
void solidForcePointPatchVectorField::autoMap
(
    const PointPatchFieldMapper& m
)
{
    //Field<vector>::autoMap(m);
#ifdef OPENFOAMFOUNDATION
    m(force_, force_);
#else
    force_.autoMap(m);
#endif
}


// Grab the values using rmap
void solidForcePointPatchVectorField::rmap
(
    const pointPatchField<vector>& ptf,
    const labelList& addr
)
{
    solidTractionPointPatchVectorField::rmap(ptf, addr);

    const solidForcePointPatchVectorField& tiptf =
      refCast<const solidForcePointPatchVectorField>(ptf);

    force_.rmap(tiptf.force_, addr);
}


void solidForcePointPatchVectorField::initEvaluate
(
    const Pstream::commsTypes commsType
)
{
    if (curTimeIndex_ != db().time().timeIndex())
    {
        curTimeIndex_ = db().time().timeIndex();

        // Called once per time-step
        if (forceFieldPtr_.valid())
        {
            // Force the traction field boundary conditions to update
            const_cast<volVectorField&>
            (
                forceFieldPtr_()
            ).correctBoundaryConditions();
        }
    }

    if (forceFieldPtr_.valid())
    {
        // Lookup the forces at the faces
        const vectorField faceForce
        (
            forceFieldPtr_().boundaryField()[patch().index()]
        );

        // Reset point forces to zero
        force_ = vector::zero;

        // Calculate the force associated with each point
        // Ideally we should use the dualMesh faces here
        // For simplicity, we will distribute the force of each face
        // between all the points in the face
        const polyPatch& ppatch =
            patch().boundaryMesh().mesh()().boundaryMesh()[patch().index()];
        const faceList& localFaces = ppatch.localFaces();
        forAll(localFaces, faceI)
        {
            const face& curFace = localFaces[faceI];
            const vector forceContribution = faceForce[faceI]/curFace.size();

            forAll(curFace, pI)
            {
                const label pointID = curFace[pI];
                force_[pointID] += forceContribution;
            }
        }
    }

    // Calculate the area associated with each point
    // Ideally we should use the dualMesh faces here
    // For simplicity, we will distribute the area of each face between all the
    // points in the face
    scalarField pointMagSf(force_.size(), 0.0);
    const polyPatch& ppatch =
        patch().boundaryMesh().mesh()().boundaryMesh()[patch().index()];
    const scalarField faceAreas(mag(ppatch.faceAreas()));
    const faceList& localFaces = ppatch.localFaces();
    forAll(localFaces, faceI)
    {
        const face& curFace = localFaces[faceI];
        const scalar areaContribution = faceAreas[faceI]/curFace.size();

        forAll(curFace, pI)
        {
            const label pointID = curFace[pI];
            pointMagSf[pointID] += areaContribution;
        }
    }

    // Set the traction at each point
    // Currently we assume a linear geometry approach
    traction() = force_/pointMagSf;

    solidTractionPointPatchVectorField::initEvaluate(commsType);
}


// Write
void solidForcePointPatchVectorField::write(Ostream& os) const
{
    solidTractionPointPatchVectorField::write(os);

    if (forceFieldPtr_.valid())
    {
        os.writeKeyword("forceField")
            << forceFieldPtr_().name() << token::END_STATEMENT << nl;
    }
    else
    {
#ifdef OPENFOAMFOUNDATION
        writeEntry(os, "force", force_);
#else
        force_.writeEntry("force", os);
#endif
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePointPatchTypeField
(
    pointPatchVectorField,
    solidForcePointPatchVectorField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
