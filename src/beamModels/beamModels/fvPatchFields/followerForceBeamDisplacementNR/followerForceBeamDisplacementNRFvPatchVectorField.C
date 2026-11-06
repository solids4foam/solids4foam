/*---------------------------------------------------------------------------*\
License
    This file is part of beamFoam.

    beamFoam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    beamFoam is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with beamFoam.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "followerForceBeamDisplacementNRFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

followerForceBeamDisplacementNRFvPatchVectorField::
followerForceBeamDisplacementNRFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    forceBeamDisplacementNRFvPatchVectorField(p, iF),
    followerForce_(p.size(), vector::zero),
    //followerForceDir_(p.size(), vector::zero),
    followerForceSeries_()
{
    // fixedFvPatchVectorField::operator=(patchInternalField());
    // gradient() = vector::zero;
}


followerForceBeamDisplacementNRFvPatchVectorField::
followerForceBeamDisplacementNRFvPatchVectorField
(
    const followerForceBeamDisplacementNRFvPatchVectorField& tdpvf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    forceBeamDisplacementNRFvPatchVectorField(tdpvf, p, iF, mapper),
    followerForce_(tdpvf.followerForce_, mapper),
   // followerForceDir_(tdpvf.followerForceDir_, mapper),
    followerForceSeries_(tdpvf.followerForceSeries_)
{}


followerForceBeamDisplacementNRFvPatchVectorField::
followerForceBeamDisplacementNRFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    forceBeamDisplacementNRFvPatchVectorField(p, iF),
    // forceBeamDisplacementNRFvPatchVectorField(p, iF, dict),
    followerForce_(p.size(), vector::zero),
   // followerForceDir_(p.size(), vector::zero),
    followerForceSeries_()
{
    if (dict.found("value"))
    {
        fixedValueFvPatchVectorField::operator==
        (
            vectorField("value", dict, p.size())
        );
    }
    else
    {
        fixedValueFvPatchVectorField::operator==
        (
            patchInternalField()
        );
    }

    if (dict.found("followerForceSeries"))
    {
        Info<< "follower force is time-varying" << endl;
	followerForceSeries_ =
            interpolationTable<vector>(dict.subDict("followerForceSeries"));

        followerForce_ =
            followerForceSeries_(this->db().time().timeOutputValue());
    }
    else
    {
        FatalErrorIn
        (
            "followerForceBeamDisplacementNRFvPatchVectorField::"
            "followerForceBeamDisplacementNRFvPatchVectorField"
            "("
            "    const fvPatch& p,"
            "    const DimensionedField<vector, volMesh>& iF,"
            "    const dictionary& dict"
            ")"
        )
            << "followerForceSeries is not defined"
            << abort(FatalError);
    }
}


followerForceBeamDisplacementNRFvPatchVectorField::
followerForceBeamDisplacementNRFvPatchVectorField
(
    const followerForceBeamDisplacementNRFvPatchVectorField& tdpvf
)
:
    forceBeamDisplacementNRFvPatchVectorField(tdpvf),
    followerForce_(tdpvf.followerForce_),
   // followerForceDir_(tdpvf.followerForceDir_),
    followerForceSeries_(tdpvf.followerForceSeries_)
{}


followerForceBeamDisplacementNRFvPatchVectorField::
followerForceBeamDisplacementNRFvPatchVectorField
(
    const followerForceBeamDisplacementNRFvPatchVectorField& tdpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    forceBeamDisplacementNRFvPatchVectorField(tdpvf, iF),
    followerForce_(tdpvf.followerForce_),
  //  followerForceDir_(tdpvf.followerForceDir_),
    followerForceSeries_(tdpvf.followerForceSeries_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void followerForceBeamDisplacementNRFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    forceBeamDisplacementNRFvPatchVectorField::autoMap(m);
    followerForce_.autoMap(m);
}


// Reverse-map the given fvPatchField onto this fvPatchField
void followerForceBeamDisplacementNRFvPatchVectorField::rmap
(
    const fvPatchVectorField& ptf,
    const labelList& addr
)
{
    forceBeamDisplacementNRFvPatchVectorField::rmap(ptf, addr);

    const followerForceBeamDisplacementNRFvPatchVectorField& dmptf =
        refCast<const followerForceBeamDisplacementNRFvPatchVectorField>(ptf);

    followerForce_.rmap(dmptf.followerForce_, addr);
}


// Update the coefficients associated with the patch field
void followerForceBeamDisplacementNRFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

        followerForce_ = followerForceSeries_
        (
            this->db().time().timeOutputValue()
        );

	/*
        const vectorField& dR0Ds =
            patch().lookupPatchField<surfaceVectorField, vector>("dR0Ds");

        const vectorField dWdS =
            db().lookupObject<volVectorField>("W")
           .boundaryField()[this->patch().index()].snGrad();

        // const vectorField dWdS =
        //     db().lookupObject<volVectorField>("W").oldTime()
        //    .boundaryField()[this->patch().index()].snGrad();

        vectorField dRdS = dR0Ds + dWdS;
        vectorField t = dRdS/(mag(dRdS) + SMALL);
        */

	const tensorField& Lambdaf =

            patch().lookupPatchField<surfaceTensorField, tensor>("Lambdaf");


	// Extraction of column vectors from the rotation matrix
	// to construct the body-attached coordinate
	// frame t1,t2,t3

	// 1st column
	const scalarField txx(Lambdaf.component(tensor::XX));
	const scalarField tyx(Lambdaf.component(tensor::YX));
	const scalarField tzx(Lambdaf.component(tensor::ZX));

	// 2nd column
	const scalarField txy(Lambdaf.component(tensor::XY));
	const scalarField tyy(Lambdaf.component(tensor::YY));
	const scalarField tzy(Lambdaf.component(tensor::ZY));

	// 3rd column
	const scalarField txz(Lambdaf.component(tensor::XZ));
	const scalarField tyz(Lambdaf.component(tensor::YZ));
	const scalarField tzz(Lambdaf.component(tensor::ZZ));


	// Fixed reference coordinate system
	const vector i(1, 0, 0);
	const vector j(0, 1, 0);
	const vector k(0, 0, 1);

	// Body-attached coordinate frame
	vectorField t1 ((i*txx + j*tyx + k*tzx));
	vectorField t2 ((i*txy + j*tyy + k*tzy));
	vectorField t3 ((i*txz + j*tyz + k*tzz));

	// t3 = t3/mag(t3);
	// followerForceDir_ = t3;

        // Current force vector
	vectorField curForce
	(
	      followerForce_.component(0)*t1
	    + followerForce_.component(1)*t2
	    + followerForce_.component(2)*t3
	);
	//vectorField followerForceDir_ = curForce/mag(curForce);

        //vectorField curForce =
        //    followerForce_.component(0)*followerForceDir_;

        force() = curForce;


    // fixedValueFvPatchVectorField::updateCoeffs();
}


// Write
void followerForceBeamDisplacementNRFvPatchVectorField::write(Ostream& os) const
{
    forceBeamDisplacementNRFvPatchVectorField::write(os);

    if (followerForceSeries_.size())
    {
        os.writeKeyword("followerForceSeries") << nl;
        os << token::BEGIN_BLOCK << nl;
        followerForceSeries_.write(os);
        os << token::END_BLOCK << nl;
    }

    followerForce_.writeEntry("followerForce", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchVectorField,
    followerForceBeamDisplacementNRFvPatchVectorField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
