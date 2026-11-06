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

#include "almDragContribution.H"
#include "addToRunTimeSelectionTable.H"
#include "HermiteSpline.H"
#include "samplingFluid.H"
#include "groundContactContribution.H"
#include "beamParallelChecks.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(almDragContribution, 0);
addToRunTimeSelectionTable
(
    beamMomentumContribution, almDragContribution, dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

almDragContribution::almDragContribution
(
    const word& name,
    const dictionary& dict
)
:
    beamMomentumContribution(name, dict),
    coeffs_(dict.subDict(name + "Coeffs")),
    Cdn_(readScalar(coeffs_.lookup("Cdn"))),
    Cdt_(readScalar(coeffs_.lookup("Cdt"))),
    forceRelaxation_
    (
        coeffs_.lookupOrDefault<scalar>("forceRelaxation", 1.0)
    ),
    samplingRadius_
    (
        coeffs_.lookupOrDefault<scalar>("samplingRadius", 5.0)
    ),
    samplingReferenceLength_
    (
        coeffs_.lookupOrDefault<scalar>("almSamplingReferenceLength", -1.0)
    ),
    dragReferenceLength_
    (
        coeffs_.lookupOrDefault<scalar>("almDragReferenceLength", -1.0)
    ),
    samplingPlaneNormal_(coeffs_.lookup("almSamplingPlaneNormal")),
    upstreamDir_(coeffs_.lookup("almUpstreamDir")),
    groundContactActive_
    (
        coeffs_.lookupOrDefault<bool>("groundContactActive", false)
    ),
    groundPlane_
    (
        coeffs_.found("groundZ")
     || coeffs_.found("groundPoint")
     || coeffs_.found("groundNormal")
      ? groundContactContribution::readGroundPlane(coeffs_)
      : plane(point::zero, vector(0, 0, 1))
    ),
    searchEnginePtr_(),
    decompositionChecked_(false)
{
    Info<< "Found beamMomentumContribution type: " << typeName << endl;

    if (mag(samplingPlaneNormal_) <= SMALL)
    {
        FatalIOErrorInFunction(coeffs_)
            << "almSamplingPlaneNormal magnitude must be positive. "
            << "Current value: " << samplingPlaneNormal_
            << exit(FatalIOError);
    }

    if (mag(upstreamDir_) <= SMALL)
    {
        FatalIOErrorInFunction(coeffs_)
            << "almUpstreamDir magnitude must be positive. "
            << "Current value: " << upstreamDir_
            << exit(FatalIOError);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void almDragContribution::preEvolve(const beamModel& bm)
{
    const fvMesh& beamMesh = bm.solutionW().mesh();

    // The sampling below builds one spline through the whole beam on the
    // master and broadcasts it, so the beam region must not be split across
    // processors. Collective, hence outside any master guard.
    if (!decompositionChecked_)
    {
        checkBeamOnMasterProcessor
        (
            beamMesh,
            "The " + typeName + " beam momentum contribution"
        );

        decompositionChecked_ = true;
    }

    const fvMesh& fluidMesh =
        beamMesh.time().db().parent().lookupObject<fvMesh>("region0");

    if (!searchEnginePtr_.valid())
    {
        searchEnginePtr_.reset(new meshSearch(fluidMesh));
    }

    // The "almForce" field is owned by the beam model and registered on the
    // beam mesh.  Look it up by name to avoid coupling to a concrete derived
    // beam type.
    volVectorField& almForce =
        beamMesh.lookupObjectRef<volVectorField>("almForce");

    // Per-cell line force at the start of this evolve() - used for
    // outer-corrector under-relaxation
    const vectorField almForceStart(almForce.internalField());

    // Beam-aligned tangent at beam cell centres. Built on the master, which
    // holds every beam cell, and left empty on the other ranks so that it
    // matches the local beam mesh size; getFluidVelocity broadcasts it for
    // the collective sampling.
    vectorField dRdScell;

    if (Pstream::master())
    {
        // The 0 argument selects the 0th beam (beam index bI = 0).
        HermiteSpline spline
        (
            bm.currentBeamPoints(0),
            bm.currentBeamTangents(0)
        );

        dRdScell = spline.midPointDerivatives();
    }

    const vectorField dRdScellHat(dRdScell/(mag(dRdScell) + SMALL));

    // Beam cell-centre global coordinates (reference + total displacement)
    const volVectorField& W = bm.solutionW();
    const vectorField beamCellCenterCoord =
        beamMesh.C() + W.mesh().lookupObject<volVectorField>("refW") + W;

    labelList seedCellIDs(beamMesh.nCells(), -1);

    const scalar samplingRefL =
        (samplingReferenceLength_ > 0) ? samplingReferenceLength_ : bm.R();

    if (samplingRefL <= SMALL)
    {
        FatalErrorInFunction
            << "almSamplingReferenceLength must be positive. "
            << "Current value: " << samplingRefL
            << abort(FatalError);
    }

    const scalar dragRefD =
        (dragReferenceLength_ > 0) ? dragReferenceLength_ : 2.0*bm.R();

    if (dragRefD <= SMALL)
    {
        FatalErrorInFunction
            << "almDragReferenceLength must be positive. "
            << "Current value: " << dragRefD
            << abort(FatalError);
    }

    std::tuple
    <
        tmp<volVectorField>,
        tmp<volScalarField>,
        tmp<volVectorField>,
        labelList
    > fluidInfo =
        getFluidVelocity
        (
            fluidMesh,
            beamMesh,
            beamCellCenterCoord,
            seedCellIDs,
            groundPlane_,
            groundContactActive_,
            dRdScell,
            samplingRefL,
            samplingRadius_,
            samplingPlaneNormal_,
            upstreamDir_,
            searchEnginePtr_()
        );

    const volVectorField fluidU(std::get<0>(fluidInfo).ref());
    const volScalarField cellMarker(std::get<1>(fluidInfo).ref());

    if (beamMesh.time().writeTime())
    {
        cellMarker.write();
    }

    // Relative velocity (fluid minus beam)
    const volVectorField& U =
        beamMesh.lookupObject<volVectorField>("U");
    const vectorField Urel(fluidU.internalField() - U.internalField());

    const vectorField Ut
    (
        (Urel & dRdScellHat)*dRdScellHat
    );

    const vectorField Un
    (
        Urel - Ut
    );

    const scalar rhoF = bm.rhoFluid().value();
    const scalar coeffN = 0.5*rhoF*Cdn_*dragRefD;
    const scalar coeffT = 0.5*rhoF*Cdt_*dragRefD;

    vectorField& almForceI = almForce.primitiveFieldRef();

    forAll(almForceI, cellI)
    {
        const vector un = Un[cellI];
        const vector ut = Ut[cellI];
        const vector vs = un + ut;
        const scalar vsMag = mag(vs);

        const vector fnLineTarget = -coeffN*vsMag*un;
        const vector ftLineTarget = -coeffT*vsMag*ut;
        const vector FLineTarget = fnLineTarget + ftLineTarget;

        almForceI[cellI] =
            forceRelaxation_*FLineTarget
          + (1.0 - forceRelaxation_)*almForceStart[cellI];
    }

    almForce.correctBoundaryConditions();

    const vector Fsum(sum(almForce*bm.L()).value());
    Info<< "sum F on the beam (relaxed) = "
        << "(" << Fsum.x() << ", "
        << Fsum.y() << ", "
        << Fsum.z() << ") N" << nl;
}


tmp<Field<scalarSquareMatrix>> almDragContribution::diagCoeff
(
    const beamModel& bm,
    const volVectorField& U,
    const volVectorField& Accl
)
{
    const fvMesh& mesh = bm.solutionW().mesh();

    return tmp<Field<scalarSquareMatrix>>
    (
        new Field<scalarSquareMatrix>
        (
            mesh.nCells(), scalarSquareMatrix(6, 0.0)
        )
    );
}


tmp<vectorField> almDragContribution::linearMomentumSource
(
    const beamModel& bm,
    const volVectorField& U,
    const volVectorField& Accl
)
{
    const fvMesh& beamMesh = bm.solutionW().mesh();

    const volVectorField& almForce =
        beamMesh.lookupObject<volVectorField>("almForce");

    const volScalarField& L = bm.L();

    tmp<vectorField> tresult(new vectorField(beamMesh.nCells(), vector::zero));
    vectorField& result = tresult.ref();

    const vectorField& almForceI = almForce.internalField();
    const scalarField& LI = L.internalField();

    forAll(result, cellI)
    {
        result[cellI] = almForceI[cellI]*LI[cellI];
    }

    return tresult;
}


tmp<vectorField> almDragContribution::angularMomentumSource
(
    const beamModel& bm,
    const volVectorField& U,
    const volVectorField& Accl
)
{
    const fvMesh& mesh = bm.solutionW().mesh();

    return tmp<vectorField>
    (
        new vectorField(mesh.nCells(), vector::zero)
    );
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
