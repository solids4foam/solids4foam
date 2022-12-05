/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
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

#include "interFluid.H"
#include "volFields.H"
#include "fvm.H"
#include "fvc.H"
#include "fvMatrices.H"
#include "addToRunTimeSelectionTable.H"
#include "findRefCell.H"
#include "adjustPhi.H"
#include "zeroGradientFvPatchFields.H"
#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace fluidModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(interFluid, 0);
addToRunTimeSelectionTable(fluidModel, interFluid, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

interFluid::interFluid
(
    Time& runTime,
    const word& region
)
:
    fluidModel(typeName, runTime, region),
    pd_
    (
        IOobject
        (
            "pd",
            runTime.timeName(),
            mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh()
    ),
    alpha1_
    (
        IOobject
        (
            "alpha1",
            runTime.timeName(),
            mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh()
    ),
    twoPhaseProperties_(U(), phi(), "alpha1"),
    rho1_(twoPhaseProperties_.rho1()),
    rho2_(twoPhaseProperties_.rho2()),
    rho_
    (
        IOobject
        (
            "rho",
            runTime.timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT
        ),
        alpha1_*rho1_ + (scalar(1) - alpha1_)*rho2_,
        alpha1_.boundaryField().types()
    ),
    rhoPhi_
    (
        IOobject
        (
            "rho*phi",
            runTime.timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        rho1_*phi()
    ),
    interface_(alpha1_, U(), twoPhaseProperties_),
    pdRefCell_(0),
    pdRefValue_(0.0),
    pRefValue_(0.0),
    correctPhi_(pimple().dict().lookupOrDefault("correctPhi", false)),
    checkMeshCourantNo_
    (
        pimple().dict().lookupOrDefault("checkMeshCourantNo", false)
    ),
    turbulence_
    (
        incompressible::turbulenceModel::New(U(), phi(), twoPhaseProperties_)
    ),
    pcorrTypes_
    (
        pd_.boundaryField().size(),
        zeroGradientFvPatchScalarField::typeName
    ),
    sumLocalContErr_(0),
    globalContErr_(0),
    cumulativeContErr_(0)
{
    UisRequired();

    // Set pcorrTypes
    for (label i = 0; i < pd_.boundaryField().size(); i++)
    {
        if (pd_.boundaryField()[i].fixesValue())
        {
            pcorrTypes_[i] = fixedValueFvPatchScalarField::typeName;
        }
    }

    // Reset p dimensions: we should allow p not to be read!
    Info<< "Resetting the dimensions of p" << endl;
    p().dimensions().reset(dimPressure);
    p() = pd_ + rho_*(g() & mesh().C());

    // Store old-time rho
    rho_.oldTime();

    setRefCell(pd_, pimple().dict(), pdRefCell_, pdRefValue_);
    mesh().schemesDict().setFluxRequired(pd_.name());

    if (pd_.needReference())
    {
        pRefValue_ = readScalar(pimple().dict().lookup("pRefValue"));

        p() += dimensionedScalar
        (
            "p",
            p().dimensions(),
            pRefValue_ - getRefCellValue(p(), pdRefCell_)
        );
    }
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<vectorField> interFluid::patchViscousForce(const label patchID) const
{
    tmp<vectorField> tvF
    (
        new vectorField(mesh().boundary()[patchID].size(), vector::zero)
    );

    tvF() =
        (
            mesh().boundary()[patchID].nf()
          & (-turbulence_->devReff()().boundaryField()[patchID])
        );

    return tvF;
}


tmp<scalarField> interFluid::patchPressureForce(const label patchID) const
{
    tmp<scalarField> tpF
    (
        new scalarField(mesh().boundary()[patchID].size(), 0)
    );

    tpF() = p().boundaryField()[patchID];

    return tpF;
}


bool interFluid::evolve()
{
    Info<< "Evolving fluid model: " << this->type() << endl;

    // Take references
    const Time& runTime = fluidModel::runTime();
    dynamicFvMesh& mesh = fluidModel::mesh();
    pimpleControl& pimple = fluidModel::pimple();
    twoPhaseMixture& twoPhaseProperties = twoPhaseProperties_;
    interfaceProperties& interface = interface_;
    autoPtr<incompressible::turbulenceModel>& turbulence = turbulence_;
    volVectorField& U = this->U();
    volScalarField& p = this->p();
    volScalarField& pd = pd_;
    volScalarField& alpha1 = alpha1_;
    surfaceScalarField& phi = this->phi();
    surfaceScalarField& rhoPhi = rhoPhi_;
    volScalarField& rho = rho_;
    const dimensionedScalar& rho1 = rho1_;
    const dimensionedScalar& rho2 = rho2_;
    const label pdRefCell = pdRefCell_;
    const scalar& pdRefValue = pdRefValue_;
    const scalar& pRefValue = pRefValue_;
    scalar& sumLocalContErr = sumLocalContErr_;
    scalar& globalContErr = globalContErr_;
    scalar& cumulativeContErr = cumulativeContErr_;
    const wordList& pcorrTypes = pcorrTypes_;

    // For now, we check for FSI mesh update here; a better way will be to
    // create a FSI dynamicFvMesh: to-do
    bool meshChanged = false;
    if (fluidModel::fsiMeshUpdate())
    {
        // The FSI interface is in charge of calling mesh.update()
        meshChanged = fluidModel::fsiMeshUpdateChanged();
    }
    else
    {
        meshChanged = mesh.update();
        reduce(meshChanged, orOp<bool>());
    }

    #include "volContinuity.H"

    // Update gh fields as the mesh may have moved
    volScalarField gh("gh", g() & mesh.C());
    surfaceScalarField ghf("ghf", g() & mesh.Cf());

    if (correctPhi_ && meshChanged)
    {
        #include "correctPhi.foamextend.H"
    }

    // Make the fluxes relative to the mesh motion
    fvc::makeRelative(phi, U);

    if (checkMeshCourantNo_)
    {
        #include "meshCourantNo.H"
    }

    // Pressure-velocity corrector
    while (pimple.loop())
    {
        twoPhaseProperties.correct();

        #include "alphaEqnSubCycle.foamextend.H"

        #include "UEqn.foamextend.H"

        // --- PISO loop
        while (pimple.correct())
        {
            #include "pEqn.foamextend.H"
        }

        p = pd + rho*gh;

        if (pd.needReference())
        {
            p += dimensionedScalar
            (
                "p",
                p.dimensions(),
                pRefValue - getRefCellValue(p, pdRefCell)
            );
        }

        turbulence->correct();

        // Update gradient fields
        gradp() = fvc::grad(p);

        gradU() = fvc::grad(U);
    }

    // Make the fluxes absolute for when runTime++ is called
    fvc::makeAbsolute(phi, U);

    return 0;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fluidModels
} // End namespace Foam

// ************************************************************************* //
