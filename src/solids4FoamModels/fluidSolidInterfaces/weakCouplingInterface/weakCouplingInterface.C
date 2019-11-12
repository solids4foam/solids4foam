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

#include "weakCouplingInterface.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace fluidSolidInterfaces
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(weakCouplingInterface, 0);
addToRunTimeSelectionTable
(
    physicsModel, weakCouplingInterface, fluidSolidInteraction
);
addToRunTimeSelectionTable
(
    fluidSolidInterface, weakCouplingInterface, dictionary
);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

weakCouplingInterface::weakCouplingInterface
(
    Time& runTime,
    const word& region
)
:
    fluidSolidInterface(typeName, runTime, region),
    solidZonesTractionPtrList_(),
    solidZonesTractionPrev_(nGlobalPatches()),
    predictedSolidZonesTractionPtrList_(),
    relaxationFactor_
    (
        fsiProperties().lookupOrDefault<scalar>("relaxationFactor", 0.01)
    )
{
    if (solidZonesTractionPtrList_.size())
    {
        FatalErrorIn
        (
            "weakCouplingInterface::weakCouplingInterface\n"
            "(\n"
            "    Time& runTime,\n"
            "    const word& region\n"
            ")\n"
        )   << "solidZonesTraction List already exists"
            << abort(FatalError);
    }

    solidZonesTractionPtrList_.setSize(nGlobalPatches());

    if (predictedSolidZonesTractionPtrList_.size())
    {
        FatalErrorIn
        (
            "weakCouplingInterface::weakCouplingInterface\n"
            "(\n"
            "    Time& runTime,\n"
            "    const word& region\n"
            ")\n"
        )   << "predictedSolidZonesTraction List already exists"
            << abort(FatalError);
    }

    predictedSolidZonesTractionPtrList_.setSize(nGlobalPatches());

    forAll(fluid().globalPatches(), interfaceI)
    {
        // Initialize solid zone traction fields
        const word solidZoneTractionName
        (
            "solidZoneTraction" + Foam::name(interfaceI)
        );

        IOobject solidZoneTractionHeader
        (
            solidZoneTractionName,
            runTime.timeName(),
            fluidMesh(),
            IOobject::MUST_READ
        );

        if (solidZoneTractionHeader.headerOk())
        {
            Info<< "Reading solidZoneTraction for global patch "
                << fluidMesh().boundary()
                   [
                       fluid().globalPatches()[interfaceI].patch().index()
                   ].name()
                << " from disk..." << endl;

            solidZonesTractionPtrList_.set
            (
                interfaceI,
                new vectorIOField
                (
                    IOobject
                    (
                        solidZoneTractionName,
                        runTime.timeName(),
                        fluidMesh(),
                        IOobject::MUST_READ,
                        IOobject::AUTO_WRITE
                    )
                )
            );
        }
        else
        {
            solidZonesTractionPtrList_.set
            (
                interfaceI,
                new vectorIOField
                (
                    IOobject
                    (
                        solidZoneTractionName,
                        runTime.timeName(),
                        fluidMesh(),
                        IOobject::NO_READ,
                        IOobject::AUTO_WRITE
                    ),
                    vectorField
                    (
                        solid().globalPatches()[interfaceI].globalPatch().size(),
                        vector::zero
                    )
                )
            );
        }

        solidZonesTractionPrev_[interfaceI] =
            vectorField
            (
                solid().globalPatches()[interfaceI].globalPatch().size(),
                vector::zero
            );

        const word predictedSolidZoneTractionName
        (
            "predictedSolidZoneTraction" + Foam::name(interfaceI)
        );

        IOobject predictedSolidZoneTractionHeader
        (
            predictedSolidZoneTractionName,
            runTime.timeName(),
            fluidMesh(),
            IOobject::MUST_READ
        );

        if (predictedSolidZoneTractionHeader.headerOk())
        {
            Info<< "Reading predictedSolidZoneTraction for global patch "
                << fluidMesh().boundary()
                   [
                       fluid().globalPatches()[interfaceI].patch().index()
                   ].name()
                << " from disk..." << endl;

            predictedSolidZonesTractionPtrList_.set
            (
                interfaceI,
                new vectorIOField
                (
                    IOobject
                    (
                        predictedSolidZoneTractionName,
                        runTime.timeName(),
                        fluidMesh(),
                        IOobject::MUST_READ,
                        IOobject::AUTO_WRITE
                    )
                )
            );
        }
        else
        {
            predictedSolidZonesTractionPtrList_.set
            (
                interfaceI,
                new vectorIOField
                (
                    IOobject
                    (
                        predictedSolidZoneTractionName,
                        runTime.timeName(),
                        fluidMesh(),
                        IOobject::NO_READ,
                        IOobject::AUTO_WRITE
                    ),
                    vectorField
                    (
                        solid().globalPatches()[interfaceI].globalPatch().size(),
                        vector::zero
                    )
                )
            );
        }
    }
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool weakCouplingInterface::evolve()
{
    initializeFields();

    updateInterpolatorAndGlobalPatches();

    solid().evolve();

    updateWeakDisplacement();

    moveFluidMesh();

    fluid().evolve();

    updateWeakTraction();

    solid().updateTotalFields();

    return 0;
}


void weakCouplingInterface::initializeFields()
{
    forAll(solid().globalPatches(), interfaceI)
    {
        predictedSolidZonesTractionPtrList_[interfaceI] =
            vectorField
            (
                solid().globalPatches()[interfaceI].globalPatch().size(),
                vector::zero
            );
    }

    fluidSolidInterface::initializeFields();
}


bool weakCouplingInterface::updateWeakDisplacement()
{
    List<scalar> minResidual(nGlobalPatches(), scalar(0));

    forAll(solid().globalPatches(), interfaceI)
    {
        const vectorField solidZonePointsDisplAtSolid =
            solid().faceZonePointDisplacementIncrement(interfaceI);

        solidZonesPointsDispls()[interfaceI] =
            ggiInterpolators()[interfaceI].slaveToMasterPointInterpolate
            (
                solidZonePointsDisplAtSolid
            );

        const vectorField solidZonePointsTotDisplAtSolid =
            solid().faceZonePointDisplacementOld(interfaceI);

        const vectorField solidZonePointsTotDispl =
            ggiInterpolators()[interfaceI].slaveToMasterPointInterpolate
            (
                solidZonePointsTotDisplAtSolid
            );

        residualsPrev()[interfaceI] = residuals()[interfaceI];

        residuals()[interfaceI] =
            solidZonesPointsDispls()[interfaceI]
          - fluidZonesPointsDispls()[interfaceI];

        fluidZonesPointsDisplsPrev()[interfaceI] =
            fluidZonesPointsDispls()[interfaceI];

        fluidZonesPointsDispls()[interfaceI] += residuals()[interfaceI];

        // Calculate residual norm
        scalar residualNorm = Foam::sqrt(gSum(magSqr(residuals()[interfaceI])));
        scalar residualNorm_2 = residualNorm;

        if (residualNorm > maxResidualsNorm()[interfaceI])
        {
            maxResidualsNorm()[interfaceI] = residualNorm;
        }

        residualNorm /= maxResidualsNorm()[interfaceI] + SMALL;

        Info<< "Current fsi relative residual norm ("
            << solidMesh().boundary()
               [
                   solid().globalPatches()[interfaceI].patch().index()
               ].name()
            << "): " << residualNorm << endl;

        interfacesPointsDisplsPrev()[interfaceI] =
            interfacesPointsDispls()[interfaceI];

        interfacesPointsDispls()[interfaceI] =
            solidZonesPointsDispls()[interfaceI];

        const vectorField intTotDispl =
            interfacesPointsDispls()[interfaceI] + solidZonePointsTotDispl;

        const scalar intTotDisplNorm = Foam::sqrt(gSum(magSqr(intTotDispl)));

        if (intTotDisplNorm > maxIntsDisplsNorm()[interfaceI])
        {
            maxIntsDisplsNorm()[interfaceI] = intTotDisplNorm;
        }

        residualNorm_2 /= maxIntsDisplsNorm()[interfaceI] + SMALL;

        Info<< "Alternative fsi residual ("
            << solidMesh().boundary()
               [
                   solid().globalPatches()[interfaceI].patch().index()
               ].name()
            << "): " << residualNorm_2 << endl;

        minResidual[interfaceI] = min(residualNorm_2, residualNorm);
    }

    // Make sure that displacement on all processors is equal to one
    // calculated on master processor
    fluidSolidInterface::syncFluidZonePointsDispl(fluidZonesPointsDispls());


    // Update elasticWallPressure boundary conditions, if found
    fluidSolidInterface::updateElasticWallPressureAcceleration();

    return max(minResidual);
}


void weakCouplingInterface::updateWeakTraction()
{
    Info<< "Update weak traction on solid patch/patches" << endl;

    List<vectorField> fluidZonesTractionAtSolid
    (
        nGlobalPatches(), vectorField()
    );

    List<vectorField> fluidZonesTraction
    (
        nGlobalPatches(), vectorField()
    );

    forAll(fluid().globalPatches(), interfaceI)
    {
        // Calc fluid traction
        const vectorField fluidZoneTraction =
            fluid().faceZoneViscousForce(interfaceI);

        const scalarField fluidZonePressure =
            fluid().faceZonePressureForce(interfaceI);

        const vectorField& p =
            fluid().globalPatches()[interfaceI].globalPatch().localPoints();

        const faceList& f =
            fluid().globalPatches()[interfaceI].globalPatch().localFaces();

        vectorField n(f.size(), vector::zero);
        forAll(n, faceI)
        {
            n[faceI] = f[faceI].normal(p);
            n[faceI] /= mag(n[faceI]);
        }

        fluidZonesTraction[interfaceI] =
            fluidZoneTraction - fluidZonePressure*n;

        fluidZonesTractionAtSolid[interfaceI] =
            vectorField
            (
                solid().globalPatches()[interfaceI].globalPatch().size(),
                vector::zero
            );

        fluidZonesTractionAtSolid[interfaceI] =
            ggiInterpolators()[interfaceI].masterToSlave
            (
                -fluidZonesTraction[interfaceI]
            );

        const scalar beta_ = relaxationFactor_;

        solidZonesTractionPrev_[interfaceI] =
            solidZonesTractionPtrList_[interfaceI];

        solidZonesTractionPtrList_[interfaceI] =
            beta_*fluidZonesTractionAtSolid[interfaceI]
          + (1.0 - beta_)*predictedSolidZonesTractionPtrList_[interfaceI];

        FatalError
            << "Check beta here!" << abort(FatalError);

        predictedSolidZonesTractionPtrList_[interfaceI] =
            2.0*solidZonesTractionPtrList_[interfaceI]
          - solidZonesTractionPrev_[interfaceI];
    }

    if (coupled())
    {
        Info<< "Setting weak traction on solid patch/patches" << endl;

        solid().setTraction
        (
            solidPatchIndices(),
            predictedSolidZonesTractionPtrList_
        );
    }

    // Total force at the fluid side of the interface
    {
        forAll(fluid().globalPatches(), interfaceI)
        {
            const vectorField& p =
                fluid().globalPatches()[interfaceI].globalPatch().localPoints();

            const faceList& f =
                fluid().globalPatches()[interfaceI].globalPatch().localFaces();

            vectorField S(f.size(), vector::zero);

            forAll(S, faceI)
            {
                S[faceI] = f[faceI].normal(p);
            }

            const vector totalTractionForce =
                sum(fluidZonesTraction[interfaceI]*mag(S));

            Info<< "Total force on interface patch "
                << fluidMesh().boundary()
                   [
                       fluid().globalPatches()[interfaceI].patch().index()
                   ].name()
                << " (fluid) = " << totalTractionForce << endl;
        }
    }

    // Total force at the solid side of the interface
    {
        forAll(fluid().globalPatches(), interfaceI)
        {
            const vectorField& p =
                solid().globalPatches()[interfaceI].globalPatch().localPoints();
            const faceList& f =
                solid().globalPatches()[interfaceI].globalPatch().localFaces();

            vectorField S(f.size(), vector::zero);

            forAll(S, faceI)
            {
                S[faceI] = f[faceI].normal(p);
            }

            const vector totalTractionForce =
                sum(fluidZonesTractionAtSolid[interfaceI]*mag(S));

            Info<< "Total force on interface patch "
                << solidMesh().boundary()
                   [
                       solid().globalPatches()[interfaceI].patch().index()
                   ].name()
                << " (solid) = " << totalTractionForce << endl;
        }
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fluidSolidInterfaces

} // End namespace Foam

// ************************************************************************* //
