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
    solidZonesTractionPrev_(fluid().globalPatches().size()),
    predictedSolidZonesTractionPtrList_(),
    relaxationFactor_
    (
        fsiProperties().lookupOrDefault<scalar>("relaxationFactor", 0.01)
    )
{
    solidZonesTractionPtrList_.setSize(fluid().globalPatches().size());

    predictedSolidZonesTractionPtrList_.setSize(fluid().globalPatches().size());

    forAll(fluid().globalPatches(), i)
    {
        if (solidZonesTractionPtrList_.set(i))
        {
            FatalErrorIn
            (
                "weakCouplingInterface::weakCouplingInterface\n"
                "(\n"
                "    Time& runTime,\n"
                "    const word& region\n"
                ")\n"
            )
                << "solidZoneTraction already set for global patch: "
                << fluidMesh().boundary()
                   [
                       fluid().globalPatches()[i].patch().index()
                   ].name()
                << "!" << abort(FatalError);
        }

        // Initialize solid zone traction fields
        const word solidZoneTractionName
        (
            "solidZoneTraction"
#if FOAMEXTEND > 40
          + name(i)
#else
          + word(std::to_string(i), false)
#endif
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
                       fluid().globalPatches()[i].patch().index()
                   ].name()
                << " from disk..." << endl;

            solidZonesTractionPtrList_.set
            (
                i,
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
                i,
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
                        solid().globalPatches()[i].globalPatch().size(),
                        vector::zero
                    )
                )
            );
        }

        solidZonesTractionPrev_[i] =
            vectorField
            (
                solid().globalPatches()[i].globalPatch().size(),
                vector::zero
            );

        if (predictedSolidZonesTractionPtrList_.set(i))
        {
            FatalErrorIn
            (
                "weakCouplingInterface::weakCouplingInterface\n"
                "(\n"
                "    Time& runTime,\n"
                "    const word& region\n"
                ")\n"
            )
                << "predictedSolidZoneTraction already set for global patch: "
                << fluidMesh().boundary()
                   [
                       fluid().globalPatches()[i].patch().index()
                   ].name()
                << "!" << abort(FatalError);
        }

        const word predictedSolidZoneTractionName
        (
            "predictedSolidZoneTraction"
#if FOAMEXTEND > 40
          + name(i)
#else
          + word(std::to_string(i), false)
#endif
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
                       fluid().globalPatches()[i].patch().index()
                   ].name()
                << " from disk..." << endl;

            predictedSolidZonesTractionPtrList_.set
            (
                i,
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
                i,
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
                        solid().globalPatches()[i].globalPatch().size(),
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
    forAll(solid().globalPatches(), i)
    {
        predictedSolidZonesTractionPtrList_[i] =
            vectorField
            (
                solid().globalPatches()[i].globalPatch().size(),
                vector::zero
            );
    }

    fluidSolidInterface::initializeFields();
}


bool weakCouplingInterface::updateWeakDisplacement()
{
    List<scalar> minResidual(solid().globalPatches().size(), scalar(0));

    const List<tmp<vectorField> > faceZonesPointDisplacementIncrement
    (
        solid().faceZonesPointDisplacementIncrement()
    );

    const List<tmp<vectorField> > faceZonesPointDisplacementOld
    (
        solid().faceZonesPointDisplacementOld()
    );

    forAll(solid().globalPatches(), i)
    {
        const vectorField& solidZonePointsDisplAtSolid =
            faceZonesPointDisplacementIncrement[i]();

        solidZonesPointsDispls()[i] =
            ggiInterpolators()[i].slaveToMasterPointInterpolate
            (
                solidZonePointsDisplAtSolid
            );

        const vectorField& solidZonePointsTotDisplAtSolid =
            faceZonesPointDisplacementOld[i]();

        const vectorField solidZonePointsTotDispl =
            ggiInterpolators()[i].slaveToMasterPointInterpolate
            (
                solidZonePointsTotDisplAtSolid
            );

        residualsPrev()[i] = residuals()[i];

        residuals()[i] = solidZonesPointsDispls()[i] - fluidZonesPointsDispls()[i];

        fluidZonesPointsDisplsPrev()[i] = fluidZonesPointsDispls()[i];

        fluidZonesPointsDispls()[i] += residuals()[i];

        // Calculate residual norm
        scalar residualNorm = Foam::sqrt(gSum(magSqr(residuals()[i])));
        scalar residualNorm_2 = residualNorm;

        if (residualNorm > maxResidualsNorm()[i])
        {
            maxResidualsNorm()[i] = residualNorm;
        }

        residualNorm /= maxResidualsNorm()[i] + SMALL;

        Info<< "Current fsi relative residual norm ("
            << solidMesh().boundary()
               [
                   solid().globalPatches()[i].patch().index()
               ].name()
            << "): " << residualNorm << endl;

        interfacesPointsDisplsPrev()[i] = interfacesPointsDispls()[i];

        interfacesPointsDispls()[i] = solidZonesPointsDispls()[i];

        const vectorField intTotDispl =
            interfacesPointsDispls()[i] + solidZonePointsTotDispl;

        const scalar intTotDisplNorm = Foam::sqrt(gSum(magSqr(intTotDispl)));

        if (intTotDisplNorm > maxIntsDisplsNorm()[i])
        {
            maxIntsDisplsNorm()[i] = intTotDisplNorm;
        }

        residualNorm_2 /= maxIntsDisplsNorm()[i] + SMALL;

        Info<< "Alternative fsi residual ("
            << solidMesh().boundary()
               [
                   solid().globalPatches()[i].patch().index()
               ].name()
            << "): " << residualNorm_2 << endl;

        minResidual[i] = min(residualNorm_2, residualNorm);
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

    const List<tmp<vectorField> > faceZonesViscousForce
    (
        fluid().faceZonesViscousForce()
    );

    const List<tmp<scalarField> > faceZonesPressureForce
    (
        fluid().faceZonesPressureForce()
    );

    List<vectorField> fluidZonesTractionAtSolid
    (
        solid().globalPatches().size(), vectorField()
    );

    List<vectorField> fluidZonesTraction
    (
        fluid().globalPatches().size(), vectorField()
    );

    forAll(fluid().globalPatches(), i)
    {
        // Calc fluid traction
        const vectorField& fluidZoneTraction = faceZonesViscousForce[i]();

        const scalarField& fluidZonePressure = faceZonesPressureForce[i]();

        const vectorField& p =
            fluid().globalPatches()[i].globalPatch().localPoints();

        const faceList& f =
            fluid().globalPatches()[i].globalPatch().localFaces();

        vectorField n(f.size(), vector::zero);
        forAll(n, faceI)
        {
            n[faceI] = f[faceI].normal(p);
            n[faceI] /= mag(n[faceI]);
        }

        fluidZonesTraction[i] =
            fluidZoneTraction - fluidZonePressure*n;

        fluidZonesTractionAtSolid[i] =
            vectorField
            (
                solid().globalPatches()[i].globalPatch().size(),
                vector::zero
            );

        fluidZonesTractionAtSolid[i] =
            ggiInterpolators()[i].masterToSlave
            (
                -fluidZonesTraction[i]
            );

        const scalar beta_ = relaxationFactor_;

        solidZonesTractionPrev_[i] = solidZonesTractionPtrList_[i];

        solidZonesTractionPtrList_[i] =
            beta_*fluidZonesTractionAtSolid[i]
          + (1.0 - beta_)*predictedSolidZonesTractionPtrList_[i];

        FatalError
            << "Check beta here!" << abort(FatalError);

        predictedSolidZonesTractionPtrList_[i] =
            2.0*solidZonesTractionPtrList_[i] - solidZonesTractionPrev_[i];
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
        forAll(fluid().globalPatches(), i)
        {
            const vectorField& p =
                fluid().globalPatches()[i].globalPatch().localPoints();

            const faceList& f =
                fluid().globalPatches()[i].globalPatch().localFaces();

            vectorField S(f.size(), vector::zero);

            forAll(S, faceI)
            {
                S[faceI] = f[faceI].normal(p);
            }

            const vector totalTractionForce = sum(fluidZonesTraction[i]*mag(S));

            Info<< "Total force on interface patch "
                << fluidMesh().boundary()
                   [
                       fluid().globalPatches()[i].patch().index()
                   ].name()
                << " (fluid) = " << totalTractionForce << endl;
        }
    }

    // Total force at the solid side of the interface
    {
        forAll(fluid().globalPatches(), i)
        {
            const vectorField& p =
                solid().globalPatches()[i].globalPatch().localPoints();
            const faceList& f =
                solid().globalPatches()[i].globalPatch().localFaces();

            vectorField S(f.size(), vector::zero);

            forAll(S, faceI)
            {
                S[faceI] = f[faceI].normal(p);
            }

            const vector totalTractionForce = sum(fluidZonesTractionAtSolid[i]*mag(S));

            Info<< "Total force on interface patch "
                << solidMesh().boundary()
                   [
                       solid().globalPatches()[i].patch().index()
                   ].name()
                << " (solid) = " << totalTractionForce << endl;
        }
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fluidSolidInterfaces

} // End namespace Foam

// ************************************************************************* //
