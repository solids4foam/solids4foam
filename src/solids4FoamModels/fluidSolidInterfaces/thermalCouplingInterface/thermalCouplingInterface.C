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

#include "thermalCouplingInterface.H"
#include "addToRunTimeSelectionTable.H"
#include "movingWallPressureFvPatchScalarField.H"
#include "elasticWallPressureFvPatchScalarField.H"
#include "OSspecific.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace fluidSolidInterfaces
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(thermalCouplingInterface, 0);
addToRunTimeSelectionTable
(
    fluidSolidInterface, thermalCouplingInterface, dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

thermalCouplingInterface::thermalCouplingInterface
(
    Time& runTime,
    const word& region
)
:
    fluidSolidInterface(typeName, runTime, region),
    relaxationFactor_
    (
        fsiProperties().lookupOrDefault<scalar>("relaxationFactor", 0.01)
    ),
    predictTemperatureAndHeatFlux_
    (
        fsiProperties().lookupOrDefault<bool>
        (
            "predictTemperatureAndHeatFlux",
            false
        )
    ),
    thermalResidualFilePtr_(),
    maxThermalResidualsNorm_(),
    oldSolidFaceZoneTemperature_(),
    oldSolidFaceZoneHeatFlux_(),
    oldOldSolidFaceZoneTemperature_(),
    oldOldSolidFaceZoneHeatFlux_(),
    timeIndex_(-1)
{
    maxThermalResidualsNorm_.setSize(nGlobalPatches());

    // Set equivalent interface heat transfer coefficient
    setEqInterHeatTransferCoeff();

    // Initialize old-times face-zone temperatures and heat-fluxes
    {
        oldSolidFaceZoneTemperature_.setSize(nGlobalPatches());
        oldSolidFaceZoneHeatFlux_.setSize(nGlobalPatches());
        oldOldSolidFaceZoneTemperature_.setSize(nGlobalPatches());
        oldOldSolidFaceZoneHeatFlux_.setSize(nGlobalPatches());
    
        forAll(oldSolidFaceZoneTemperature_, interI)
        {
            const standAlonePatch& solidZone =
                solid().globalPatches()[interI].globalPatch();
        
#ifdef OPENFOAMESIORFOUNDATION
            oldSolidFaceZoneTemperature_.set
            (
                interI,
                new scalarField(solidZone.size(), 0)
            );
#else
            oldSolidFaceZoneTemperature_.set
            (
                interI,
                scalarField(solidZone.size(), 0)
            );
#endif
            
            oldSolidFaceZoneHeatFlux_.set
            (
                interI,
                new scalarField(solidZone.size(), 0)
            );
            
            oldOldSolidFaceZoneTemperature_.set
            (
                interI,
                new scalarField(solidZone.size(), 0)
            );
            
            oldOldSolidFaceZoneHeatFlux_.set
            (
                interI,
                new scalarField(solidZone.size(), 0)
            );
        }
    }
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool thermalCouplingInterface::evolve()
{
    initializeFields();

    updateInterpolatorAndGlobalPatches();

    scalar residualNorm = 0;

    do
    {
        outerCorr()++;

        Info<< nl << "Time = " << fluid().runTime().timeName()
            << ", iteration: " << outerCorr() << endl;

        // Transfer temperature and heat flux from the solid to the fluid
        updateHeatFluxAndTemperatureOnFluidInterface();

        // Solve fluid
        fluid().evolve();

        // Transfer temperature and heat flux from the fluid to the solid
        updateHeatFluxAndTemperatureOnSolidInterface();

        // Solve solid
        solid().evolve();

        // Calculate the thermal residual
        residualNorm = calcThermalResidual();

        // Optional: write residuals to file
        if (writeResidualsToFile() && Pstream::master())
        {
            thermalResidualFile()
                << runTime().value() << " "
                << outerCorr() << " "
                << residualNorm << endl;
        }
    }
    while (residualNorm > outerCorrTolerance() && outerCorr() < nOuterCorr());

    return 0;
}



void thermalCouplingInterface::updateHeatFluxAndTemperatureOnFluidInterface()
{
    Info<< "Setting heat flux and temperature on fluid interfaces" << endl;

    for (label interfaceI = 0; interfaceI < nGlobalPatches(); interfaceI++)
    {
        // Take references to zones
        const standAlonePatch& fluidZone =
            fluid().globalPatches()[interfaceI].globalPatch();
        const standAlonePatch& solidZone =
            solid().globalPatches()[interfaceI].globalPatch();

        // Calculate temperature of solid zone
        scalarField solidZoneTemperature
        (
            solid().faceZoneTemperature(interfaceI)
        );

        // Calculate heat flux of solid zone
        scalarField solidZoneHeatFlux
        (
            solid().faceZoneHeatFlux(interfaceI)
        );

        if (timeIndex_ < runTime().timeIndex())
        {
            oldOldSolidFaceZoneTemperature_[interfaceI] =
                oldSolidFaceZoneTemperature_[interfaceI];
            oldOldSolidFaceZoneHeatFlux_[interfaceI] =
                oldSolidFaceZoneHeatFlux_[interfaceI];

            oldSolidFaceZoneTemperature_[interfaceI] =
                solidZoneTemperature;
            oldSolidFaceZoneHeatFlux_[interfaceI] =
                solidZoneHeatFlux;

            if
            (
                predictTemperatureAndHeatFlux_
             && (runTime().timeIndex() > 2)
            )
            {
                Info << "Predicting temperature and heat-flux" << endl;

                // Linear extrapolation of temp. and heat-flux in time
                solidZoneTemperature =
                    2*oldSolidFaceZoneTemperature_[interfaceI]
                  - oldOldSolidFaceZoneTemperature_[interfaceI];
                solidZoneHeatFlux =
                    2*oldSolidFaceZoneHeatFlux_[interfaceI]
                  - oldOldSolidFaceZoneHeatFlux_[interfaceI];
            }
            
            timeIndex_ = runTime().timeIndex();
        }
        
        // Initialise the fluid zone temperature field
        // that is to be interpolated from the solid zone
        scalarField fluidZoneTemperature(fluidZone.size(), 0);

        // Initialise the fluid zone heat flux field
        // that is to be interpolated from the solid zone
        scalarField fluidZoneHeatFlux(fluidZone.size(), 0);

        // Transfer the field from the solid interface to the fluid interface
        interfaceToInterfaceList()[interfaceI].transferFacesZoneToZone
        (
            solidZone,                 // from zone
            fluidZone,                 // to zone
            solidZoneTemperature,      // from field
            fluidZoneTemperature       // to field
        );

        interfaceToInterfaceList()[interfaceI].transferFacesZoneToZone
        (
            solidZone,                 // from zone
            fluidZone,                 // to zone
            solidZoneHeatFlux,         // from field
            fluidZoneHeatFlux          // to field
        );
        fluidZoneHeatFlux *= -1;

        // Set temperature on fluid interface
        if (coupled())
        {
            fluid().setTemperatureAndHeatFlux
            (
                interfaceI,
                fluidPatchIndices()[interfaceI],
                fluidZoneTemperature,
                fluidZoneHeatFlux
            );
        }

        // Print total heat flux on solid and fluid interfaces
        Info<< "Heat flow rate on solid interface " << interfaceI << ": "
            << heatFlowRateOnInterface(solidZone, solidZoneHeatFlux)
            << nl << endl;
    }
}



void thermalCouplingInterface::updateHeatFluxAndTemperatureOnSolidInterface()
{
    Info<< "Setting heat flux and temperature on solid interfaces" << endl;

    for (label interfaceI = 0; interfaceI < nGlobalPatches(); interfaceI++)
    {
        // Take references to zones
        const standAlonePatch& fluidZone =
            fluid().globalPatches()[interfaceI].globalPatch();
        const standAlonePatch& solidZone =
            solid().globalPatches()[interfaceI].globalPatch();

        // Calculate temperature of fluid zone
        scalarField fluidZoneTemperature
        (
            fluid().faceZoneTemperature(interfaceI)
        );

        // Calculate heat flux of solid zone
        scalarField fluidZoneHeatFlux
        (
            fluid().faceZoneHeatFlux(interfaceI)
        );

        // Initialise the fluid zone temperature field
        // that is to be interpolated from the solid zone
        scalarField solidZoneTemperature(solidZone.size(), 0);

        // Initialise the fluid zone heat flux field
        // that is to be interpolated from the solid zone
        scalarField solidZoneHeatFlux(solidZone.size(), 0);

        // Transfer the field from the fluid interface
        // to the solid interface
        interfaceToInterfaceList()[interfaceI].transferFacesZoneToZone
        (
            fluidZone,                 // from zone
            solidZone,                 // to zone
            fluidZoneTemperature,      // from field
            solidZoneTemperature       // to field
        );
        interfaceToInterfaceList()[interfaceI].transferFacesZoneToZone
        (
            fluidZone,                 // from zone
            solidZone,                 // to zone
            fluidZoneHeatFlux,         // from field
            solidZoneHeatFlux          // to field
        );
        
        solidZoneHeatFlux *= -1;

        // Set temperature on fluid interface
        if (coupled())
        {
            solid().setTemperatureAndHeatFlux
            (
                interfaceI,
                solidPatchIndices()[interfaceI],
                solidZoneTemperature,
                solidZoneHeatFlux
            );
        }

        // Print total heat flux on solid and fluid interfaces
        Info<< "Heat flow rate on fluid interface " << interfaceI << ": "
            << heatFlowRateOnInterface(fluidZone, fluidZoneHeatFlux)
            << nl << endl;
    }
}


scalar thermalCouplingInterface::calcThermalResidual()
{
    // Maximum residual for all interfaces
    scalar maxResidual = 0;

    for (label interfaceI = 0; interfaceI < nGlobalPatches(); interfaceI++)
    {
        // Take references to zones
        const standAlonePatch& fluidZone =
            fluid().globalPatches()[interfaceI].globalPatch();
        const standAlonePatch& solidZone =
            solid().globalPatches()[interfaceI].globalPatch();

        // Calculate the face temperature of the solid interface
        const scalarField solidZoneTemperatureAtSolid
        (
            solid().faceZoneTemperature(interfaceI)
        );

        // Initialise the solid zone temperature field
        // interpolated to the fluid zone
        scalarField solidZoneTemperature(fluidZone.size(), 0);

        // Transfer displacement field from the solid to the fluid
        interfaceToInterfaceList()[interfaceI].transferFacesZoneToZone
        (
            solidZone,                              // from zone
            fluidZone,                              // to zone
            solidZoneTemperatureAtSolid,            // from field
            solidZoneTemperature                    // to field
        );

        // Calculate the face temperature of the fluid interface
        const scalarField fluidZoneTemperature
        (
            fluid().faceZoneTemperature(interfaceI)
        );

        // Calculate thermal residual
        scalarField residual =
            solidZoneTemperature
          - fluidZoneTemperature;


        // Calculate thermal resudal norm
        scalar residualNorm =
            Foam::sqrt(gSum(magSqr(residual)));

        if (residualNorm > maxThermalResidualsNorm_[interfaceI])
        {
            maxThermalResidualsNorm_[interfaceI] = residualNorm;
        }

        residualNorm /= maxThermalResidualsNorm_[interfaceI] + SMALL;

        Info<< "CHT relative residual norm for interface " << interfaceI
            << ": " << residualNorm << endl;

        // Update the maximum residual for all interfaces
        maxResidual = max(maxResidual, residualNorm);
    }

    return maxResidual;
}


OFstream& thermalCouplingInterface::thermalResidualFile()
{
    if (Pstream::parRun())
    {
        if (!Pstream::master())
        {
            FatalErrorIn
            (
                "OFstream& thermalCouplingInterface::thermalResidualFile()"
            )   << "Only the master processor can call this functon!"
                << abort(FatalError);
        }
    }

    if (thermalResidualFilePtr_.empty())
    {
        const fileName historyDir = runTime().path()/"residuals";
        mkDir(historyDir);
        thermalResidualFilePtr_.set
        (
            new OFstream(historyDir/"thermalResiduals.dat")
        );
        thermalResidualFilePtr_()
            << "Time outerCorrector residual" << endl;
    }

    return thermalResidualFilePtr_();
}


scalar thermalCouplingInterface::heatFlowRateOnInterface
(
    const standAlonePatch& zone, const scalarField& zoneHeatFlux
) const
{
    const vectorField& localPoints = zone.localPoints();
    const faceList& localFaces = zone.localFaces();

    // Calculate face area vectors
    vectorField S(localFaces.size(), vector::zero);
    forAll(S, faceI)
    {
#ifdef OPENFOAMFOUNDATION
        S[faceI] = localFaces[faceI].area(localPoints);
#else
        S[faceI] = localFaces[faceI].normal(localPoints);
#endif
    }

    // No need for global sum as the zone is already global
    return sum(zoneHeatFlux*mag(S));
}

void thermalCouplingInterface::setEqInterHeatTransferCoeff()
{
    Info<< "Setting equivalent interface heat transfeor coefficient" << endl;

    for (label interfaceI = 0; interfaceI < nGlobalPatches(); interfaceI++)
    {
        // Take references to zones
        const standAlonePatch& fluidZone =
            fluid().globalPatches()[interfaceI].globalPatch();
        const standAlonePatch& solidZone =
            solid().globalPatches()[interfaceI].globalPatch();

        // Calculate HTC of fluid zone
        scalarField fluidZoneHTC
        (
            fluid().faceZoneHeatTransferCoeff(interfaceI)
        );

        // Initialise the solid zone traction field that is to be interpolated
        // from the fluid zone
        scalarField solidZoneHTC(solidZone.size(), 0);

        // Transfer the field frm the fluid interface to the solid interface
        interfaceToInterfaceList()[interfaceI].transferFacesZoneToZone
        (
            fluidZone,                 // from zone
            solidZone,                 // to zone
            fluidZoneHTC,              // from field
            solidZoneHTC               // to field
        );

        // Add HTC from solid side and get equivalent HTC
        solidZoneHTC +=
            solid().faceZoneHeatTransferCoeff(interfaceI);

        // Set eq interface HTC on solid
        if (coupled())
        {
            solid().setEqInterHeatTransferCoeff
            (
                interfaceI,
                solidPatchIndices()[interfaceI],
                solidZoneHTC
            );
        }

        // Transfer the field frm the fluid interface to the solid interface
        interfaceToInterfaceList()[interfaceI].transferFacesZoneToZone
        (
            solidZone,                 // from zone
            fluidZone,                 // to zone
            solidZoneHTC,              // from field
            fluidZoneHTC               // to field
        );

        // Set eq interface HTC on fluid
        if (coupled())
        {
            fluid().setEqInterHeatTransferCoeff
            (
                interfaceI,
                fluidPatchIndices()[interfaceI],
                fluidZoneHTC
            );
        }

        // Print average eq interface HTC on fluid and solid
        Info<< "Avg eq. inter HTC on fluid interface " << interfaceI << ": "
            << average(fluidZoneHTC) << nl
            << "Avg eq. inter HTC on solid interface " << interfaceI << ": "
            << average(solidZoneHTC) << nl << endl;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fluidSolidInterfaces

} // End namespace Foam

// ************************************************************************* //
