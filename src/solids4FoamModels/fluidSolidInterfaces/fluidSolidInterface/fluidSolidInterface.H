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

Class
    fluidSolidInterface

Description
    Virtual base class for fluid-solid interface coupling.

Author
    Zeljko Tukovic, FSB Zagreb.  All rights reserved.
    Philip Cardiff, UCD. All rights reserved.

SourceFiles
    fluidSolidInterface.C
    newFluidSolidInterface.C

\*---------------------------------------------------------------------------*/

#ifndef fluidSolidInterface_H
#define fluidSolidInterface_H

#include "physicsModel.H"
#include "IOdictionary.H"
#include "dynamicFvMesh.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "fluidModel.H"
#include "solidModel.H"
#include "dynamicFvMesh.H"
#include "ggiInterpolation.H"
#include "RBFInterpolation.H"

using namespace rbf;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class fluidSolidInterface Declaration
\*---------------------------------------------------------------------------*/

class fluidSolidInterface
:
    public physicsModel,
    public IOdictionary
{
    // Private data

        //- FSI properties dictionary
        dictionary fsiProperties_;

        //- Flow solver
        autoPtr<fluidModel> fluid_;

        //- Solid solver
        autoPtr<solidModel> solid_;

        //- Solid interface patch index
        label solidPatchIndex_;

        //- Fluid interface patch index
        label fluidPatchIndex_;

        //- Solid face zone current points
        mutable vectorField* currentSolidZonePointsPtr_;

        //- Solid face zone primitive patch in current configuration
        mutable standAlonePatch* currentSolidZonePatchPtr_;

        //- Fluid zone to solid zone interpolator
        mutable std::shared_ptr<RBFInterpolation> rbfFluidToSolidPtr_;

        //- Solid zone to fluid zone interpolator
        mutable std::shared_ptr<RBFInterpolation> rbfSolidToFluidPtr_;

        //- Ggi zone-to-zone interpolation
        mutable GGIInterpolation
        <
            standAlonePatch,
            standAlonePatch
        >* ggiInterpolatorPtr_;

        //- Face map for conformal fluid and solid interfaces
        //  For each solid face, give the index of the corresponding fluid face
        mutable autoPtr<labelIOList> fluidToSolidFaceMapPtr_;

        //- Face map for conformal fluid and solid interfaces
        //  For each fluid face, give the index of the corresponding solid face
        mutable autoPtr<labelIOList> solidToFluidFaceMapPtr_;

        //- Point map for conformal fluid and solid interfaces
        //  For each solid point, give the index of the corresponding fluid
        //  point
        mutable autoPtr<labelIOList> fluidToSolidPointMapPtr_;

        //- Point map for conformal fluid and solid interfaces
        //  For each fluid point, give the index of the corresponding solid
        //  point
        mutable autoPtr<labelIOList> solidToFluidPointMapPtr_;

        //- Outer correction loop stoping tolerance
        const scalar outerCorrTolerance_;

        //- Maximal number of outer correctors
        const label nOuterCorr_;

        //- Fsi coupling on/off
        Switch coupled_;

        //- Optional: enable coupling after specified time
        const scalar couplingStartTime_;

        //- Predictor on/off
        const Switch predictor_;

        //- Interface deformation limit
        const scalar interfaceDeformationLimit_;

        //- Fluid zone point displacement
        vectorField fluidZonePointsDispl_;

        //- Fluid zone ref. point displacement
        vectorField fluidZonePointsDisplRef_;

        //- Fluid zone previous point displacement
        vectorField fluidZonePointsDisplPrev_;

        //- Solid zone point displacement
        vectorField solidZonePointsDispl_;

        //- Solid zone ref. point displacement
        vectorField solidZonePointsDisplRef_;

        //- Solid zone point displacement
        vectorField interfacePointsDispl_;

        //- Solid zone ref. point displacement
        vectorField interfacePointsDisplPrev_;

        //- Current FSI residual
        vectorField residual_;

        //- Previous FSI residual
        vectorField residualPrev_;

        //- Maximal resudual norm
        scalar maxResidualNorm_;

        //- Maximal interface displacement norm
        scalar maxIntDisplNorm_;

        //- Outer corrector
        label outerCorr_;

        //- Switch to enable writing o FSI residual to file
        const Switch writeResidualsToFile_;

        //- Residual file pointer
        autoPtr<OFstream> residualFilePtr_;

        //- Interpolator update frequency
        label interpolatorUpdateFrequency_;

        //- Accumulated fluid side interface displacement
        mutable vectorIOField* accumulatedFluidInterfaceDisplacementPtr_;

        //- Min edge length for interface points at fluid side
        mutable scalarField* minEdgeLengthPtr_;

    // Private Member Functions

        //- Update the coupled switch based on the couplingStartTime
        bool updateCoupled();

        //- Calculate current solid zone points
        void calcCurrentSolidZonePoints() const;

        //- Calculate current solid zone primitive patch
        void calcCurrentSolidZonePatch() const;

        //- Calculate RBF fluid to solid interpolator
        void calcRbfFluidToSolidInterpolator() const;

        //- Calculate RBF fluid to solid interpolator
        void calcRbfSolidToFluidInterpolator() const;

        //- Calculate fluid-to-solid / solid-to-fluid ggi interpolator
        void calcGgiInterpolator() const;

        //- Calculate fluid to solid face map
        void calcFluidToSolidFaceMap() const;

        //- Return a reference to the fluid to solid face map
        const labelList& fluidToSolidFaceMap() const;

        //- Calculate solid to fluid face map
        void calcSolidToFluidFaceMap() const;

        //- Return a reference to the solid to fluid face map
        const labelList& solidToFluidFaceMap() const;

        //- Calculate fluid to solid point map
        void calcFluidToSolidPointMap() const;

        //- Return a reference to the fluid to solid point map
        const labelList& fluidToSolidPointMap() const;

        //- Calculate solid to fluid point map
        void calcSolidToFluidPointMap() const;

        //- Return a reference to the solid to fluid point map
        const labelList& solidToFluidPointMap() const;

        //- Accumulated fluid interface displacement
        void calcAccumulatedFluidInterfaceDisplacement() const;

        //- Calculate minimal edge lengths
        void calcMinEdgeLength() const;

        //- Return accumulated interface displacement
        vectorIOField& accumulatedFluidInterfaceDisplacement();

        //- Return minimal edge length
        const scalarField& minEdgeLength() const;

        //- Disallow default bitwise copy construct
        fluidSolidInterface(const fluidSolidInterface&);

        //- Disallow default bitwise assignment
        void operator=(const fluidSolidInterface&);


public:

    //- Runtime type information
    TypeName("fluidSolidInterface");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            fluidSolidInterface,
            dictionary,
            (
                Time& runTime,
                const word& region
            ),
            (runTime, region)
        );


    // Constructors

        //- Construct from components
        fluidSolidInterface
        (
            const word& type,
            Time& runTime,
            const word& region = dynamicFvMesh::defaultRegion
        );


    // Selectors

        //- Select constructed from fluid and solid meshes
        static autoPtr<fluidSolidInterface> New
        (
            Time& runTime,
            const word& region
        );


    // Destructor

        virtual ~fluidSolidInterface();


    // Public enumerations

        //- Method for transferring information between the solid and the fluid
        //  interfaces
        enum interfaceTransferMethod
        {
            directMap,
            RBF,
            GGI
        };

        //- interfaceTransferMethod names
        static const NamedEnum<interfaceTransferMethod, 3>
            interfaceTransferMethodNames_;

private:
        //- Transfer method
        interfaceTransferMethod transferMethod_;

public:


    // Member Functions

        // Access

            //- Return reference to the fluid mesh
            dynamicFvMesh& fluidMesh()
            {
                return fluid_().mesh();
            }

            //- Return const reference to the fluid mesh
            const dynamicFvMesh& fluidMesh() const
            {
                return fluid_().mesh();
            }

            //- Return reference to the solid mesh
            dynamicFvMesh& solidMesh()
            {
                return solid_().mesh();
            }

            //- Return const reference to the solid mesh
            const dynamicFvMesh& solidMesh() const
            {
                return solid_().mesh();
            }

            //- Return time
            // to-do: consider storing a reference to time
            const Time& runTime() const
            {
                return fluid_().mesh().time();
            }

            //- Return sub-dictionary
            const dictionary& fsiProperties() const
            {
                return fsiProperties_;
            }

            //- Return fluid model
            const fluidModel& fluid() const
            {
                return fluid_();
            }

            //- Return fluid model
            fluidModel& fluid()
            {
                return fluid_();
            }

            //- Return solid model
            const solidModel& solid() const
            {
                return solid_();
            }

            //- Return solid model
            solidModel& solid()
            {
                return solid_();
            }

            //- Return solid patch index
            label solidPatchIndex() const
            {
                return solidPatchIndex_;
            }

            //- Return fluid patch index
            label fluidPatchIndex() const
            {
                return fluidPatchIndex_;
            }

            //- Return current solid zone points
            const vectorField& currentSolidZonePoints() const;

            //- Return current solid zone patch
            const standAlonePatch& currentSolidZonePatch() const;

            //- Return fluid to solid interpolator
            const std::shared_ptr<RBFInterpolation>& rbfFluidToSolid() const;

            //- Return fluid to solid interpolator
            const GGIInterpolation
            <
                standAlonePatch, standAlonePatch
            >& ggiInterpolator() const;

            //- Return fluid to solid interpolator
            const std::shared_ptr<RBFInterpolation>& rbfSolidToFluid() const;

            //- Transfer/map/interpolate from one zone faces to another zone
            //  faces
            template<class Type>
            void transferFacesZoneToZone
            (
                const word& fromRegion,          // from region name
                const word& toRegion,            // to region name
                const standAlonePatch& fromZone, // from zone
                const standAlonePatch& toZone,   // to zone
                const Field<Type>& fromField,    // from field
                Field<Type>& toField             // to field
            ) const;

            //- Transfer/map/interpolate from one zone points to another zone
            //  points
            template<class Type>
            void transferPointsZoneToZone
            (
                const word& fromRegion,          // from region name
                const word& toRegion,            // to region name
                const standAlonePatch& fromZone, // from zone
                const standAlonePatch& toZone,   // to zone
                const Field<Type>& fromField,    // from field
                Field<Type>& toField             // to field
            ) const;

            //- Return outer corrector loop tolerance
            scalar outerCorrTolerance() const
            {
                return outerCorrTolerance_;
            }

            //- Return max numter of outer correctors
            label nOuterCorr() const
            {
                return nOuterCorr_;
            }

            //- Return reference to writeResidualToFile switch
            const Switch& writeResidualsToFile() const
            {
                return writeResidualsToFile_;
            }

            //- Return reference to the residualFile
            OFstream& residualFile();

            //- Is it fluid and solid coupled
            const Switch& coupled() const
            {
                return coupled_;
            }

            //- Is it fluid and solid coupled
            Switch& coupled()
            {
                return coupled_;
            }

            //- Return predictor switch
            const Switch& predictor() const
            {
                return predictor_;
            }

            //- Return relaxation factor
            scalar interfaceDeformationLimit() const
            {
                return interfaceDeformationLimit_;
            }

            //- Return fluid zone point displacement
            vectorField& fluidZonePointsDispl()
            {
                return fluidZonePointsDispl_;
            }

            //- Return fluid zone ref. point displacement
            vectorField& fluidZonePointsDisplRef()
            {
                return fluidZonePointsDisplRef_;
            }

            //- Return fluid zone previous point displacement
            vectorField& fluidZonePointsDisplPrev()
            {
                return fluidZonePointsDisplPrev_;
            }

            //- Return solid zone point displacement
            vectorField& solidZonePointsDispl()
            {
                return solidZonePointsDispl_;
            }

            //- Return solid zone ref. point displacement
            vectorField& solidZonePointsDisplRef()
            {
                return solidZonePointsDisplRef_;
            }

            //- Return solid zone point displacement
            vectorField& interfacePointsDispl()
            {
                return interfacePointsDispl_;
            }

            //- Return solid zone ref. point displacement
            vectorField& interfacePointsDisplPrev()
            {
                return interfacePointsDisplPrev_;
            }

            //- Return current fsi residual
            vectorField& residual()
            {
                return residual_;
            }

            //- Return previous fsi residual
            vectorField& residualPrev()
            {
                return residualPrev_;
            }

            //- Return current outer iteration
            label& outerCorr()
            {
                return outerCorr_;
            }

            //- Return const current outer iteration
            const label& outerCorr() const
            {
                return outerCorr_;
            }

            //- Return maximal resudual norm
            scalar& maxResidualNorm()
            {
                return maxResidualNorm_;
            }

            //- Return const maximal resudual norm
            const scalar& maxResidualNorm() const
            {
                return maxResidualNorm_;
            }

            //- Return interface displacement norm
            scalar& maxIntDisplNorm()
            {
                return maxIntDisplNorm_;
            }

            //- Return const interface displacement norm
            const scalar& maxIntDisplNorm() const
            {
                return maxIntDisplNorm_;
            }


        // Edit

            //- Update the size of the time-step
            virtual void setDeltaT(Time& runTime);

            //- Evolve the interface
            virtual bool evolve() = 0;

            //- Initialize fields
            virtual void initializeFields();

            //- Update interpolator
            virtual void updateInterpolatorAndGlobalPatches();

            //- Move fluid mesh
            virtual void moveFluidMesh();

            //- Update interface force
            virtual void updateForce();

            //- Update interface force
            virtual scalar updateResidual();

            //- Update movingWallPressure boundary condition acceleration, if
            //  found
            virtual void updateMovingWallPressureAcceleration();

            //- Update elasticWallPressure boundary condition acceleration, if
            //  found
            virtual void updateElasticWallPressureAcceleration();

            //- Sync the fluidZonePointsDispl field across processors
            virtual void syncFluidZonePointsDispl
            (
                vectorField& fluidZonePointsDispl
            );

            //- Write fields
            virtual void writeFields(const Time& runTime);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#ifdef NoRepository
#   include "fluidSolidInterfaceTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
