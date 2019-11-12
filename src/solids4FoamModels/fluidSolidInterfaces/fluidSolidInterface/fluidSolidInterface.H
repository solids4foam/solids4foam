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
    Extended to multiple interfaces by Danial Khazaei

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

        //- Solid interface patch names
        wordList solidPatchNames_;

        //- Fluid interface patch Names
        wordList fluidPatchNames_;

        //- Solid interface patch indices
        labelList solidPatchIndices_;

        //- Fluid interface patch indices
        labelList fluidPatchIndices_;

        //- Number of global poly patches
        label nGlobalPatches_;

        //- List of solid face zones current points
        mutable PtrList<vectorField> currentSolidZonesPointsList_;

        //- List of solid face zones primitive patches in current configuration
        mutable PtrList<standAlonePatch> currentSolidZonesPatchesList_;

        //- List of fluid zone to solid zone interpolators
        mutable List<std::shared_ptr<RBFInterpolation> > rbfFluidToSolidList_;

        //- List of solid zone to fluid zone interpolators
        mutable List<std::shared_ptr<RBFInterpolation> > rbfSolidToFluidList_;

        //- List of ggi zone-to-zone interpolators
        mutable PtrList
        <
            GGIInterpolation<standAlonePatch, standAlonePatch>
        > ggiInterpolatorsList_;

        //- List of face map for conformal fluid and solid interfaces
        //  For each solid face, give the index of the corresponding fluid face
        mutable PtrList<labelIOList> fluidToSolidFaceMapsList_;

        //- List of face map for conformal fluid and solid interfaces
        //  For each fluid face, give the index of the corresponding solid face
        mutable PtrList<labelIOList> solidToFluidFaceMapsList_;

        //- List of point map for conformal fluid and solid interfaces
        //  For each solid point, give the index of the corresponding fluid
        //  point
        mutable PtrList<labelIOList> fluidToSolidPointMapsList_;

        //- List of point map for conformal fluid and solid interfaces
        //  For each fluid point, give the index of the corresponding solid
        //  point
        mutable PtrList<labelIOList> solidToFluidPointMapsList_;

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

        //- List of fluid zones point displacements
        List<vectorField> fluidZonesPointsDispls_;

        //- List of fluid zones ref. point displacements
        List<vectorField> fluidZonesPointsDisplsRef_;

        //- List of fluid zones previous point displacements
        List<vectorField> fluidZonesPointsDisplsPrev_;

        //- List of solid zones point displacements
        List<vectorField> solidZonesPointsDispls_;

        //- List of solid zones ref. point displacements
        List<vectorField> solidZonesPointsDisplsRef_;

        //- List of interfaces point displacements
        List<vectorField> interfacesPointsDispls_;

        //- List of interfaces ref. point displacements
        List<vectorField> interfacesPointsDisplsPrev_;

        //- List of current FSI residuals
        List<vectorField> residuals_;

        //- List of Previous FSI residuals
        List<vectorField> residualsPrev_;

        //- List of maximal residuals norm
        List<scalar> maxResidualsNorm_;

        //- List of maximal interfaces displacements norm
        List<scalar> maxIntsDisplsNorm_;

        //- Outer corrector
        label outerCorr_;

        //- Switch to enable writing o FSI residual to file
        const Switch writeResidualsToFile_;

        //- Residual file pointer
        autoPtr<OFstream> residualFilePtr_;

        //- Interpolator update frequency
        label interpolatorUpdateFrequency_;

        //- List of accumulated fluid side interfaces displacements
        mutable PtrList<vectorIOField>
        accumulatedFluidInterfacesDisplacementsList_;

        //- List of min edge lengths for interfaces points at fluid side
        mutable PtrList<scalarField> minEdgeLengthsList_;

    // Private Member Functions

        //- Update the coupled switch based on the couplingStartTime
        bool updateCoupled();

        //- Calculate current solid zones points
        void calcCurrentSolidZonesPoints() const;

        //- Calculate current solid zone primitive patches
        void calcCurrentSolidZonesPatches() const;

        //- Calculate RBF fluid to solid interpolators
        void calcRbfFluidToSolidInterpolators() const;

        //- Calculate RBF fluid to solid interpolators
        void calcRbfSolidToFluidInterpolators() const;

        //- Calculate fluid-to-solid / solid-to-fluid ggi interpolators
        void calcGgiInterpolators() const;

        //- Calculate fluid to solid face maps
        void calcFluidToSolidFaceMaps() const;

        //- Return a reference to the fluid to solid face maps
        const PtrList<labelIOList>& fluidToSolidFaceMaps() const;

        //- Calculate solid to fluid face maps
        void calcSolidToFluidFaceMaps() const;

        //- Return a reference to the solid to fluid face maps
        const PtrList<labelIOList>& solidToFluidFaceMaps() const;

        //- Calculate fluid to solid point maps
        void calcFluidToSolidPointMaps() const;

        //- Return a reference to the fluid to solid point maps
        const PtrList<labelIOList>& fluidToSolidPointMaps() const;

        //- Calculate solid to fluid point maps
        void calcSolidToFluidPointMaps() const;

        //- Return a reference to the solid to fluid point maps
        const PtrList<labelIOList>& solidToFluidPointMaps() const;

        //- Accumulated fluid interfaces displacements
        void calcAccumulatedFluidInterfacesDisplacements() const;

        //- Return accumulated interfaces displacements
        PtrList<vectorIOField>& accumulatedFluidInterfacesDisplacements();

        //- Calculate minimal edge lengths
        void calcMinEdgeLengths() const;

        //- Return minimal edge lengths
        const PtrList<scalarField>& minEdgeLengths() const;

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

            //- Return solid patch indices
            const labelList& solidPatchIndices() const
            {
                return solidPatchIndices_;
            }

            //- Return fluid patch indices
            const labelList& fluidPatchIndices() const
            {
                return fluidPatchIndices_;
            }

            //- Return number of global poly patches
            const label& nGlobalPatches() const
            {
                return nGlobalPatches_;
            }

            //- Return current solid zones points
            const PtrList<vectorField>& currentSolidZonesPoints() const;

            //- Return current solid zones patches
            const PtrList<standAlonePatch>& currentSolidZonesPatches() const;

            //- Return fluid to solid interpolators
            const List<std::shared_ptr<RBFInterpolation> >&
            rbfFluidToSolid() const;

            //- Return fluid to solid interpolators
            const List<std::shared_ptr<RBFInterpolation> >&
            rbfSolidToFluid() const;

            //- Return fluid to solid interpolators
            const PtrList
            <
                GGIInterpolation<standAlonePatch, standAlonePatch>
            >& ggiInterpolators() const;

            //- Transfer/map/interpolate from one zone faces to another zone
            //  faces
            template<class Type>
            void transferFacesZonesToZones
            (
                const word& fromRegion,                    // from region name
                const word& toRegion,                      // to region name
                const PtrList<globalPolyPatch>& fromZones, // from zones
                const PtrList<globalPolyPatch>& toZones,   // to zones
                const List<Field<Type> >& fromFields,      // from fields
                List<Field<Type> >& toFields               // to fields
            ) const;

            //- Transfer/map/interpolate from one zone points to another zone
            //  points
            template<class Type>
            void transferPointsZonesToZones
            (
                const word& fromRegion,                    // from region name
                const word& toRegion,                      // to region name
                const PtrList<globalPolyPatch>& fromZones, // from zones
                const PtrList<globalPolyPatch>& toZones,   // to zones
                const List<Field<Type> >& fromFields,      // from fields
                List<Field<Type> >& toFields               // to fields
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

            //- Return fluid zones points displacements
            List<vectorField>& fluidZonesPointsDispls()
            {
                return fluidZonesPointsDispls_;
            }

            //- Return fluid zones ref. points displacements
            List<vectorField>& fluidZonesPointsDisplsRef()
            {
                return fluidZonesPointsDisplsRef_;
            }

            //- Return fluid zones previous points displacements
            List<vectorField>& fluidZonesPointsDisplsPrev()
            {
                return fluidZonesPointsDisplsPrev_;
            }

            //- Return solid zones points displacements
            List<vectorField>& solidZonesPointsDispls()
            {
                return solidZonesPointsDispls_;
            }

            //- Return solid zones ref. points displacements
            List<vectorField>& solidZonesPointsDisplsRef()
            {
                return solidZonesPointsDisplsRef_;
            }

            //- Return interfaces points displacements
            List<vectorField>& interfacesPointsDispls()
            {
                return interfacesPointsDispls_;
            }

            //- Return interfaces previous points displacements
            List<vectorField>& interfacesPointsDisplsPrev()
            {
                return interfacesPointsDisplsPrev_;
            }

            //- Return current fsi residuals
            List<vectorField>& residuals()
            {
                return residuals_;
            }

            //- Return previous fsi residuals
            List<vectorField>& residualsPrev()
            {
                return residualsPrev_;
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

            //- Return maximal residuals norm
            List<scalar>& maxResidualsNorm()
            {
                return maxResidualsNorm_;
            }

            //- Return const maximal residuals norm
            const List<scalar>& maxResidualsNorm() const
            {
                return maxResidualsNorm_;
            }

            //- Return interfaces displacements norm
            List<scalar>& maxIntsDisplsNorm()
            {
                return maxIntsDisplsNorm_;
            }

            //- Return const interfaces displacements norm
            const List<scalar>& maxIntsDisplsNorm() const
            {
                return maxIntsDisplsNorm_;
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
                List<vectorField>& fluidZonesPointsDispls
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
