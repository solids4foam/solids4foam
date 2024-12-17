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
#include "interfaceToInterfaceMapping.H"


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
        dictionary& fsiProperties_;

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

        //- Mapping objects for passing fields between the pairs of solid and
        //  fluid interfaces
        //  e.g. this could be RBF or GGI or a direct map or other
        mutable PtrList<interfaceToInterfaceMapping> interfaceToInterfaceList_;

        //- Outer correction loop stoping tolerance
        const scalar outerCorrTolerance_;

        //- Maximal number of outer correctors
        const label nOuterCorr_;

        //- Optional additional mesh correction step after FSI loop
        const Switch additionalMeshCorrection_;

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

        //- Incremental or absolute residual calculation approach
        const Switch incrementalResiduals_;

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

    // Private Member Functions

        //- Calculate interface-to-interface list
        void calcInterfaceToInterfaceList() const;

        //- Accumulated fluid interfaces displacements
        void calcAccumulatedFluidInterfacesDisplacements() const;

        //- Return accumulated interfaces displacements
        PtrList<vectorIOField>& accumulatedFluidInterfacesDisplacements();

        //- Disallow default bitwise copy construct
        fluidSolidInterface(const fluidSolidInterface&);

        //- Disallow default bitwise assignment
        void operator=(const fluidSolidInterface&);

protected:
    // Protected Member Functions

        //- Update the coupled switch based on the couplingStartTime
        bool updateCoupled();

        //- Return sub-dictionary
        dictionary& fsiProperties()
        {
            return fsiProperties_;
        }

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

            //- Return const access to the interface-to-interface list
            const PtrList<interfaceToInterfaceMapping>&
                interfaceToInterfaceList() const;

            //- Return const access to the additional mesh correction switch
            const Switch& additionalMeshCorrection() const
            {
                return additionalMeshCorrection_;
            }


        // Edit

            //- Return the current force on the interface given the zone and
            //  zone traction field
            vector totalForceOnInterface
            (
                const standAlonePatch& zone, const vectorField& zoneTraction
            ) const;

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

            //- Update separately interface viscous force and pressure
            virtual void updateViscousForceAndPressure();

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

            //- This function is called at the end of the simulation
            virtual void end();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
