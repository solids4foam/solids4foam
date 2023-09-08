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
    fluidModel

Description
    Virtual base class for fluid solvers.

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved.
    Zeljko Tukovic, FSB Zagreb.  All rights reserved.
    Philip Cardiff, UCD. All rights reserved.

SourceFiles
    fluidModel.C

\*---------------------------------------------------------------------------*/

#ifndef fluidModel_H
#define fluidModel_H

#include "physicsModel.H"
#include "fvMesh.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "IOdictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "pimpleControl.H"
#include "pisoControl.H"
#ifdef OPENFOAM_COM
        #include "gravityMeshObject.H"
#else
    #include "uniformDimensionedFields.H"
#endif
#include "globalPolyPatch.H"
#if FOAMEXTEND
    #include "oversetMesh.H"
#endif
#ifdef OPENFOAM_ORG
    #include "fvModels.H"
    #include "fvConstraints.H"
#elif OPENFOAM_COM
    #include "fvOptions.H"
#endif


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class fluidModel Declaration
\*---------------------------------------------------------------------------*/

class fluidModel
:
    public physicsModel,
    public IOdictionary
{
    // Private data

        //- Mesh
        autoPtr<dynamicFvMesh> meshPtr_;

        //- Fluid properties dictionary
        dictionary fluidProperties_;

        //- Piso control
        mutable autoPtr<pisoControl> pisoPtr_;

        //- Pimple control
        mutable autoPtr<pimpleControl> pimplePtr_;

        //- Wave properties needed to use waves2Foam
        const IOdictionary waveProperties_;

        //- Gravitional acceleration
#ifdef OPENFOAM_COM
        const meshObjects::gravity g_;
#else
        const uniformDimensionedVectorField g_;
#endif

        //- Velocity field header
        IOobject Uheader_;

        //- Pressure field header
        IOobject pheader_;

        //- Velocity field
        autoPtr<volVectorField> UPtr_;

        //- Pressure field
        autoPtr<volScalarField> pPtr_;

        //- Velocity gradient
        autoPtr<volTensorField> gradUPtr_;

        //- Pressure gradient
        autoPtr<volVectorField> gradpPtr_;

        //- Flux field
        autoPtr<surfaceScalarField> phiPtr_;

        //- Adjustable time-step switch
        const Switch adjustTimeStep_;

        //- For adjustable time-step, this is the maximum Courant number
        const scalar maxCo_;

        //- For adjustable time-step, this is the maximum time-stpe
        const scalar maxDeltaT_;

        //- Pressure minimum bound
        dimensionedScalar pMin_;

        //- Pressure maximum bound
        dimensionedScalar pMax_;

        //- Velocity bound
        dimensionedScalar UMax_;

        //- Small velocity for bounding
        dimensionedScalar smallU_;

        //- Cumulative continuity error
        scalar cumulativeContErr_;

#ifdef OPENFOAM_ORG
        //- fvModels
        fvModels& fvModels_;

        //- fvConstraints
        fvConstraints& fvConstraints_;
#elif OPENFOAM_COM
        //- fvOptions
        fv::options& fvOptions_;
#endif

        //- Flag to let the fluid model know that the fsi interface is in charge
        //  of calling fluid mesh.update()
        bool fsiMeshUpdate_;

        //- Flag to let the fluid model if the mesh changed when the FSI
        //  interface called mesh.updated()
        bool fsiMeshUpdateChanged_;

        //- Ptr list of global poly patches
        //  This is a copy of a patch which is reconstructed in full on each
        //  processor. This is convenient for interpolating between two patches
        //  in parallel, as used, for example, in the fluid-solid interaction
        //  procedure
        mutable PtrList<globalPolyPatch> globalPatchesPtrList_;


    // Private Member Functions

        //- Make the pimple control
        void makePisoControl() const;

        //- Make the pimple control
        void makePimpleControl() const;

        //- Read g field from disk, if it exists, otherwise return zero field
#ifdef OPENFOAM_COM
        meshObjects::gravity readG() const;
#else
        uniformDimensionedVectorField readG() const;
#endif

        //- Disallow default bitwise copy construct
        fluidModel(const fluidModel&);

        //- Disallow default bitwise assignment
        void operator=(const fluidModel&);


protected:

    // Protected member functions

        //- Return fluid properties dictionary
        const dictionary& fluidProperties() const
        {
            return fluidProperties_;
        }

        //- Update Robin boundary conditions, if found
        void updateRobinFsiInterface
        (
            const volScalarField& p,
            const volVectorField& U,
            surfaceScalarField& phi,
            surfaceScalarField& rAUf
        );

        //- Calculate the Courant number
        void CourantNo(scalar& CoNum, scalar& meanCoNum, scalar& velMag) const;

        //- Calculate the Courant number
        void CourantNo() const;

#if FOAMEXTEND
        //- Calculate the Courant number
        void oversetCourantNo
        (
            scalar& CoNum, scalar& meanCoNum, scalar& velMag
        ) const;

        //- Calculate the Courant number
        void oversetCourantNo() const;
#endif

        //- Continuity errors
        void continuityErrs();

#if FOAMEXTEND
        //- Continuity errors
        void oversetContinuityErrs();
#endif

        //- Bound the pressure and velocity fields
        void boundPU
        (
            volScalarField& p,
            volVectorField& U
        ) const;


public:

    //- Runtime type information
    TypeName("fluidModel");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            fluidModel,
            dictionary,
            (
                Time& runTime,
                const word& region
            ),
            (runTime, region)
        );


    // Constructors

        //- Construct from components
        fluidModel
        (
            const word& type,
            Time& runTime,
            const word& region = dynamicFvMesh::defaultRegion,
            const bool constructNull = false
        );


    // Selectors

        //- Select constructed from mesh
        static autoPtr<fluidModel> New
        (
            Time& runTime, const word& region = "region0"
        );


    // Destructor

        virtual ~fluidModel();


    // Member Functions

        // Access

            //- Return reference to the mesh
            dynamicFvMesh& mesh()
            {
                return meshPtr_();
            }

            //- Return const reference to the mesh
            const dynamicFvMesh& mesh() const
            {
                return meshPtr_();
            }

            //- Return time
            const Time& runTime() const
            {
                return meshPtr_().time();
            }

            //- Return non-const acess to the pimple control
            pisoControl& piso();

            //- Return non-const acess to the pimple control
            pimpleControl& pimple();

#if FOAMEXTEND
            //- Return non-const acess to the overset mesh
            const oversetMesh& osMesh() const;
#endif

            //- Return const gravity field
#ifdef OPENFOAM_COM
            virtual const meshObjects::gravity& g() const
#else
            virtual const uniformDimensionedVectorField& g() const
#endif
            {
                return g_;
            }

            //- Return const velocity field
            virtual const volVectorField& U() const
            {
                return UPtr_();
            }

            //- Return velocity field
            virtual volVectorField& U()
            {
                return UPtr_();
            }

            //- Return const pressure field
            virtual const volScalarField& p() const
            {
                return pPtr_();
            }

            //- Return pressure field
            virtual volScalarField& p()
            {
                return pPtr_();
            }

            //- Return const velocity gradient field
            virtual const volTensorField& gradU() const
            {
                return gradUPtr_();
            }

            //- Return velocity gradient field
            virtual volTensorField& gradU()
            {
                return gradUPtr_();
            }

            //- Return const pressure gradient field
            virtual const volVectorField& gradp() const
            {
                return gradpPtr_();
            }

            //- Return pressure gradient field
            virtual volVectorField& gradp()
            {
                return gradpPtr_();
            }

            //- Return const reference to flux field
            virtual const surfaceScalarField& phi() const
            {
                return phiPtr_();
            }

            //- Return reference to flux field
            virtual surfaceScalarField& phi()
            {
                return phiPtr_();
            }

#ifdef OPENFOAM_ORG
            //- Return fvModels
            const fvModels& models() const
            {
                return fvModels_;
            }

            //- Access fvModels
            fvModels& models()
            {
                return fvModels_;
            }

            //- Return  fvConstraints
            const fvConstraints& constraints() const
            {
                return fvConstraints_;
            }

            //- Access fvConstraints
            fvConstraints& constraints()
            {
                return fvConstraints_;
            }
#elif OPENFOAM_COM
            //- Return fvOptions
            const fv::options& options() const
            {
                return fvOptions_;
            }

            //- Access fvOptions
            fv::options& options()
            {
                return fvOptions_;
            }
#endif

            //- Non-const access to the fsiMeshUpdate flag
            virtual bool& fsiMeshUpdate()
            {
                return fsiMeshUpdate_;
            }

            //- Non-const access to the fsiMeshUpdateChanged flag
            virtual bool& fsiMeshUpdateChanged()
            {
                return fsiMeshUpdateChanged_;
            }

            //- Patch viscous force (N/m2)
            virtual tmp<vectorField> patchViscousForce
            (
                const label patchID
            ) const = 0;

            //- Patch pressure force (N/m2)
            virtual tmp<scalarField> patchPressureForce
            (
                const label patchID
            ) const = 0;

            //- Patch temperature
            virtual tmp<scalarField> patchTemperature(const label patchID) const
            {
                NotImplemented;
                return tmp<scalarField>(new scalarField());
            }

            //- Patch heat flux (W/m2)
            virtual tmp<scalarField> patchHeatFlux(const label patchID) const
            {
                NotImplemented;
                return tmp<scalarField>(new scalarField());
            }

            //- Patch heat transfer coefficient
            virtual tmp<scalarField> patchHeatTransferCoeff
            (
                const label patchID
            ) const
            {
                NotImplemented;
                return tmp<scalarField>(new scalarField());
            }

            //- Face zone viscous force (N/m2)
            virtual tmp<vectorField> faceZoneViscousForce
            (
                const label interfaceI
            ) const;

            //- Face zone pressure force (N/m2)
            virtual tmp<scalarField> faceZonePressureForce
            (
                const label interfaceI
            ) const;

            //- Face zone temperature
            virtual tmp<scalarField> faceZoneTemperature
            (
                const label interfaceI
            ) const;

            //- Face zone heat flux (W/m2)
            virtual tmp<scalarField> faceZoneHeatFlux
            (
                const label interfaceI
            ) const;

            //- Face zone heat transfer coefficient
            virtual tmp<scalarField> faceZoneHeatTransferCoeff
            (
                const label interfaceI
            ) const;

            //- This function will check that the U field was read from disk
            virtual void UisRequired();

            //- This function will check that the p field was read from disk
            virtual void pisRequired();

            //- Make the global poly patches
            void makeGlobalPatches(const wordList& patchNames) const;

            //- Return const-access to the global poly patches
            const PtrList<globalPolyPatch>& globalPatches() const;

            //- Clear global poly patches
            void clearGlobalPatches() const;

        // Edit

            //- Set temperature and heat flux at specified patch
            virtual void setTemperatureAndHeatFlux
            (
                const label interfaceI,
                const label patchID,
                const scalarField& faceZoneTemperature,
                const scalarField& faceZoneHeatFlux
            )
            {
                NotImplemented;
            }

            //- Set equivalent heat transfer coeff at specified patch
            virtual void setEqInterHeatTransferCoeff
            (
                const label interfaceI,
                const label patchID,
                const scalarField& faceZoneHeatTransferCoeff
            )
            {
                NotImplemented;
            }

            //- Update the size of the time-step
            virtual void setDeltaT(Time& runTime);

            //- Evolve the fluid model
            virtual bool evolve() = 0;

            //- Read fluidProperties dictionary
            virtual bool read();

            //- Write fields
            virtual void writeFields(const Time& runTime);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
