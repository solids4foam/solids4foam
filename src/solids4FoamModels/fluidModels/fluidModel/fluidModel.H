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
    fluidModel

Description
    Virtual base class for fluid solvers.

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved.
    Zeljko Tukovic, FSB Zagreb.  All rights reserved.
    Philip Cardiff, UCD. All rights reserved.

SourceFiles
    fluidModel.C
    newFluidModel.C

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

        //- Pimple control
        mutable autoPtr<pimpleControl> pimplePtr_;

        //- Global face zones
        mutable labelList* globalFaceZonesPtr_;

        //- Global to local face zone point maping
        mutable labelListList* globalToLocalFaceZonePointMapPtr_;

        //- Velocity field
        volVectorField U_;

        //- Pressure field
        volScalarField p_;

        //- Velocity gradient
        volTensorField gradU_;

        //- Pressure gradient
        volVectorField gradp_;

        //- Flux field
        surfaceScalarField phi_;

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

        //- Flag to let the fluid model know that the fsi interface is in charge
        //  of calling fluid mesh.update()
        bool fsiMeshUpdate_;

        //- Flag to let the fluid model if the mesh changed when the FSI
        //  interface called mesh.updated()
        bool fsiMeshUpdateChanged_;

    // Private Member Functions

        //- Make the pimple control
        void makePimpleControl() const;

        //- Find global face zones
        void calcGlobalFaceZones() const;

        // //- Calculate global to local point addressing
        void calcGlobalToLocalFaceZonePointMap() const;

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

        //- Continuity errors
        void continuityErrs();

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
            const word& region = dynamicFvMesh::defaultRegion
        );


    // Selectors

        //- Select constructed from mesh
        static autoPtr<fluidModel> New(Time& runTime, const word& region);


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
            pimpleControl& pimple();

            //- Return const velocity field
            virtual const volVectorField& U() const
            {
                return U_;
            }

            //- Return velocity field
            virtual volVectorField& U()
            {
                return U_;
            }

            //- Return const pressure field
            virtual const volScalarField& p() const
            {
                return p_;
            }

            //- Return pressure field
            virtual volScalarField& p()
            {
                return p_;
            }

            //- Return const velocity gradient field
            virtual const volTensorField& gradU() const
            {
                return gradU_;
            }

            //- Return velocity gradient field
            virtual volTensorField& gradU()
            {
                return gradU_;
            }

            //- Return const pressure gradient field
            virtual const volVectorField& gradp() const
            {
                return gradp_;
            }

            //- Return pressure gradient field
            virtual volVectorField& gradp()
            {
                return gradp_;
            }

            //- Return const reference to flux field
            virtual const surfaceScalarField& phi() const
            {
                return phi_;
            }

            //- Return reference to flux field
            virtual surfaceScalarField& phi()
            {
                return phi_;
            }

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

            //- Patch viscous force (N/m2)
            virtual tmp<vectorField> faceZoneViscousForce
            (
                const label zoneID,
                const label patchID
            ) const;

            //- Patch pressure force (N/m2)
            virtual tmp<scalarField> faceZonePressureForce
            (
                const label zoneID,
                const label patchID
            ) const;

            // //- Face zone effective dynamic viscosity
            virtual tmp<scalarField> faceZoneMuEff
            (
                const label zoneID,
                const label patchID
            ) const = 0;

            //- Get global face zones list
            const labelList& globalFaceZones() const;

            //- Get global to local face zone point addressing
            const labelListList& globalToLocalFaceZonePointMap() const;

        // Edit

            //- Update the size of the time-step
            virtual void setDeltaT(Time& runTime);

            //- Evolve the fluid model
            virtual bool evolve() = 0;

            //- Read fluidProperties dictionary
            virtual bool read();

            //- Write fields
            virtual void writeFields(const Time& runTime)
            {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
