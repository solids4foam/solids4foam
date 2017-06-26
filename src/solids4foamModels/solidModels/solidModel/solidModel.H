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
    solidModel

Description
    Virtual base class for solid mechanics models

Author
    Zeljko Tukovic, FSB Zagreb.  All rights reserved.
    Philip Cardiff, UCD. All rights reserved.

SourceFiles
    solidModel.C
    newSolidModel.C

\*---------------------------------------------------------------------------*/

#ifndef solidModel_H
#define solidModel_H

#include "dynamicFvMesh.H"
#include "IOdictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "mechanicalModel.H"
#include "nonLinearGeometry.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class solidModel Declaration
\*---------------------------------------------------------------------------*/

class solidModel
:
    public IOdictionary
{
    // Private data

        //- Non-const reference to mesh
        dynamicFvMesh& mesh_;

        //- Solid properties dictionary
        dictionary solidProperties_;

        //- Mechanical model
        mutable autoPtr<mechanicalModel> mechanicalPtr_;

        //- Global face zones
        mutable labelList* globalFaceZonesPtr_;

        //- Global to local face zone point mapping
        mutable labelListList* globalToLocalFaceZonePointMapPtr_;

        //- Switch to enable/disable enforce linear to help convergence
        Switch enforceLinear_;

    // Private Member Functions

        //- Find global face zones
        void calcGlobalFaceZones() const;

        //- Find global to local point addressing
        void calcGlobalToLocalFaceZonePointMap() const;

        //- Make the mechanical model
        void makeMechanicalModel() const;

        //- Disallow default bitwise copy construct
        solidModel(const solidModel&);

        //- Disallow default bitwise assignment
        void operator=(const solidModel&);


protected:

    // Protected member functions

        //- Return nonlinear geometry enumerator
        virtual nonLinearGeometry::nonLinearType nonLinGeom() const = 0;

        //- Return solid properties dictionary
        dictionary& solidProperties()
        {
            return solidProperties_;
        }

        //- Return mechanical model
        mechanicalModel& mechanical();

public:

    //- Runtime type information
    TypeName("solidModel");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            solidModel,
            dictionary,
            (
                dynamicFvMesh& mesh
            ),
            (mesh)
        );


    // Constructors

        //- Construct from components
        solidModel
        (
            const word& type,
            dynamicFvMesh& mesh
        );


    // Selectors

        //- Select constructed from mesh
        static autoPtr<solidModel> New(dynamicFvMesh& mesh);


    // Destructor

        virtual ~solidModel();


    // Member Functions

        // Access

            //- Return mesh
            const dynamicFvMesh& mesh() const
            {
                return mesh_;
            }

            //- Return non-const mesh
            dynamicFvMesh& mesh()
            {
                return mesh_;
            }

            //- Return time
            const Time& runTime() const
            {
                return mesh_.time();
            }

            //- Return solid properties dictionary
            const dictionary& solidProperties() const
            {
                return solidProperties_;
            }

            //- Return mechanical model
            const mechanicalModel& mechanical() const;

            //- Face zone point displacement increment
            virtual tmp<vectorField> faceZonePointDisplacementIncrement
            (
                 const label zoneID
            ) const = 0;

            //- Face zone velocity gradient
            virtual tmp<tensorField> faceZoneSurfaceGradientOfVelocity
            (
                const label zoneID,
                const label patchID
            ) const = 0;

            //- Face zone current normal
            virtual tmp<vectorField> faceZoneNormal
            (
                const label zoneID,
                const label patchID
            ) const = 0;

            //- Face zone primitive patch in current configuration
            virtual tmp<vectorField> currentFaceZonePoints
            (
                const label zoneID
            ) const = 0;

            //- Get global face zones list
            const labelList& globalFaceZones() const;

            //- Get global to local face zone point addressing
            const labelListList& globalToLocalFaceZonePointMap() const;

            //- Return enforceLinear switch
            const Switch& enforceLinear() const
            {
                return enforceLinear_;
            }

            //- Return enforceLinear switch
            Switch& enforceLinear()
            {
                return enforceLinear_;
            }

        // Edit

            //- Update the points on global face zones when moving the mesh
            void updateGlobalFaceZoneNewPoints
            (
                const pointField& pointDDI,
                pointField& newPoints
            );

            //- Set traction at specified patch
            virtual void setTraction
            (
                const label patchID,
                const vectorField& traction
            ) = 0;

            //- Set traction at specified patch
            virtual void setPressure
            (
                const label patchID,
                const scalarField& pressure
            ) = 0;

            //- Set traction at specified patch
            virtual void setTraction
            (
                const label patchID,
                const label zoneID,
                const vectorField& faceZoneTraction
            );

            //- Set traction at specified patch
            virtual void setPressure
            (
                const label patchID,
                const label zoneID,
                const scalarField& faceZonePressure
            );

            //- Evolve the solid model
            virtual bool evolve() = 0;

            //- Traction boundary surface normal gradient
            //  Given the user-specified boundary patch traction, this function
            //  should return the surface normal gradient
            //  This function is used by traction-type boundary conditions to
            //  set the boundary gradient
            virtual tmp<Foam::vectorField> tractionBoundarySnGrad
            (
                const vectorField& traction,
                const scalarField& pressure,
                const fvPatch& patch
            ) const = 0;

            //- Update total accumulated fields
            virtual void updateTotalFields() = 0;

            //- Write fields
            virtual void writeFields(const Time& runTime);

            //- Return the desired new time-step size
            virtual scalar newDeltaT();

            //- Move the mesh to the deformed configuration
            virtual void moveMesh
            (
                const pointField& oldPoints,
                const volVectorField& DD,
                pointVectorField& pointDD
            );

            //- End of the program
            virtual void end()
            {}

            //- Read solidProperties dictionary
            virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
