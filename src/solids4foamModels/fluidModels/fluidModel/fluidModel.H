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

#include "fvMesh.H"
#include "IOdictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class fluidModel Declaration
\*---------------------------------------------------------------------------*/

class fluidModel
:
    public IOdictionary
{
    // Private data

        //- Reference to mesh
        const fvMesh& mesh_;

        //- Fluid properties dictionary
        dictionary fluidProperties_;

        //- Global face zones
        mutable labelList* globalFaceZonesPtr_;

        //- Global to local face zone point maping
        mutable labelListList* globalToLocalFaceZonePointMapPtr_;

    // Private Member Functions

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
                const fvMesh& mesh
            ),
            (mesh)
        );


    // Constructors

        //- Construct from components
        fluidModel
        (
            const word& type,
            const fvMesh& mesh
        );


    // Selectors

        //- Select constructed from mesh
        static autoPtr<fluidModel> New(const fvMesh& mesh);


    // Destructor

        virtual ~fluidModel();


    // Member Functions

        // Access

            //- Return mesh
            const fvMesh& mesh() const
            {
                return mesh_;
            }

            //- Return time
            const Time& runTime() const
            {
                return mesh_.time();
            }

            //- Return velocity field
            //virtual const volVectorField& U() const = 0;

            //- Return pressure field
            //virtual const volScalarField& p() const = 0;

            //- Patch viscous force (N/m2)
            // virtual tmp<vectorField> patchViscousForce
            // (
            //     const label patchID
            // ) const = 0;

            // //- Patch pressure force (N/m2)
             virtual tmp<scalarField> patchPressureForce
             (
                 const label patchID
             ) const = 0;

            // //- Patch viscous force (N/m2)
             virtual tmp<vectorField> faceZoneViscousForce
             (
                 const label zoneID,
                 const label patchID
             ) const = 0;

            // //- Patch pressure force (N/m2)
             virtual tmp<scalarField> faceZonePressureForce
             (
                 const label zoneID,
                 const label patchID
             ) const = 0;

            // //- Face zone effective dynamic viscosity
             virtual tmp<scalarField> faceZoneMuEff
             (
                 const label zoneID,
                 const label patchID
             ) const = 0;

            // //- Get global face zones list
             const labelList& globalFaceZones() const;

            // //- Get global to local face zone point addressing
             const labelListList& globalToLocalFaceZonePointMap() const;

            //- Return point displacement field
            virtual const pointVectorField& pointD() const = 0;

            //- Face zone primitive patch in current configuration
             virtual tmp<vectorField>
             currentFaceZonePoints(const label zoneID) const = 0;

        // Edit

            //- Evolve the fluid model
            virtual void evolve() = 0;

            //- Read fluidProperties dictionary
            virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
