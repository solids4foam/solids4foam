/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2004-2007 Hrvoje Jasak
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
    Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

Class
    normalContactModel

Description
    Model for the behaviour of the normal contact

SourceFiles
    normalContactModel.C
    newNormalContactModel.C

Author
    Philip Cardiff, UCD. All rights reserved

\*---------------------------------------------------------------------------*/

#ifndef normalContactModel_H
#define normalContactModel_H

#include "IOdictionary.H"
#include "typeInfo.H"
#include "runTimeSelectionTables.H"
#include "volFields.H"
#include "tmp.H"
#include "autoPtr.H"
#include "extendedGgiInterpolation.H"
#include "standAlonePatch.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------* \
                         Class normalContactModel Declaration
\*---------------------------------------------------------------------------*/

class normalContactModel
{
    // Private data

        //- Name
        const word name_;

        //- Reference to primitivePatch
        const fvPatch& patch_;

        //- Master patch index
        const label masterPatchID_;

        //- Slave patch index
        const label slavePatchID_;

        //- Master face zone index
        const label masterFaceZoneID_;

        //- Slave face zone index
        const label slaveFaceZoneID_;

        //- Contact point gap
        scalarField slaveContactPointGap_;

    // Private Member Functions

        //- Disallow copy construct
        //normalContactModel(const normalContactModel&);

        //- Disallow default bitwise assignment
        void operator=(const normalContactModel&);


protected:

        //- Return reference to mesh
        const fvMesh& mesh() const
        {
            return patch_.boundaryMesh().mesh();
        }


public:

    //- Runtime type information
    TypeName("normalContactModel");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            normalContactModel,
            dictionary,
            (
                const word name,
                const fvPatch& patch,
                const dictionary& dict,
                const label masterPatchID,
                const label slavePatchID,
                const label masterFaceZoneID,
                const label slaveFaceZoneID,
                const standAlonePatch& masterFaceZonePatch,
                const standAlonePatch& slaveFaceZonePatch
            ),
            (
                name,
                patch,
                dict,
                masterPatchID,
                slavePatchID,
                masterFaceZoneID,
                slaveFaceZoneID,
                masterFaceZonePatch,
                slaveFaceZonePatch
            )
        );


    // Selectors

        //- Return a reference to the selected mechanical model
        static autoPtr<normalContactModel> New
        (
            const word& name,
            const fvPatch& patch,
            const dictionary& dict,
            const label masterPatchID,
            const label slavePatchID,
            const label masterFaceZoneID,
            const label slaveFaceZoneID,
            const standAlonePatch& masterFaceZonePatch,
            const standAlonePatch& slaveFaceZonePatch
        );


    // Constructors

        //- Construct from dictionary
        normalContactModel
        (
            const word& name,
            const fvPatch& patch,
            const dictionary& dict,
            const label masterPatchID,
            const label slavePatchID,
            const label masterFaceZoneID,
            const label slaveFaceZoneID,
            const standAlonePatch& masterFaceZonePatch,
            const standAlonePatch& slaveFaceZonePatch
        );

        //- Construct as copy
        normalContactModel(const normalContactModel&);

        //- Construct and return a clone
        virtual autoPtr<normalContactModel> clone() const = 0;


    // Destructor

        virtual ~normalContactModel()
        {}


    // Member Functions

        //- Return name
        const word& name() const
        {
            return name_;
        }

        //- Return reference to patch
        const fvPatch& patch() const
        {
            return patch_;
        }

        //- Correct contact model
        virtual void correct
        (
            const vectorField& slavePatchFaceNormals,
            const extendedGgiStandAlonePatchInterpolation& zoneToZone
        ) = 0;

        //- It is the start of a new time step
        virtual void newTimeStep() const
        {}

        //- Return slave pressure
        virtual const vectorField& slavePressure() const = 0;

        //- Return master patch ID
        virtual label masterPatchID() const
        {
            return masterPatchID_;
        }

        //- Return master patch ID
        virtual label slavePatchID() const
        {
            return slavePatchID_;
        }

        //- Return master face zone ID
        virtual label masterFaceZoneID() const
        {
            return masterFaceZoneID_;
        }

        //- Return master face zone ID
        virtual label slaveFaceZoneID() const
        {
            return slaveFaceZoneID_;
        }

        //- Return point gap field
        virtual scalarField& slaveContactPointGap()
        {
            return slaveContactPointGap_;
        }

        //- Return area in contact of slave patch
        virtual const scalarField& areaInContact() const = 0;

        //- Write model dictionary
        virtual void writeDict(Ostream& os) const {};
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
