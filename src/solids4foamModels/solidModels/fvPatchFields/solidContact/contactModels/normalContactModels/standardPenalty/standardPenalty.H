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
    standardPenalty

Description
    Standard penalty method for normal contact law.

    Using a fixed penalty factor, force is applied to the contact faces based on
    the gap.

    The penalty method is effectively like attaching springs to any face that
    comes into contact.

SourceFiles
    standardPenalty.C

Author
    Philip Cardiff, UCD. All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef standardPenalty_H
#define standardPenalty_H

#include "normalContactModel.H"
#include "typeInfo.H"
#include "runTimeSelectionTables.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class standardPenalty Declaration
\*---------------------------------------------------------------------------*/

class standardPenalty
:
    public normalContactModel
{
    // Private data

        //- Contact law dictionary
        dictionary normalContactModelDict_;

        //- Const reference to mesh
        const fvMesh& mesh_;

        //- Write debug file switch
        const Switch writeDebugFile_;

        //- Slave pressure
        vectorField slavePressure_;

        //- Global and local slave point penetrations
        scalarField globalSlavePointPenetration_;

        //- Local slave point penetrations
        scalarField slavePointPenetration_;

        //- Area in contact
        scalarField areaInContact_;

        //- Penalty factor
        mutable scalar* penaltyFactorPtr_;

        //- Penalty scale factor
        const scalar penaltyScale_;

        //- Under-relaxation factor
        const scalar relaxFac_;

        //- Current slave point contact pressure
        scalarField totalSlavePointPressure_;

        //- Contact iteration number
        label contactIterNum_;

        //- Frequency for printing info to the contact debug file
        const int infoFreq_;

        //- Contact info file
        OFstream* contactFilePtr_;


    // Private Member Functions

        //- Calculate penalty factor
        void calcPenaltyFactor() const;

        //- Disallow copy construct
        //standardPenalty(const standardPenalty&);

        //- Disallow default bitwise assignment
        void operator=(const standardPenalty&);



public:

    //- Runtime type information
    TypeName("standardPenalty");


    // Constructors

        //- Construct from dictionary
        standardPenalty
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
        standardPenalty(const standardPenalty&);

        //- Construct and return a clone
        virtual autoPtr<normalContactModel> clone() const
        {
            return autoPtr<normalContactModel>(new standardPenalty(*this));
        }


    // Destructor

        virtual ~standardPenalty()
        {}


    // Member Functions


        //- Correct contact model
        virtual void correct
        (
            const vectorField& slavePatchFaceNormals,
            const extendedGgiStandAlonePatchInterpolation& zoneToZone
        );

        //- Return slave pressure
        virtual const vectorField& slavePressure() const
        {
            return slavePressure_;
        }

        //- Return slave pressure
        virtual vectorField& slavePressure()
        {
            return slavePressure_;
        }

        //- Penalty factor
        virtual scalar penaltyFactor() const;

        //- Return area in contact of slave patch
        virtual const scalarField& areaInContact() const
        {
            return areaInContact_;
        };

        //- Write model dictionary
        virtual void writeDict(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
