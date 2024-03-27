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
        mutable dictionary normalContactModelDict_;

        //- Const reference to mesh
        const fvMesh& mesh_;

        //- Pressure fields
        volVectorField pressureVolField_;

        //- Area in contact
        volScalarField areaInContactVolField_;

        //- Penalty factor
        mutable scalar penaltyFactor_;

        //- Penalty scale factor
        mutable scalar penaltyScale_;

        //- Under-relaxation factor
        const scalar relaxFac_;

        //- Current average penetration in the contact region
        //  This is zero if there is no contact
        scalar averagePenetration_;

        //- Current minimum (greatest) penetration in the contact region
        //  Remember that penetration is negative for points in contact
        scalar minPenetration_;

        //- Reference distance spring force calculation
        const scalar epsilon0_;

        //- Contact iteration number
        label contactIterNum_;


    // Private Member Functions

        //- Calculate penalty factor
        void calcPenaltyFactor() const;

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
            const scalarField& slavePointPenetration,
            const vectorField& slaveDU,
            const vectorField& masterDUInterpToSlave
        );

        //- Correct contact model
        virtual void correct
        (
            const vectorField& patchFaceNormals,
            const scalarField& faceVolPenetration,
            const scalarField& faceContactArea,
            const bool master
        );

        //- It is the start of a new time step e.g. update penalty factor
        virtual void newTimeStep() const
        {}

        //- Return slave pressure
        virtual const vectorField& slavePressure() const
        {
            return pressureVolField_.boundaryField()[slavePatchID()];
        }

        //- Return master pressure
        virtual const vectorField& masterPressure() const
        {
            return pressureVolField_.boundaryField()[masterPatchID()];
        }

        //- Return slave pressure
        virtual vectorField& slavePressure()
        {
#ifdef OPENFOAM_NOT_EXTEND
            return pressureVolField_.boundaryFieldRef()[slavePatchID()];
#else
            return pressureVolField_.boundaryField()[slavePatchID()];
#endif
        }

        //- Return master pressure
        virtual vectorField& masterPressure()
        {
#ifdef OPENFOAM_NOT_EXTEND
            return pressureVolField_.boundaryFieldRef()[masterPatchID()];
#else
            return pressureVolField_.boundaryField()[masterPatchID()];
#endif
        }

        //- Penalty factor
        virtual scalar penaltyFactor() const;

        //- Return area in contact of slave patch
        virtual const scalarField& slaveAreaInContact() const
        {
            return areaInContactVolField_.boundaryField()[slavePatchID()];
        };

        //- Return area in contact of master patch
        virtual const scalarField& masterAreaInContact() const
        {
            return areaInContactVolField_.boundaryField()[masterPatchID()];
        };


        //- Return area in contact of slave patch
        virtual scalarField& slaveAreaInContact()
        {
#ifdef OPENFOAM_NOT_EXTEND
            return areaInContactVolField_.boundaryFieldRef()[slavePatchID()];
#else
            return areaInContactVolField_.boundaryField()[slavePatchID()];
#endif
        };

        //- Return area in contact of master patch
        virtual scalarField& masterAreaInContact()
        {
#ifdef OPENFOAM_NOT_EXTEND
            return areaInContactVolField_.boundaryFieldRef()[masterPatchID()];
#else
            return areaInContactVolField_.boundaryField()[masterPatchID()];
#endif
        };

        //- Update penalty scale factors (for penalty laws)
        virtual scalar updatePenaltyScale(const scalar previousPenaltyScale);

        //- Perform any reauired field mapping e.g. after a topo change
        virtual void autoMap(const fvPatchFieldMapper& m);

        //- Write model dictionary
        virtual void writeDict(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
