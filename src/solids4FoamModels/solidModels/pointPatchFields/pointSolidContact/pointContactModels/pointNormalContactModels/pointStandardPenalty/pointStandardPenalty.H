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
    pointStandardPenalty

Description
    Point version of StandardPenalty method for normal contact law.

SourceFiles
    pointStandardPenalty.C

Author
    Philip Cardiff, UCD. All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef pointStandardPenalty_H
#define pointStandardPenalty_H

#include "pointNormalContactModel.H"
#include "typeInfo.H"
#include "runTimeSelectionTables.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class pointStandardPenalty Declaration
\*---------------------------------------------------------------------------*/

class pointStandardPenalty
:
    public pointNormalContactModel
{
    // Private data

        //- Contact law dictionary
        mutable dictionary pointNormalContactModelDict_;

        //- Slave point pressure
        vectorField slavePressure_;
        //volVectorField slavePressureVolField_;

        //- Area in contact
        // scalarField areaInContact_;
        // volScalarField areaInContactVolField_;

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
        void operator=(const pointStandardPenalty&);


public:

    //- Runtime type information
    TypeName("standardPenalty");


    // Constructors

        //- Construct from dictionary
        pointStandardPenalty
        (
            const word& name,
            const fvMesh& mesh, //const fvPatch& patch,
            const dictionary& dict,
            const label masterPatchID,
            const label slavePatchID,
            const standAlonePatch& masterFaceZonePatch,
            const standAlonePatch& slaveFaceZonePatch
        );

        //- Construct as copy
        pointStandardPenalty(const pointStandardPenalty&);

        //- Construct and return a clone
        virtual autoPtr<pointNormalContactModel> clone() const
        {
            return autoPtr<pointNormalContactModel>(new pointStandardPenalty(*this));
        }


    // Destructor

        virtual ~pointStandardPenalty()
        {}


    // Member Functions


        //- Correct contact model
        virtual void correct
        (
            const vectorField& slavePatchPointNormals,
            const scalarField& slavePointPenetration //,
            // const vectorField& slaveDU,
            // const vectorField& masterDUInterpToSlave
        );

        //- It is the start of a new time step e.g. update penalty factor
        virtual void newTimeStep() const
        {}

        //- Return slave pressure
        virtual const vectorField& slavePressure() const
        {
            return slavePressure_;
            // return slavePressureVolField_.boundaryField()[slavePatchID()];
        }

        //- Return slave pressure
        virtual vectorField& slavePressure()
        {
            return slavePressure_;
// #ifdef OPENFOAM_NOT_EXTEND
//             return slavePressureVolField_.boundaryFieldRef()[slavePatchID()];
// #else
//             return slavePressureVolField_.boundaryField()[slavePatchID()];
// #endif
        }

        //- Penalty factor
        virtual scalar penaltyFactor() const;

        //- Return area in contact of slave patch
        // virtual const scalarField& areaInContact() const
        // {
        //     return areaInContact_;
        //     // return areaInContactVolField_.boundaryField()[slavePatchID()];
        // };

        //- Return area in contact of slave patch
//         virtual scalarField& areaInContact()
//         {
//             return areaInContact_;
// // #ifdef OPENFOAM_NOT_EXTEND
// //             return areaInContactVolField_.boundaryFieldRef()[slavePatchID()];
// // #else
// //             return areaInContactVolField_.boundaryField()[slavePatchID()];
// // #endif
//         };

        //- Update penalty scale factors (for penalty laws)
        virtual scalar updatePenaltyScale(const scalar previousPenaltyScale);

        //- Perform any reauired field mapping e.g. after a topo change
        // virtual void autoMap(const fvPatchFieldMapper& m);

        //- Write model dictionary
        virtual void writeDict(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
