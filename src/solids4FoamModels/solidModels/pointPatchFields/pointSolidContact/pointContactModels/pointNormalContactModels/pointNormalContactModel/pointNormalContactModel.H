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
    pointNormalContactModel

Description
    Model for the behaviour of the normal contact

SourceFiles
    pointNormalContactModel.C
    newPointNormalContactModel.C

Author
    Philip Cardiff, UCD. All rights reserved

\*---------------------------------------------------------------------------*/

#ifndef pointNormalContactModel_H
#define pointNormalContactModel_H

#include "IOdictionary.H"
#include "typeInfo.H"
#include "runTimeSelectionTables.H"
#include "volFields.H"
#include "tmp.H"
#include "autoPtr.H"
#include "standAlonePatch.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------* \
                         Class pointNormalContactModel Declaration
\*---------------------------------------------------------------------------*/

class pointNormalContactModel
{
    // Private data

        //- Name
        const word name_;

        //- Reference to mesh
        const fvMesh& mesh_;

        //- Master patch index
        const label masterPatchID_;

        //- Slave patch index
        const label slavePatchID_;

        //- Contact point gap
        //scalarField slaveContactPointGap_;

    // Private Member Functions

        //- Disallow default bitwise assignment
        void operator=(const pointNormalContactModel&);


protected:

        //- Return reference to mesh
        const fvMesh& mesh() const
        {
            return mesh_;
        }


public:

    //- Runtime type information
    TypeName("pointNormalContactModel");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            pointNormalContactModel,
            dictionary,
            (
                const word name,
                const fvMesh& mesh, //const fvPatch& patch,
                const dictionary& dict,
                const label masterPatchID,
                const label slavePatchID,
                const standAlonePatch& masterFaceZonePatch,
                const standAlonePatch& slaveFaceZonePatch
            ),
            (
                name,
                mesh, //patch,
                dict,
                masterPatchID,
                slavePatchID,
                masterFaceZonePatch,
                slaveFaceZonePatch
            )
        );


    // Selectors

        //- Return a reference to the selected mechanical model
        static autoPtr<pointNormalContactModel> New
        (
            const word& name,
            const fvMesh& mesh, //const fvPatch& patch,
            const dictionary& dict,
            const label masterPatchID,
            const label slavePatchID,
            const standAlonePatch& masterFaceZonePatch,
            const standAlonePatch& slaveFaceZonePatch
        );


    // Constructors

        //- Construct from dictionary
        pointNormalContactModel
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
        pointNormalContactModel(const pointNormalContactModel&);

        //- Construct and return a clone
        virtual autoPtr<pointNormalContactModel> clone() const = 0;


    //- Destructor
    virtual ~pointNormalContactModel()
    {}


    // Member Functions

        //- Return name
        const word& name() const
        {
            return name_;
        }

        //- Correct contact model
        virtual void correct
        (
            const vectorField& slavePatchPointNormals,
            const scalarField& slavePointPenetration //,
            // const vectorField& slaveDU,
            // const vectorField& masterDUInterpToSlave
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

        //- Return point gap field
        // virtual scalarField& slaveContactPointGap()
        // {
        //     return slaveContactPointGap_;
        // }

        // //- Return area in contact of slave patch
        // virtual const scalarField& areaInContact() const = 0;

        //- Update penalty scale factors (for penalty laws)
        virtual scalar updatePenaltyScale(const scalar previousPenaltyScale)
        {
            return previousPenaltyScale;
        }

        //- Perform any reauired field mapping e.g. after a topo change
        // virtual void autoMap(const fvPatchFieldMapper& m)
        // {}

        //- Write model dictionary
        virtual void writeDict(Ostream& os) const {};
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
