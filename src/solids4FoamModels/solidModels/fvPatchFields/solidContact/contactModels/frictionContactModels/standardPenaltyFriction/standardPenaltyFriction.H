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
    standardPenaltyFriction

Description
    Friction contact law using the penslty method to enforce the friction
    constraints.

SourceFiles
    standardPenaltyFriction.C

Author
    Philip Cardiff, UCD. All rights reserved.
    Matthias Schnaithmann, Universitat Stuttgart. All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef standardPenaltyFriction_H
#define standardPenaltyFriction_H

#include "frictionContactModel.H"
#include "frictionLaw.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class standardPenaltyFriction Declaration
\*---------------------------------------------------------------------------*/

class standardPenaltyFriction
:
    public frictionContactModel
{
    // Private data

         //- The contact model dictionary
        dictionary frictionContactModelDict_;

         //- Friction law
        autoPtr<frictionLaw> frictionLawPtr_;

         //- Const reference to the mesh
        const fvMesh& mesh_;

        //- Friction traction on slave or on slave and master
        volVectorField tractionVolField_;

        //- Previous slave traction used for under-relaxation
        //vectorField prevSlaveTraction_;

        //- Slave slip field
        vectorField slipOnSlave_;

        //- Master slip field
        vectorField slipOnMaster_;

        //- Frictional penalty factor and scaling factor
        scalar frictionPenaltyFactor_;

        //- Frictional penalty scaling factor
        scalar frictionPenaltyScale_;

        //- Under-relaxation factor
        scalar relaxFac_;

        //- Contact iteration number
        label contactIterNum_;


    // Private Member Functions

        //- Return penalty factor
        virtual scalar frictionPenaltyFactor();

        //- Disallow default bitwise copy construct
        //standardPenaltyFriction(const standardPenaltyFriction&);

        //- Disallow default bitwise assignment
        void operator=(const standardPenaltyFriction&);

public:

    //- Runtime type information
    TypeName("standardPenalty");


    // Constructors

        //- Construct from dictionary
        standardPenaltyFriction
        (
            const word& name,
            const fvPatch& patch,
            const dictionary& dict,
            const label masterPatchID,
            const label slavePatchID
        );

        //- Construct as copy
        standardPenaltyFriction(const standardPenaltyFriction&);

        //- Construct and return a clone
        virtual autoPtr<frictionContactModel> clone() const
        {
            return
                autoPtr<frictionContactModel>
                (
                    new standardPenaltyFriction(*this)
                );
        }


    //- Destructor
    virtual ~standardPenaltyFriction();


    // Member Functions

        //- Correct the friction model
        virtual void correct
        (
            const vectorField& slavePressure,
            const vectorField& slaveFaceNormals,
            const scalarField& slavePatchAreaInContact,
            const vectorField& slaveDU,
            const vectorField& masterDUInterpToSlave
        );

        //- Correct the friction model
        virtual void correct
        (
            const vectorField& patchPressure,
            const vectorField& patchFaceNormals,
            const scalarField& patchAreaInContact,
            const vectorField& DD,
            const vectorField& shadowDDInterpToPatch,
            const bool master
        );

        //- Return const reference to slave friction traction
        virtual const vectorField& slaveTraction() const
        {
            return tractionVolField_.boundaryField()[slavePatchID()];
        }

        //- Return const reference to master friction traction
        virtual const vectorField& masterTraction() const
        {
            return tractionVolField_.boundaryField()[masterPatchID()];
        }

        //- Return non-const reference to slave friction traction
        virtual vectorField& slaveTraction()
        {
#ifdef OPENFOAM_NOT_EXTEND
            return tractionVolField_.boundaryFieldRef()[slavePatchID()];
#else
            return tractionVolField_.boundaryField()[slavePatchID()];
#endif
        }

        //- Return non-const reference to master friction traction
        virtual vectorField& masterTraction()
        {
#ifdef OPENFOAM_NOT_EXTEND
            return tractionVolField_.boundaryFieldRef()[masterPatchID()];
#else
            return tractionVolField_.boundaryField()[masterPatchID()];
#endif
        }

        //- calculate penalty factor
        void calcFrictionPenaltyFactor();

        //- Return slip field on slave
        virtual const vectorField& slipOnSlave() const
        {
            // Note: slip is updated inside the correct function
            return slipOnSlave_;
        }

        //- Return slip field on slave
        virtual const vectorField& slipOnMaster() const
        {
            // Note: slip is updated inside the correct function
            return slipOnMaster_;
        }

        //- Return slip field
        virtual vectorField& slip(bool master)
        {
            if (master)
            {
                return slipOnMaster_;
            }
            else
            {
                return slipOnSlave_;
            }
        }


        //- Perform any reauired field mapping e.g. after a topo change
        virtual void autoMap(const fvPatchFieldMapper& m);

        //- Write law dictionary
        virtual void writeDict(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
