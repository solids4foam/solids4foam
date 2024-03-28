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
    frictionless

Description
    frictionless contact law

SourceFiles
    frictionless.C

Author
    Philip Cardiff, UCD. All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef frictionless_H
#define frictionless_H

#include "frictionContactModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class frictionless Declaration
\*---------------------------------------------------------------------------*/

class frictionless
:
    public frictionContactModel
{
    // Private data

       //- Slave friction traction - always zero for frictionless
       vectorField slaveTraction_;

       //- Master friction traction - always zero for frictionless
       vectorField masterTraction_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        //frictionless(const frictionless&);

        //- Disallow default bitwise assignment
        void operator=(const frictionless&);

public:

    //- Runtime type information
    TypeName("frictionless");

    // Static data members


    // Constructors

        //- Construct from dictionary
        frictionless
        (
            const word& name,
            const fvPatch& patch,
            const dictionary& dict,
            const label masterPatchID,
            const label slavePatchID
        );

        //- Construct as copy
        frictionless(const frictionless&);

        //- Construct and return a clone
        virtual autoPtr<frictionContactModel> clone() const
        {
            return autoPtr<frictionContactModel>
            (
                new frictionless(*this)
            );
        }


    //- Destructor
    virtual ~frictionless();


    // Member Functions

        //- Correct the friction model
        virtual void correct
        (
            const vectorField& slavePressure,
            const vectorField& slaveFaceNormals,
            const scalarField& slavePatchAreaInContact,
            const vectorField& slaveDU,
            const vectorField& masterDUInterpToSlave
        )
        {}

        //- Correct the friction model
        virtual void correct
        (
            const vectorField& patchPressure,
            const vectorField& patchFaceNormals,
            const scalarField& patchAreaInContact,
            const vectorField& DD,
            const vectorField& shadowDDInterpToPatch,
            const bool master
        )
        {}

        //- Return slave friction traction
        virtual const vectorField& slaveTraction() const
        {
            return slaveTraction_;
        }

        //- Return master friction traction
        virtual const vectorField& masterTraction() const
        {
            return masterTraction_;
        }

        //- Return slip field
        virtual const vectorField& slipOnSlave() const
        {
            // slaveTraction is always zero
            return slaveTraction_;
        }

        //- Return slip field
        virtual const vectorField& slipOnMaster() const
        {
            // masterTraction is always zero
            return masterTraction_;
        }
        //- Perform any reauired field mapping e.g. after a topo change
        virtual void autoMap(const fvPatchFieldMapper& m);

        //- Write law dictionary
        virtual void writeDict(Ostream& os) const
        {};
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
