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

Class
    Foam::mechanicalConstitutiveLawStateSpec

Description
    What a mechanical constitutive law declares about the state it keeps.

    A law states the name, type, role and default of every state field it
    uses, and the manager does the rest: it allocates the fields at the right
    size for the topology, applies the defaults, and for a prescribed field
    looks for a volField of the same name on the mesh and broadcasts it to the
    law's integration points.

    This exists so that a law never touches an IOobject, a mesh or a registry.
    The legacy laws construct their own fields, which is exactly the coupling
    this framework removes.

    Two roles, which differ in direction:

      - persistent: history the law writes and reads back, such as the plastic
        strain. Owned by the law, and on restart it must come back.

      - prescribed: something the user supplies and the law only reads, such as
        an initial stress or a fibre direction. Never written back. Absent
        means the declared default is used, which is what makes adding one to
        an existing law a no-op for every case that does not supply it.

SourceFiles
    None: header only

\*---------------------------------------------------------------------------*/

#ifndef mechanicalConstitutiveLawStateSpec_H
#define mechanicalConstitutiveLawStateSpec_H

#include "word.H"
#include "DynamicList.H"
#include "scalar.H"
#include "vector.H"
#include "tensor.H"
#include "symmTensor.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
            Class mechanicalConstitutiveLawStateSpec Declaration
\*---------------------------------------------------------------------------*/

class mechanicalConstitutiveLawStateSpec
{
public:

    //- Role of a declared state field
    enum class stateRole
    {
        //- History the law writes and reads back
        persistent,

        //- Supplied by the user as a field of the same name, read only,
        //  never written back
        prescribed,

        //- Constant at its declared value. No field is consulted, which is
        //  how a law says that something else - typically an explicit entry
        //  in its own dictionary - has already settled the value
        fixed
    };

    //- Name of a role, for messages
    static word stateRoleName(const stateRole r)
    {
        switch (r)
        {
            case stateRole::persistent: return "persistent";
            case stateRole::prescribed: return "prescribed";
            case stateRole::fixed: return "fixed";
        }

        return "unknown";
    }

    //- One declared field. The default is held as each possible type rather
    //  than a variant, because a state field is one of exactly four types and
    //  four members cost less than the machinery to avoid them
    struct entry
    {
        word name;
        word typeName;
        stateRole role;

        scalar scalarDefault;
        vector vectorDefault;
        tensor tensorDefault;
        symmTensor symmTensorDefault;
    };


private:

        //- The declarations, in the order the law made them
        DynamicList<entry> entries_;

        //- Add a declaration, rejecting a repeated name
        inline void append(const entry& e)
        {
            forAll(entries_, i)
            {
                if (entries_[i].name == e.name)
                {
                    FatalErrorInFunction
                        << "State field '" << e.name << "' is declared twice "
                        << "by the same law." << nl
                        << "A repeated name would silently share storage "
                        << "between two things the law believes are separate."
                        << exit(FatalError);
                }
            }

            entries_.append(e);
        }


public:

    // Constructors

        mechanicalConstitutiveLawStateSpec()
        :
            entries_()
        {}


    // Member Functions

        //- Declare a scalar state field
        inline void addScalar
        (
            const word& name,
            const stateRole role,
            const scalar deflt = 0.0
        )
        {
            entry e;
            e.name = name;
            e.typeName = "scalar";
            e.role = role;
            e.scalarDefault = deflt;
            append(e);
        }

        //- Declare a vector state field, such as a fibre direction
        inline void addVector
        (
            const word& name,
            const stateRole role,
            const vector& deflt = vector::zero
        )
        {
            entry e;
            e.name = name;
            e.typeName = "vector";
            e.role = role;
            e.vectorDefault = deflt;
            append(e);
        }

        //- Declare a tensor state field
        inline void addTensor
        (
            const word& name,
            const stateRole role,
            const tensor& deflt = tensor::zero
        )
        {
            entry e;
            e.name = name;
            e.typeName = "tensor";
            e.role = role;
            e.tensorDefault = deflt;
            append(e);
        }

        //- Declare a symmetric tensor state field, such as an initial stress
        inline void addSymmTensor
        (
            const word& name,
            const stateRole role,
            const symmTensor& deflt = symmTensor::zero
        )
        {
            entry e;
            e.name = name;
            e.typeName = "symmTensor";
            e.role = role;
            e.symmTensorDefault = deflt;
            append(e);
        }

        //- The declarations
        inline const UList<entry>& entries() const
        {
            return entries_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
