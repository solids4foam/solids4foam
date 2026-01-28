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
    fieldWithHistory

Description
    Lightweight container for storing a Field with a configurable number of
    previous time times (current, old, old-old, etc.).

    Intended for non-geometric fields (e.g. vectorField) where oldTime()
    functionality is required but GeometricField is not appropriate.

    Time history is advanced explicitly by calling advance(), which should
    be done once per time step (not per nonlinear iteration).

    Currently only implemented for Field<vector> - this can easily be
    templated.

Author
    Philip Cardiff, UCD.

\*---------------------------------------------------------------------------*/

#ifndef fieldWithHistory_H
#define fieldWithHistory_H

#include "Field.H"
#include "List.H"
#include "vector.H"
#include "fvPatchFieldMapper.H"
 
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class fieldWithHistory Declaration
\*---------------------------------------------------------------------------*/

class fieldWithHistory
:
    public Field<vector>
{
    //- Private data

        //- List of fields
        List<vectorField> oldTimes_;  // [0]=old, [1]=oldOld, ...

public:

    // Constructors

        //- Empty constructor
        fieldWithHistory() = default;

        //- Construct given the field size and number of old times
        fieldWithHistory
        (
            const label size,
            const vector defaultVal = vector::zero,
            const label nOldTimes = 0   // 0, 1, 2
        )
        :
            Field<vector>(size, defaultVal),
            oldTimes_(nOldTimes)
        {
            forAll(oldTimes_, i)
            {
                oldTimes_[i].setSize(size);
                oldTimes_[i] = defaultVal;
            }
        }

        //- Construct using a mapper
        fieldWithHistory
        (
            const fieldWithHistory& fld,
            const fvPatchFieldMapper& mapper
        )
        :
            Field<vector>(fld, mapper),
            oldTimes_(fld.nOldTimes())
        {
            forAll(oldTimes_, i)
            {
                oldTimes_[i] = fld.old(i);
            }
        }

    // Member functions

        // Access

            //- Current-time field
            vectorField& current()
            {
                return *this;
            }
            const vectorField& current() const
            {
                return *this;
            }

            //- Old-time field
            vectorField& old(label i = 0)
            {
                if (i >= oldTimes_.size())
                {
                    FatalErrorInFunction
                        << "old(" << i << ") requested but only "
                        << oldTimes_.size() << " old times stored"
                        << abort(FatalError);
                }

                return oldTimes_[i];
            }
            const vectorField& old(label i = 0) const
            {
                if (i >= oldTimes_.size())
                {
                    FatalErrorInFunction
                        << "old(" << i << ") requested but only "
                        << oldTimes_.size() << " old times stored"
                        << abort(FatalError);
                }

                return oldTimes_[i];
            }

            //- Old-old-time field
            vectorField& oldOld()
            {
                return old(1);
            }
            const vectorField& oldOld() const
            {
                return old(1);
            }

            //- Number of stored old times
            label nOldTimes() const
            {
                return oldTimes_.size();
            }

        // Edit

            //- Set the number of old times
            void storeOldTimes(const label i = 1)
            {
                oldTimes_.resize(i);

                forAll(oldTimes_, i)
                {
                    // Initialise with the current field
                    oldTimes_[i] = *this;
                }
            }

            //- Update the old fields: call once per time step
            void advance()
            {
                for (label i = oldTimes_.size() - 1; i > 0; --i)
                {
                    oldTimes_[i] = oldTimes_[i - 1];
                }

                if (oldTimes_.size() > 0)
                {
                    // Update old-time from the current field
                    oldTimes_[0] = *this;
                }
            }

            //- Map fields
            void autoMap
            (
                const fvPatchFieldMapper& mapper,
                const bool applyFlip = true
            )
            {
                Field<vector>::autoMap(mapper, true);

                forAll(oldTimes_, i)
                {
                    oldTimes_[i].autoMap(mapper, true);
                }
            }

        // Operators

            //- Copy assignment
            void operator=(const Field<vector>& rhs)
            {
                Field<vector>::operator=(rhs);
            }
            void operator=(const tmp<Field<vector>>& trhs)
            {
                Field<vector>::operator=(trhs);
            }            
            void operator=(const vector& rhs)
            {
                Field<vector>::operator=(rhs);
            }            
};

    
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
