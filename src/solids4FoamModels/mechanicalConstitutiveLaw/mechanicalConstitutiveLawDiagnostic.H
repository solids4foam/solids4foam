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
    Foam::mechanicalConstitutiveLawDiagnostic

Description
    One quantity a law reports at the end of a time step, together with the
    way it must be combined.

    A law reports what it did locally and reduces nothing. The manager gathers
    the reports from every state the law owns - the internal points and each
    boundary patch - and performs the collectives itself, once per quantity.

    The reason is that a law cannot safely reduce. The number of boundary
    patches is not the same on every rank, because processor patches exist
    only where the mesh was cut, so a law that called reduce() once per state
    would call a collective a different number of times on different ranks and
    the run would hang rather than fail. That is why boundary states went
    unvisited for as long as the laws did their own reducing.

    The combining operation travels with the value because it is not the
    manager's to know: a count of yielding points is a sum, a largest plastic
    increment is a maximum, and nothing about the number distinguishes them.

SourceFiles
    mechanicalConstitutiveLawDiagnostic.H

\*---------------------------------------------------------------------------*/

#ifndef mechanicalConstitutiveLawDiagnostic_H
#define mechanicalConstitutiveLawDiagnostic_H

#include "scalar.H"
// word.H rather than the string header: foam-extend calls that one
// foamString.H, and word.H pulls string in on every fork
#include "word.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class mechanicalConstitutiveLawDiagnostic Declaration
\*---------------------------------------------------------------------------*/

class mechanicalConstitutiveLawDiagnostic
{
public:

    //- How this quantity combines across states and across ranks
    enum class combineOperation
    {
        sum,
        maximum
    };


private:

    // Private data

        //- What the quantity is called when it is reported
        string name_;

        //- The value on this state, on this rank
        scalar value_;

        //- How to combine it with the same quantity from elsewhere
        combineOperation op_;


public:

    // Constructors

        //- Construct null. Only so that a DynamicList of these can be
        //  resized; a default-constructed one carries nothing
        mechanicalConstitutiveLawDiagnostic()
        :
            name_("undefined"),
            value_(0.0),
            op_(combineOperation::sum)
        {}

        //- Construct from name, value and operation
        mechanicalConstitutiveLawDiagnostic
        (
            const string& name,
            const scalar value,
            const combineOperation op
        )
        :
            name_(name),
            value_(value),
            op_(op)
        {}


    // Member Functions

        //- The quantity's name
        const string& name() const
        {
            return name_;
        }

        //- The value
        scalar value() const
        {
            return value_;
        }

        //- Non-const access, for the manager to combine into
        scalar& value()
        {
            return value_;
        }

        //- How it combines
        combineOperation op() const
        {
            return op_;
        }

        //- The value that leaves a combination unchanged
        scalar identity() const
        {
            // Both quantities reported so far are non-negative, so zero is
            // the identity for each. A signed quantity that had to be
            // maximised would need -VGREAT here instead, and would have to
            // say so rather than inherit this
            return 0.0;
        }
};


} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
