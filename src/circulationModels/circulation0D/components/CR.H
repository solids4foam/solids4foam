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
    Foam::CR

Description
    compliance, then resistor: in -> C -> R -> out. Windkessel-type vessel.

        type CR;
        C   20.5;  R   0.075;
        V0 0; P0 0;                 // optional, see compliance
        initialPressure 100;        // optional, see compliance

    Its parts are named <name>.C, <name>.R (and <name>.L), so their results
    appear as <name>.C.P, <name>.C.V, <name>.R.Q and so on.

Author
    David Kelly, UCD.

\*---------------------------------------------------------------------------*/

#ifndef CR_H
#define CR_H

#include "circulationComponent.H"

namespace Foam
{

class CR
:
    public circulationComponent
{
public:

    TypeName("CR");

    CR
    (
        const word& name,
        const dictionary& dict,
        const circulationContext& ctx
    );

    virtual role kind() const { return role::composite; }
    virtual void expand(circulationBuilder& builder) const;
};

} // End namespace Foam

#endif

// ************************************************************************* //
