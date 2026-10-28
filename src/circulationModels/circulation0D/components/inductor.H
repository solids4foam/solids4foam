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
    Foam::inductor

Description
    Blood inertia, L dq/dt = Pin - Pout. The flow q is a state.

        type inductor;
        L    0.0017;        // pressure*s^2/volume
        initialFlow 0;      // optional, volume/s

Author
    David Kelly, UCD.

\*---------------------------------------------------------------------------*/

#ifndef inductor_H
#define inductor_H

#include "circulationComponent.H"

namespace Foam
{

class inductor
:
    public circulationComponent
{
    scalar L_;
    scalar q0_;

public:

    TypeName("inductor");

    inductor
    (
        const word& name,
        const dictionary& dict,
        const circulationContext& ctx
    );

    virtual void parameters(DynamicList<parameter>&) const;

    virtual role kind() const { return role::inertia; }
    virtual scalar inertance() const { return L_; }
    virtual scalar initialFlow() const { return q0_; }
};

} // End namespace Foam

#endif

// ************************************************************************* //
