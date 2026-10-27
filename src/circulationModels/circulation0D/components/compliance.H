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
    Foam::compliance

Description
    A compliant vessel: a node that stores volume, P = (V - V0)/C + P0. The
    volume is a state, dV/dt = flow in - flow out.

        type compliance;
        C    1.6;           // volume/pressure
        V0   0;             // unstressed volume, optional (default 0)
        P0   0;             // pressure offset, optional (default 0)
        initialPressure 100;  // or initialVolume; default is V = V0

Author
    David Kelly, UCD.

\*---------------------------------------------------------------------------*/

#ifndef compliance_H
#define compliance_H

#include "circulationComponent.H"

namespace Foam
{

class compliance
:
    public circulationComponent
{
    scalar C_;
    scalar V0_;
    scalar P0_;
    scalar Vinit_;

public:

    TypeName("compliance");

    compliance
    (
        const word& name,
        const dictionary& dict,
        const circulationContext& ctx
    );

    virtual void parameters(DynamicList<parameter>&) const;

    virtual role kind() const { return role::storage; }

    virtual scalar pressure(const scalar, const scalar V) const
    {
        return (V - V0_)/C_ + P0_;
    }

    virtual scalar initialVolume() const { return Vinit_; }
};

} // End namespace Foam

#endif

// ************************************************************************* //
