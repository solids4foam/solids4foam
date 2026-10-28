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
    Foam::coupledDomain

Description
    The 3D (or any external) replacement for a valvedChamber. It has the
    same two terminals, "in" and "out", so swapping a lumped ventricle for a
    CFD ventricle is a change of "type" and nothing else:

        LV { type valvedChamber; chamber {...}; inletValve {...}; ... }
        LV { type coupledDomain; provider external;
             inlet  { patch mitral; valve { CQ 400; reverseFraction 0.005; } }
             outlet { patch aorta;  valve { CQ 350; reverseFraction 0.005; } } }

    The chamber itself is outside the network; what the network keeps is the
    two valves and everything beyond them. At each of the two patches the
    domain gives the flow and receives the pressure to impose, through an
    externalPort (<name>.inlet, <name>.outlet) with the valve law of the
    case between the lumped node and the patch.

    provider external   the coupling code supplies the port flows each step
    provider table      each port replays a waveform table (verification of
                        the seam without a CFD run)

Author
    David Kelly, UCD.

\*---------------------------------------------------------------------------*/

#ifndef coupledDomain_H
#define coupledDomain_H

#include "circulationComponent.H"

namespace Foam
{

class coupledDomain
:
    public circulationComponent
{
public:

    TypeName("coupledDomain");

    coupledDomain
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
