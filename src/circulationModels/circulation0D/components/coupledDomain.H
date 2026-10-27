/*---------------------------------------------------------------------------*\
  Author: David Kelly
  License: GNU General Public License v3 (as solids4foam)

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

    coupledDomain(const word& name, const dictionary& dict, const circulationContext& ctx);

    virtual role kind() const { return role::composite; }
    virtual void expand(circulationBuilder& builder) const;
};

} // End namespace Foam

#endif
