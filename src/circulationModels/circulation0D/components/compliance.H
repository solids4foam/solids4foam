/*---------------------------------------------------------------------------*\
  Author: David Kelly
  License: GNU General Public License v3 (as solids4foam)

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

    compliance(const word& name, const dictionary& dict, const circulationContext& ctx);

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
