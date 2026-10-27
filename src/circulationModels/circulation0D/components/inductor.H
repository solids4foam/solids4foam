/*---------------------------------------------------------------------------*\
  Author: David Kelly
  License: GNU General Public License v3 (as solids4foam)

Class
    Foam::inductor

Description
    Blood inertia, L dq/dt = Pin - Pout. The flow q is a state.

        type inductor;
        L    0.0017;        // pressure*s^2/volume
        initialFlow 0;      // optional, volume/s

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

    inductor(const word& name, const dictionary& dict, const circulationContext& ctx);

    virtual void parameters(DynamicList<parameter>&) const;

    virtual role kind() const { return role::inertia; }
    virtual scalar inertance() const { return L_; }
    virtual scalar initialFlow() const { return q0_; }
};

} // End namespace Foam

#endif
