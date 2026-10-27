/*---------------------------------------------------------------------------*\
  Author: David Kelly
  License: GNU General Public License v3 (as solids4foam)

Class
    Foam::CR

Description
    compliance, then resistor: in -> C -> R -> out (thesis CR). Windkessel-type vessel.

        type CR;
        C   20.5;  R   0.075;
        V0 0; P0 0;                 // optional, see compliance
        initialPressure 100;        // optional, see compliance

    Its parts are named <name>.C, <name>.R (and <name>.L), so their results
    appear as <name>.C.P, <name>.C.V, <name>.R.Q and so on.

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

    CR(const word& name, const dictionary& dict, const circulationContext& ctx);

    virtual role kind() const { return role::composite; }
    virtual void expand(circulationBuilder& builder) const;
};

} // End namespace Foam

#endif
