/*---------------------------------------------------------------------------*\
  Author: David Kelly
  License: GNU General Public License v3 (as solids4foam)

Class
    Foam::CRL

Description
    compliance, resistor, inductor: in -> C -> R -> L -> out (thesis CRL). A vessel segment with storage, loss and inertia.

        type CRL;
        C   0.08;  R   0.003;  L   6.2e-5;
        V0 0; P0 0;                 // optional, see compliance
        initialPressure 100;        // optional, see compliance
        initialFlow 0;              // optional, see inductor

    Its parts are named <name>.C, <name>.R (and <name>.L), so their results
    appear as <name>.C.P, <name>.C.V, <name>.R.Q and so on.

\*---------------------------------------------------------------------------*/

#ifndef CRL_H
#define CRL_H

#include "circulationComponent.H"

namespace Foam
{

class CRL
:
    public circulationComponent
{
public:

    TypeName("CRL");

    CRL(const word& name, const dictionary& dict, const circulationContext& ctx);

    virtual role kind() const { return role::composite; }
    virtual void expand(circulationBuilder& builder) const;
};

} // End namespace Foam

#endif
