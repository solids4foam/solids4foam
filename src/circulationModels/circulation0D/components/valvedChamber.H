/*---------------------------------------------------------------------------*\
  Author: David Kelly
  License: GNU General Public License v3 (as solids4foam)

Class
    Foam::valvedChamber

Description
    A heart chamber with its inlet and outlet valves, in -> inletValve ->
    chamber -> outletValve -> out: the unit that is swapped for a 3D domain.
    Its two terminals are the same as those of coupledDomain, so replacing
    one by the other changes nothing in "connections".

        type valvedChamber;
        chamber     { V0 5; P0 1; Emin 0.1; Emax 2.5; Tes 0.3; Ted 0.45;
                      initialVolume 500; }          // see elastanceChamber
        inletValve  { CQ 400; }                     // see orificeValve
        outletValve { CQ 350; }

    Its parts are named <name>.inletValve, <name>.chamber and
    <name>.outletValve.

\*---------------------------------------------------------------------------*/

#ifndef valvedChamber_H
#define valvedChamber_H

#include "circulationComponent.H"

namespace Foam
{

class valvedChamber
:
    public circulationComponent
{
public:

    TypeName("valvedChamber");

    valvedChamber(const word& name, const dictionary& dict, const circulationContext& ctx);

    virtual role kind() const { return role::composite; }
    virtual void expand(circulationBuilder& builder) const;
};

} // End namespace Foam

#endif
