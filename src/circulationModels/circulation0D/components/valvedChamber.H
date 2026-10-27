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

Author
    David Kelly, UCD.

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

    valvedChamber
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
