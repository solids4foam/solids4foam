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
    Foam::resistor

Description
    Ohm's law, Q = (Pin - Pout)/R.

        type resistor;
        R    0.5;           // pressure*s/volume

Author
    David Kelly, UCD.

\*---------------------------------------------------------------------------*/

#ifndef resistor_H
#define resistor_H

#include "circulationComponent.H"

namespace Foam
{

class resistor
:
    public circulationComponent
{
    scalar R_;      //!< [Pa s/m3]

public:

    TypeName("resistor");

    resistor
    (
        const word& name,
        const dictionary& dict,
        const circulationContext& ctx
    );

    virtual void parameters(DynamicList<parameter>&) const;

    virtual role kind() const { return role::branch; }
    virtual bool linear() const { return true; }

    virtual scalar flow(const scalar, const scalar dP) const
    {
        return dP/R_;
    }

    virtual scalar dFlowdDP(const scalar, const scalar) const
    {
        return 1/R_;
    }
};

} // End namespace Foam

#endif

// ************************************************************************* //
