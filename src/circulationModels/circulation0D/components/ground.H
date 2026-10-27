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
    Foam::ground

Description
    A node held at a fixed pressure (default zero), the reference of a
    Windkessel.

        type ground;
        P    0;             // optional

Author
    David Kelly, UCD.

\*---------------------------------------------------------------------------*/

#ifndef ground_H
#define ground_H

#include "circulationComponent.H"

namespace Foam
{

class ground
:
    public circulationComponent
{
    scalar P_;

public:

    TypeName("ground");

    ground
    (
        const word& name,
        const dictionary& dict,
        const circulationContext& ctx
    );

    virtual role kind() const { return role::fixedPressure; }
    virtual scalar fixedPressure(const scalar) const { return P_; }
};

} // End namespace Foam

#endif

// ************************************************************************* //
