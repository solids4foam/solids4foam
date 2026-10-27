/*---------------------------------------------------------------------------*\
  Author: David Kelly
  License: GNU General Public License v3 (as solids4foam)

Class
    Foam::ground

Description
    A node held at a fixed pressure (default zero), the reference of a
    Windkessel.

        type ground;
        P    0;             // optional

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

    ground(const word& name, const dictionary& dict, const circulationContext& ctx);

    virtual role kind() const { return role::fixedPressure; }
    virtual scalar fixedPressure(const scalar) const { return P_; }
};

} // End namespace Foam

#endif
