/*---------------------------------------------------------------------------*\
  Author: David Kelly
  License: GNU General Public License v3 (as solids4foam)

Class
    Foam::resistor

Description
    Ohm's law, Q = (Pin - Pout)/R.

        type resistor;
        R    0.5;           // pressure*s/volume

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

    resistor(const word& name, const dictionary& dict, const circulationContext& ctx);

    virtual void parameters(DynamicList<parameter>&) const;

    virtual role kind() const { return role::branch; }
    virtual bool linear() const { return true; }
    virtual scalar flow(const scalar, const scalar dP) const { return dP/R_; }
    virtual scalar dFlowdDP(const scalar, const scalar) const { return 1/R_; }
};

} // End namespace Foam

#endif
