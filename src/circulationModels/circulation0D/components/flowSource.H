/*---------------------------------------------------------------------------*\
  Author: David Kelly
  License: GNU General Public License v3 (as solids4foam)

Class
    Foam::flowSource

Description
    Prescribes the flow from "in" to "out" whatever the pressure; see
    circulationWaveform. The thesis DrivenFlow.

        type flowSource;
        waveform { type sine; mean 100; amplitude 100; period 0.8333; phase -1.5708; }

\*---------------------------------------------------------------------------*/

#ifndef flowSource_H
#define flowSource_H

#include "circulationComponent.H"
#include "circulationWaveform.H"

namespace Foam
{

class flowSource
:
    public circulationComponent
{
    circulationWaveform wave_;

public:

    TypeName("flowSource");

    flowSource(const word& name, const dictionary& dict, const circulationContext& ctx);

    virtual role kind() const { return role::fixedFlow; }
    virtual scalar fixedFlow(const scalar t) const { return wave_.value(t); }
};

} // End namespace Foam

#endif
