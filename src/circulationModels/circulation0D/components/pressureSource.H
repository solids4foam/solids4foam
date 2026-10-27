/*---------------------------------------------------------------------------*\
  Author: David Kelly
  License: GNU General Public License v3 (as solids4foam)

Class
    Foam::pressureSource

Description
    A node whose pressure is prescribed in time; see circulationWaveform.

        type pressureSource;
        waveform { type sine; mean 80; amplitude 40; period 0.8; }

\*---------------------------------------------------------------------------*/

#ifndef pressureSource_H
#define pressureSource_H

#include "circulationComponent.H"
#include "circulationWaveform.H"

namespace Foam
{

class pressureSource
:
    public circulationComponent
{
    circulationWaveform wave_;

public:

    TypeName("pressureSource");

    pressureSource(const word& name, const dictionary& dict, const circulationContext& ctx);

    virtual role kind() const { return role::fixedPressure; }
    virtual scalar fixedPressure(const scalar t) const { return wave_.value(t); }
};

} // End namespace Foam

#endif
