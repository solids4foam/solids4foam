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
    Foam::pressureSource

Description
    A node whose pressure is prescribed in time; see circulationWaveform.

        type pressureSource;
        waveform { type sine; mean 80; amplitude 40; period 0.8; }

Author
    David Kelly, UCD.

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

    pressureSource
    (
        const word& name,
        const dictionary& dict,
        const circulationContext& ctx
    );

    virtual role kind() const { return role::fixedPressure; }
    virtual scalar fixedPressure(const scalar t) const
    {
        return wave_.value(t);
    }
};

} // End namespace Foam

#endif

// ************************************************************************* //
