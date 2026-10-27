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
    Foam::flowSource

Description
    Prescribes the flow from "in" to "out" whatever the pressure; see
    circulationWaveform. The thesis DrivenFlow.

        type flowSource;
        waveform { type sine; mean 100; amplitude 100; period 0.8333;
                   phase -1.5708; }

Author
    David Kelly, UCD.

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

    flowSource
    (
        const word& name,
        const dictionary& dict,
        const circulationContext& ctx
    );

    virtual role kind() const { return role::fixedFlow; }
    virtual scalar fixedFlow(const scalar t) const { return wave_.value(t); }
};

} // End namespace Foam

#endif

// ************************************************************************* //
