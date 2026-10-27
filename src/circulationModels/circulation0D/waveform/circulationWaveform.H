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
    Foam::circulationWaveform

Description
    A prescribed function of time for sources: flowSource, pressureSource and
    the table provider of externalPort. Read from a "waveform" sub-dictionary
    of the component, in the units of the case:

        waveform { type constant; value 5; }
        waveform { type sine;     mean 100; amplitude 100;
                   period 0.8333; phase -1.5708; }
                   // mean + A sin(2 pi t/T + phase)
        waveform { type halfSine; amplitude 400; period 1;
                   systolicDuration 0.3; }
        waveform { type table;    file "constant/Qin.dat"; column 1;
                   periodic yes; }

    A table file is whitespace separated with the time in column 0 and "#"
    comment lines; "column" counts from the time column. Linear
    interpolation; constant extrapolation, or wrapped in time if periodic.

    If the component has no "waveform" entry but has "value", it is a
    constant.

Author
    David Kelly, UCD.

\*---------------------------------------------------------------------------*/

#ifndef circulationWaveform_H
#define circulationWaveform_H

#include "circulationContext.H"
#include "dictionary.H"
#include "DynamicList.H"
#include "scalarList.H"

namespace Foam
{

class circulationWaveform
{
public:

    enum class kind { constant, sine, halfSine, table };

private:

    kind kind_;
    scalar value_;          //!< constant value / mean   [SI]
    scalar amplitude_;      //!< [SI]
    scalar period_;         //!< [s]
    scalar phase_;          //!< [rad]
    scalar systolicDuration_;
    bool periodic_;
    scalarList tTable_;
    scalarList vTable_;

public:

    circulationWaveform
    (
        const dictionary& componentDict,
        const circulationContext& ctx,
        const circulationContext::quantity q
    );

    //- Value at time t [SI]
    scalar value(const scalar t) const;
};

} // End namespace Foam

#endif

// ************************************************************************* //
