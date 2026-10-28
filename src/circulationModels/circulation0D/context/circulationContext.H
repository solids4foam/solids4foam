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
    Foam::circulationContext

Description
    Global settings shared by every component of a circulation network: the
    unit system that case dictionaries are written in, and the cardiac period.

    Everything inside the library is SI (Pa, m3, s). Dictionaries may be
    written in the units of the literature instead, selected by

        units   SI;           // Pa, m3 (default)
        units   mmHg_mL;      // mmHg, mL; time is always in seconds

    or a custom system

        units   { pressure 133.322; volume 1e-6; }   // SI per unit

    Values are converted on reading, with the factor of their physical
    quantity (a resistance is pressure*time/volume, and so on), and results
    are written back out in the same units the case was written in.

Author
    David Kelly, UCD.

\*---------------------------------------------------------------------------*/

#ifndef circulationContext_H
#define circulationContext_H

#include "dictionary.H"
#include "scalar.H"
#include "word.H"
#include "string.H"
#include "dimensionSet.H"

namespace Foam
{

//- Built-in sub-circuit modules (windkessel2/3/4, healthyHuman, pipeLoopHeart)
dictionary circulationPresets();

class circulationContext
{
public:

    //- Physical quantities that carry units
    enum class quantity
    {
        pressure,           //!< P
        volume,             //!< V
        flow,               //!< Q = V/t
        resistance,         //!< R = P t / V
        compliance,         //!< C = V / P
        inertance,          //!< L = P t^2 / V
        elastance,          //!< E = P / V
        valveCoefficient,   //!< CQ = V / (t sqrt(P))
        quadraticResistance,//!< Rhat = P t^2 / V^2   (dP = Rhat Q^2)
        elastanceQuadratic, //!< E2 = P / V^2
        inverseVolume,      //!< K = 1 / V
        time,               //!< seconds, never converted
        none
    };


private:

    //- Sub-circuit templates: the built-in presets, overridden by the case's
    //  "modules" dictionary
    dictionary modules_;

    word name_;
    scalar pressureFactor_;     //!< Pa per pressure unit
    scalar volumeFactor_;       //!< m3 per volume unit
    word pressureLabel_;
    word volumeLabel_;
    scalar period_;             //!< cardiac period [s]


public:

    //- Read "units" and "Tcycle" from the top-level circulationProperties
    explicit circulationContext(const dictionary& dict);

    //- SI per user unit
    scalar factor(const quantity q) const;

    scalar toSI(const quantity q, const scalar v) const
    {
        return v*factor(q);
    }

    scalar fromSI(const quantity q, const scalar v) const
    {
        return v/factor(q);
    }

    //- The dimensions of a quantity in SI base units
    static dimensionSet dimensions(const quantity q);

    //- Read a required entry as SI. A plain number is in the case's units
    //  (converted); an entry with dimensions, "R [1 -4 -1 0 0 0 0] 10", is
    //  checked against the dimensions the parameter must have, and is SI.
    scalar read(const dictionary& d, const word& key, const quantity q) const;

    //- Read an optional entry (default given in user units) as SI
    scalar read
    (
        const dictionary& d,
        const word& key,
        const quantity q,
        const scalar defaultUserValue
    ) const
    {
        return
            d.found(key)
          ? read(d, key, q)
          : defaultUserValue*factor(q);
    }

    //- Unit label used in output headings, e.g. "mmHg", "mL/s"
    string label(const quantity q) const;

    const word& name() const { return name_; }

    //- Cardiac period [s]
    scalar period() const { return period_; }

    //- The module definition of that name, or null
    const dictionary* findModule(const word& type) const
    {
        return modules_.isDict(type) ? &modules_.subDict(type) : nullptr;
    }

    //- Names of all modules available
    wordList moduleNames() const { return modules_.sortedToc(); }
};

} // End namespace Foam

#endif

// ************************************************************************* //
