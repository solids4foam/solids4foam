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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

//- Built-in sub-circuit modules (windkessel2/3/4, healthyHuman, pipeLoopHeart)
//  as a dictionary of module definitions; defined in circulationPresets.C
dictionary circulationPresets();


/*---------------------------------------------------------------------------*\
                      Class circulationContext Declaration
\*---------------------------------------------------------------------------*/

class circulationContext
{
public:

    // Public Data Types

        //- Physical quantities that carry units
        enum class quantity
        {
            //- P
            pressure,
            //- V
            volume,
            //- Q = V/t
            flow,
            //- R = P t / V
            resistance,
            //- C = V / P
            compliance,
            //- L = P t^2 / V
            inertance,
            //- E = P / V
            elastance,
            //- CQ = V / (t sqrt(P))
            valveCoefficient,
            //- Rhat = P t^2 / V^2   (dP = Rhat Q^2)
            quadraticResistance,
            //- E2 = P / V^2
            elastanceQuadratic,
            //- K = 1 / V
            inverseVolume,
            //- seconds, never converted
            time,
            //- Dimensionless, never converted
            none
        };


private:

    // Private Data

        //- Sub-circuit templates: the built-in presets, overridden by the
        //  case's "modules" dictionary
        dictionary modules_;

        //- Name of the unit system: SI, mmHg_mL or custom
        word name_;

        //- Pa per pressure unit
        scalar pressureFactor_;

        //- m3 per volume unit
        scalar volumeFactor_;

        //- Pressure unit label used in the output, e.g. "mmHg"
        word pressureLabel_;

        //- Volume unit label used in the output, e.g. "mL"
        word volumeLabel_;

        //- Cardiac period [s]
        scalar period_;


public:

    // Constructors

        //- Read "units", "Tcycle" and "modules" from the top-level
        //  circulationProperties
        explicit circulationContext(const dictionary& dict);


    // Member Functions

        //- SI per user unit of quantity q
        scalar factor(const quantity q) const;

        //- Convert a value of quantity q from the case's units to SI
        scalar toSI(const quantity q, const scalar v) const
        {
            return v*factor(q);
        }

        //- Convert an SI value of quantity q to the case's units
        scalar fromSI(const quantity q, const scalar v) const
        {
            return v/factor(q);
        }

        //- The dimensions of a quantity in SI base units
        static dimensionSet dimensions(const quantity q);

        //- Read a required entry as SI. A plain number is in the case's
        //  units (converted); an entry with dimensions,
        //  "R [1 -4 -1 0 0 0 0] 10", is checked against the dimensions the
        //  parameter must have, and is SI.
        scalar read
        (
            const dictionary& d,
            const word& key,
            const quantity q
        ) const;

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

        //- Name of the unit system: SI, mmHg_mL or custom
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


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
