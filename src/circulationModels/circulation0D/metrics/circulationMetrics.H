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
    Foam::circulationMetrics

Description
    Per-cardiac-cycle summary quantities of a circulation run, so nobody has to
    re-write the same post-processing. Configured in circulationProperties by
    naming output columns (the "<element>.<P|V|Q>" headings of the results):

        metrics
        {
            ventricle { type chamber;  volume LV.V; pressure LV.P; }
            aorta     { type pressure; column SAT.C.P; }
            aortic    { type flow;     column AV.Q; }
        }

    Cycles are [k T, (k+1) T] for the cardiac period T = Tcycle. A cycle is
    reported when it has been run from its first instant to its last; a run
    that starts or stops in mid-cycle does not report the partial one.

    \table
        type      | per-cycle quantities
        chamber   | EDV (max V), ESV (min V), SV, EF = SV/EDV, CO = SV 60/T
                  | [volume/min], Pmax, Pmin
        pressure  | mean (time average), systolic (max), diastolic (min), pulse
        flow      | mean, max, min, net volume, forward volume, backward
                  | volume, regurgitant fraction
    \endtable

    One row per cycle is appended to metrics.dat; for each chamber the
    pressure-volume loop of the latest cycle is written to
    <name>.pvLoop.dat. Quantities are in the units of the case.

Author
    David Kelly, UCD.

SourceFiles
    circulationMetrics.C

\*---------------------------------------------------------------------------*/

#ifndef circulationMetrics_H
#define circulationMetrics_H

#include "circulationContext.H"
#include "dictionary.H"
#include "DynamicList.H"
#include "fileName.H"

#include <fstream>

namespace Foam
{

class circulationMetrics
{
    struct item
    {
        word name;
        word type;
        label cVolume = -1;
        label cPressure = -1;
        label cFlow = -1;
    };

    List<item> items_;
    scalar period_;
    fileName dir_;
    std::ofstream out_;
    bool headerWritten_;

    //- Samples of the cycle in progress: time, then the monitored columns
    label cycle_;                    //!< index k of the cycle in progress
    bool cycleStartedOnBoundary_;
    DynamicList<scalar> t_;
    List<DynamicList<scalar>> v_;    //!< one list per monitored column
    labelList monitored_;            //!< column indices of the samples kept

    //- Unit labels for the header
    string volumeLabel_, pressureLabel_, flowLabel_;

    label slot(const label column) const;

    void finishCycle();

    void startCycle(const label k);

public:

    circulationMetrics
    (
        const dictionary& dict,
        const wordList& columns,
        const circulationContext& ctx,
        const fileName& outDir
    );

    //- Add a sample (t, the output values in the order of the columns)
    void add(const scalar t, const scalarField& values);
};

} // End namespace Foam

#endif

// ************************************************************************* //
