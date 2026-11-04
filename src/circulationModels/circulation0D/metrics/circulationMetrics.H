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
    that starts or stops in mid-cycle does not report the partial one. A cycle
    boundary between two samples is given values interpolated linearly to it,
    so the output interval need not divide T.

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
#include "OFstream.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class circulationMetrics Declaration
\*---------------------------------------------------------------------------*/

class circulationMetrics
{
    // Private Data Types

        //- A metric from the dictionary: its name and type, and the output
        //  columns it monitors
        struct item
        {
            //- Name of the metric, the prefix of its output headings
            word name;

            //- Type of the metric: chamber, pressure or flow
            word type;

            //- Column of the volume (chamber), or -1
            label cVolume = -1;

            //- Column of the pressure (chamber, pressure), or -1
            label cPressure = -1;

            //- Column of the flow (flow), or -1
            label cFlow = -1;
        };


    // Private Data

        //- The metrics to evaluate
        List<item> items_;

        //- Cardiac period [s]
        scalar period_;

        //- Directory that metrics.dat and the pressure-volume loops are
        //  written to
        fileName dir_;

        //- Stream of metrics.dat
        autoPtr<OFstream> out_;

        //- Whether the header of metrics.dat has been written
        bool headerWritten_;

        //- Index k of the cycle in progress (-1 before the first sample)
        label cycle_;

        //- Whether the cycle in progress was sampled from its first instant
        //  (only such cycles are reported)
        bool cycleStartedOnBoundary_;

        //- Sample times of the cycle in progress
        DynamicList<scalar> t_;

        //- Samples of the cycle in progress, one list per monitored column
        List<DynamicList<scalar>> v_;

        //- Column indices of the samples kept
        labelList monitored_;

        //- Unit labels for the header
        string volumeLabel_, pressureLabel_, flowLabel_;


    // Private Member Functions

        //- Index in v_ of the samples of an output column (-1 if the
        //  column is not monitored)
        label slot(const label column) const;

        //- Evaluate the metrics of the cycle in progress, append its row to
        //  metrics.dat and write the pressure-volume loops of the chambers;
        //  a cycle not sampled from its first instant is skipped
        void finishCycle();

        //- Start cycle k, clearing the samples
        void startCycle(const label k);


public:

    // Constructors

        //- Construct from the metrics dictionary, the names of the output
        //  columns, the context and the output directory
        circulationMetrics
        (
            const dictionary& dict,
            const wordList& columns,
            const circulationContext& ctx,
            const fileName& outDir
        );


    // Member Functions

        //- Add a sample (t, the output values in the order of the columns)
        void add(const scalar t, const scalarField& values);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
