/*---------------------------------------------------------------------------*\
License
    This file is part of beamFoam.

    beamFoam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    beamFoam is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with beamFoam.  If not, see <http://www.gnu.org/licenses/>.

Description
    Run-time checks on the parallel decomposition of a beam region.

Author
    Amirhossein Taran, UCD

SourceFiles
    beamParallelChecks.C

\*---------------------------------------------------------------------------*/

#ifndef beamParallelChecks_H
#define beamParallelChecks_H

#include "fvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    //- Check that every cell of the given beam region lives on processor 0,
    //  and exit with a FatalError if it does not
    //
    //  Couplings between a beam region and a fluid region (the actuator-line
    //  drag contribution, the actuator-line fvOption and the mooring-line
    //  restraint) gather the beam geometry, tangents and forces on the master
    //  processor and broadcast them to all ranks, i.e. they treat the master's
    //  local slice of the beam as if it were the whole beam.  A beam region
    //  that is decomposed across processors is therefore sampled and forced
    //  using the master's slice alone.  That failure is silent rather than
    //  loud - Field::operator= resizes on assignment and the field size checks
    //  are only compiled in FULLDEBUG builds - hence this explicit check.
    //
    //  Note that the requirement is processor 0 specifically, not merely a
    //  single processor: if the beam were placed entirely on a non-master
    //  rank, the master would broadcast an empty beam to everyone.
    //
    //  This is a collective operation and must be called on all ranks, i.e.
    //  never from inside a Pstream::master() guard.
    //
    //  \param beamMesh the beam region mesh
    //  \param context  what requires the constraint, used in the error message
    void checkBeamOnMasterProcessor
    (
        const fvMesh& beamMesh,
        const string& context
    );

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
