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

\*---------------------------------------------------------------------------*/

#include "beamParallelChecks.H"
#include "Pstream.H"
#include "PstreamReduceOps.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::checkBeamOnMasterProcessor
(
    const fvMesh& beamMesh,
    const string& context
)
{
    if (!Pstream::parRun())
    {
        return;
    }

    const label nLocal = beamMesh.nCells();

    // Collective reductions: all ranks must reach both of these
    const label nGlobal = returnReduce(nLocal, sumOp<label>());
    const label nOnMaster =
        returnReduce(Pstream::master() ? nLocal : label(0), sumOp<label>());

    if (nOnMaster != nGlobal)
    {
        FatalErrorInFunction
            // c_str() so the context reads as prose rather than a quoted
            // Foam::string
            << context.c_str() << " requires every cell of beam region '"
            << beamMesh.name() << "' to be on processor 0, but only "
            << nOnMaster << " of " << nGlobal << " beam cells are." << nl
            << nl
            << "The beam geometry, tangents and forces are gathered on the"
            << " master processor and broadcast to all ranks, so a beam"
            << " region that is split across processors would be sampled and"
            << " forced using the master's slice alone." << nl
            << nl
            << "Decompose the beam region with:" << nl
            << "    method       manual;" << nl
            << "    manualCoeffs { dataFile \"cellDecomposition\"; }" << nl
            << "and a constant/cellDecomposition labelList holding one entry"
            << " per beam control volume, all of them 0. Keep"
            << " numberOfSubdomains equal to that of the fluid region. See"
            << " tutorials/beamTunnel for a worked example." << nl
            << abort(FatalError);
    }
}


// ************************************************************************* //
