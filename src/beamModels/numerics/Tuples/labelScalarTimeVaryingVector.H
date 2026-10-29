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

Typedef
    Foam::labelScalarTimeVaryingVector

Description
    Tuple2 labelScalarTimeVaryingVector

SourceFiles

\*---------------------------------------------------------------------------*/

#ifndef labelScalarTimeVaryingVector_H
#define labelScalarTimeVaryingVector_H

#include "labelScalar.H"
#include "timeVaryingVector.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    typedef Tuple2<labelScalar, timeVaryingVector> labelScalarTimeVaryingVector;
    typedef List<labelScalarTimeVaryingVector> labelScalarTimeVaryingVectorList;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
