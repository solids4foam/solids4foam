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

#include "denseMatrixHelperFunctions.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

tmp<scalarField> magSqr(const Field<scalarRectangularMatrix>& mats)
{
    tmp<scalarField> tresult(new scalarField(mats.size(), 0.0));
    scalarField& result = tresult.ref();

    forAll(mats, i)
    {
        result[i] = mats[i].norm(true);
    }

    return tresult;
}


} // End namespace Foam


// ************************************************************************* //
