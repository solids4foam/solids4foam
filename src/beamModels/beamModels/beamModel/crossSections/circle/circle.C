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

InClass
    circle

\*---------------------------------------------------------------------------*/

#include "circle.H"

#include "volFields.H"
#include "fvc.H"
#include "addToRunTimeSelectionTable.H"
#include "primitivePatchInterpolation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

  defineTypeNameAndDebug(circle, 0);
  addToRunTimeSelectionTable(crossSectionModel, circle, dictionary);

// * * * * * * * * * * * * Private Member Functions * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

circle::circle
(
    const word& name,
    const dictionary& dict
)
:
    crossSectionModel
    (
        name,
        dict
    ),
    crossSectionModelDict_(dict.subDict(name + "CrossSectionModelDict")),
    radius_
    (
        readScalar(crossSectionModelDict_.lookup("radius"))
    ),
    scalingArea_
    (
        crossSectionModelDict_.lookupOrDefault<scalar>("scalingArea", 1.0)
    ),
    scalingSecondMomentArea_
    (
        crossSectionModelDict_
        .lookupOrDefault<scalar>("scalingSecondMomentArea", 1.0)
    ),
    points_(0, vector::zero)
{}

// * * * * * * * * * * * * * * Member Functions * * * * * * * * * * * * * * * //

// tmp<vectorField> circle::greenLagrangianStrain
// (
//     const vector& Gamma,
//     const vector& K
// ) const
// {
//     tmp<vectorField> tE
//     (
//         new vectorField(points_.size(), vector::zero)
//     );

//     return tE;
// }


} // End namespace Foam

// ************************************************************************* //
