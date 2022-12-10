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

Description
    Dugdale simple cohesive law.

\*---------------------------------------------------------------------------*/

#include "DugdaleSimpleCohesiveZoneLaw.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(DugdaleSimpleCohesiveZoneLaw, 0);
    addToRunTimeSelectionTable
    (
        simpleCohesiveZoneLaw,
        DugdaleSimpleCohesiveZoneLaw,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::DugdaleSimpleCohesiveZoneLaw::DugdaleSimpleCohesiveZoneLaw
(
    const word& cohesiveLawName,
    const dictionary& dict
)
:
    simpleCohesiveZoneLaw(cohesiveLawName, dict),
    deltaC_(GIc()/sigmaMax())
{}


Foam::DugdaleSimpleCohesiveZoneLaw::DugdaleSimpleCohesiveZoneLaw
(
    const DugdaleSimpleCohesiveZoneLaw& dcl
)
:
    simpleCohesiveZoneLaw(dcl),
    deltaC_(dcl.deltaC_)
{}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::DugdaleSimpleCohesiveZoneLaw::~DugdaleSimpleCohesiveZoneLaw()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

//- Return current holding traction
Foam::scalar Foam::DugdaleSimpleCohesiveZoneLaw::traction(scalar delta) const
{
    if (delta > deltaC().value())
    {
        return 0.0;
    }
    else if (delta < 0)
    {
        return sigmaMax().value();
    }

    return sigmaMax().value();
}

// ************************************************************************* //
