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
    Linear cohesive law.

\*---------------------------------------------------------------------------*/

#include "linearSimpleCohesiveZoneLaw.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(linearSimpleCohesiveZoneLaw, 0);
    addToRunTimeSelectionTable
    (
        simpleCohesiveZoneLaw,
        linearSimpleCohesiveZoneLaw,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::linearSimpleCohesiveZoneLaw::linearSimpleCohesiveZoneLaw
(
    const word& cohesiveLawName,
    const dictionary& dict
)
:
    simpleCohesiveZoneLaw(cohesiveLawName, dict),
    deltaC_(2.0*GIc()/sigmaMax())
{}


Foam::linearSimpleCohesiveZoneLaw::linearSimpleCohesiveZoneLaw
(
    const linearSimpleCohesiveZoneLaw& lcl
)
:
    simpleCohesiveZoneLaw(lcl),
    deltaC_(lcl.deltaC_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::linearSimpleCohesiveZoneLaw::~linearSimpleCohesiveZoneLaw()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

//- Return current holding traction
Foam::scalar Foam::linearSimpleCohesiveZoneLaw::traction(scalar delta) const
{
    if (delta > deltaC().value())
    {
        return 0.0;
    }
    else if (delta < 0)
    {
        return sigmaMax().value();
    }

    return sigmaMax().value()*(1.0 - delta/deltaC().value());
}

// ************************************************************************* //
