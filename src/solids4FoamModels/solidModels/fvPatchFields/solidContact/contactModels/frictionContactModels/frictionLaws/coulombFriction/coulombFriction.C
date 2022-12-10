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

\*---------------------------------------------------------------------------*/

#include "coulombFriction.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(coulombFriction, 0);
    addToRunTimeSelectionTable(frictionLaw, coulombFriction, dictionary);



// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from dictionary
coulombFriction::coulombFriction
(
    const word& name,
    const frictionContactModel& fricModel,
    const dictionary& dict
)
:
  frictionLaw(name, fricModel, dict),
  frictionLawDict_(dict.subDict("frictionLawDict")),
  frictionCoeff_(readScalar(frictionLawDict_.lookup("frictionCoeff")))
{}


// Construct as a copy
coulombFriction::coulombFriction
(
    const coulombFriction& fricLaw
)
:
  frictionLaw(fricLaw),
  frictionLawDict_(fricLaw.frictionLawDict_),
  frictionCoeff_(fricLaw.frictionCoeff_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

coulombFriction::~coulombFriction()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

scalar coulombFriction::slipTraction(const scalar pressure)
{
    return frictionCoeff_*pressure;
}


scalar coulombFriction::slipTraction(const scalar pressure, const vector&)
{
    return frictionCoeff_*pressure;
}


void coulombFriction::writeDict(Ostream& os) const
{
    word keyword("frictionLawDict");

    os.writeKeyword(keyword)
        << frictionLawDict_;
}

// ************************************************************************* //

} // end of namespace foam
