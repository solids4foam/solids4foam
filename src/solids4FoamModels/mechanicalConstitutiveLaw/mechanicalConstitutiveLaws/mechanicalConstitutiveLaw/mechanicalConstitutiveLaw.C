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

InClass
    Foam::mechanicalConstitutiveLaw

\*---------------------------------------------------------------------------*/

#include "mechanicalConstitutiveLaw.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(mechanicalConstitutiveLaw, 0);
    defineRunTimeSelectionTable
    (
        mechanicalConstitutiveLaw, mechanicalConstitutiveLaw
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::mechanicalConstitutiveLaw::mechanicalConstitutiveLaw
(
    const dictionary& dict
)
{
    // Suppress unused-parameter warning
    (void)dict;
}


// ************************************************************************* //
