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
    mechanicalConstitutiveLaw

\*---------------------------------------------------------------------------*/

#include "mechanicalConstitutiveLaw.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

autoPtr<mechanicalConstitutiveLaw>
mechanicalConstitutiveLaw::New(const dictionary& dict)
{
    const word lawType(dict.lookup("type"));

    Info<< "Selecting mechanical constitutive law: "
        << lawType << endl;

    auto* ctorPtr =
        mechanicalConstitutiveLawConstructorTable(lawType);

    if (!ctorPtr)
    {
        FatalErrorInFunction
            << "Unknown mechanicalConstitutiveLaw type "
            << lawType << nl << nl
            << "Valid mechanicalConstitutiveLaw types are:" << nl
            << mechanicalConstitutiveLawConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<mechanicalConstitutiveLaw>(ctorPtr(dict));
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
