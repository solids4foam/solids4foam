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

#include "mechanicalConstitutiveLawStateIO.H"

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class Type>
bool Foam::stateIOFieldProxy<Type>::writeData(Ostream& os) const
{
    // Gather the parts into one list in write order. The state is the
    // authoritative copy and is read here rather than mirrored, so there is no
    // second copy that could be stale by the time the registry writes
    Field<Type> all(mechanicalConstitutiveLawStateIO::totalSize(parts_));

    label n = 0;
    forAll(parts_, partI)
    {
        const Field<Type>& f =
            stateFieldAccess<Type>::get(*parts_[partI], variableName_);

        forAll(f, i)
        {
            all[n++] = f[i];
        }
    }

    os << all;

    return os.good();
}


// ************************************************************************* //
