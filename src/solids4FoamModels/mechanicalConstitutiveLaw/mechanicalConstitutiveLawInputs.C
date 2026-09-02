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

Class

\*---------------------------------------------------------------------------*/

#include "mechanicalConstitutiveLawInputs.H"
#include "error.H"

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

const Foam::scalarField* Foam::mechanicalConstitutiveLawInputs::findScalar
(
    const word& name
) const
{
    HashTable<const scalarField*>::const_iterator iter =
        scalarFields_.find(name);

    if (iter == scalarFields_.end())
    {
        return nullptr;
    }

    return *iter;
}


const Foam::vectorField* Foam::mechanicalConstitutiveLawInputs::findVector
(
    const word& name
) const
{
    HashTable<const vectorField*>::const_iterator iter =
        vectorFields_.find(name);

    if (iter == vectorFields_.end())
    {
        return nullptr;
    }

    return *iter;
}


const Foam::scalarField& Foam::mechanicalConstitutiveLawInputs::getScalar
(
    const word& name
) const
{
    const scalarField* ptr = findScalar(name);

    if (!ptr)
    {
        FatalErrorInFunction
            << "The scalar input '" << name << "' was not supplied." << nl
            << "A law that requires an input must not read zero when it is "
            << "missing, so this is an error rather than a default. The solid "
            << "model must supply it before the evaluation." << nl
            << "Available scalar inputs: " << scalarFields_.toc()
            << exit(FatalError);
    }

    return *ptr;
}


const Foam::vectorField& Foam::mechanicalConstitutiveLawInputs::getVector
(
    const word& name
) const
{
    const vectorField* ptr = findVector(name);

    if (!ptr)
    {
        FatalErrorInFunction
            << "The vector input '" << name << "' was not supplied." << nl
            << "Available vector inputs: " << vectorFields_.toc()
            << exit(FatalError);
    }

    return *ptr;
}


void Foam::mechanicalConstitutiveLawInputs::setScalar
(
    const word& name,
    const scalarField& fld
)
{
    scalarFields_.set(name, &fld);
}


void Foam::mechanicalConstitutiveLawInputs::setVector
(
    const word& name,
    const vectorField& fld
)
{
    vectorFields_.set(name, &fld);
}


// ************************************************************************* //
