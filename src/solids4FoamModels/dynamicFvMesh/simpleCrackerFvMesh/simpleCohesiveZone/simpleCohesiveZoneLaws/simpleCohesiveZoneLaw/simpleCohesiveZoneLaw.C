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
    Virtual base class for cohesive law.

\*---------------------------------------------------------------------------*/

#include "simpleCohesiveZoneLaw.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(simpleCohesiveZoneLaw, 0);
    defineRunTimeSelectionTable(simpleCohesiveZoneLaw, dictionary);
}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::simpleCohesiveZoneLaw> Foam::simpleCohesiveZoneLaw::New
(
    const word& modelType,
    const dictionary& dict
)
{
    if (debug)
    {
        Info << "Selecting cohesive law: " << modelType << endl;
    }

#if (OPENFOAM >= 2112)
    auto* ctorPtr = dictionaryConstructorTable(modelType);

    if (!ctorPtr)
    {
        FatalIOErrorInLookup
        (
            dict,
            "simpleCohesiveZoneLaw",
            modelType,
            *dictionaryConstructorTablePtr_
        ) << exit(FatalIOError);
    }

#else
    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(modelType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "simpleCohesiveZoneLaw::New(const word& simpleCohesiveZoneLawName, "
            "const dictionary& dict)"
        )   << "Unknown cohesive law " << modelType
            << endl << endl
            << "Valid cohesive laws are :" << endl
            << dictionaryConstructorTablePtr_->toc()
            << exit(FatalError);
    }

    auto* ctorPtr = cstrIter();
#endif

    return autoPtr<simpleCohesiveZoneLaw>(ctorPtr(modelType, dict));
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::simpleCohesiveZoneLaw::simpleCohesiveZoneLaw
(
    const word& simpleCohesiveZoneLawName,
    const dictionary& dict
)
:
    simpleCohesiveZoneLawCoeffs_
    (
        dict.subDict(simpleCohesiveZoneLawName + "Coeffs")
    ),
    GIc_(simpleCohesiveZoneLawCoeffs_.lookup("GIc")),
    sigmaMax_(simpleCohesiveZoneLawCoeffs_.lookup("sigmaMax"))
{}


Foam::simpleCohesiveZoneLaw::simpleCohesiveZoneLaw
(
    const simpleCohesiveZoneLaw& cl
)
:
    refCount(),
    simpleCohesiveZoneLawCoeffs_(cl.simpleCohesiveZoneLawCoeffs_),
    GIc_(cl.GIc_),
    sigmaMax_(cl.sigmaMax_)
{}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::simpleCohesiveZoneLaw::~simpleCohesiveZoneLaw()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::simpleCohesiveZoneLaw::writeDict(Ostream& os) const
{
    os.writeKeyword(word(type() + "Coeffs"))
        << simpleCohesiveZoneLawCoeffs();
}


// ************************************************************************* //
