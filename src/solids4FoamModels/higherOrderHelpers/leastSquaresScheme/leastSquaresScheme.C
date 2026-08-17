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

#include "leastSquaresScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(leastSquaresScheme, 0);
defineRunTimeSelectionTable(leastSquaresScheme, dictionary);


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

autoPtr<leastSquaresScheme> leastSquaresScheme::New
(
    const fvMesh& mesh,
    const boolList& includePatchInStencils,
    const dictionary& dict
)
{
    const word schemeType
    (
        dict.lookupOrDefault<word>("type", "movingLeastSquares")
    );

    Info<< "Selecting least-squares reconstruction scheme "
        << schemeType << endl;

#if (OPENFOAM >= 2112)
    auto* ctorPtr = dictionaryConstructorTable(schemeType);

    if (!ctorPtr)
    {
        FatalIOErrorInLookup
        (
            dict,
            "leastSquaresScheme",
            schemeType,
            *dictionaryConstructorTablePtr_
        )   << exit(FatalIOError);
    }
#else
    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(schemeType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn("leastSquaresScheme::New(...)")
            << "Unknown least-squares reconstruction scheme type "
            << schemeType << nl << nl
            << "Valid least-squares reconstruction scheme types are:"
            << nl << dictionaryConstructorTablePtr_->toc()
            << exit(FatalError);
    }

    auto* ctorPtr = cstrIter();
#endif

    return autoPtr<leastSquaresScheme>
    (
        ctorPtr(mesh, includePatchInStencils, dict)
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * //

leastSquaresScheme::leastSquaresScheme(const fvMesh& mesh)
:
    mesh_(mesh)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * //

leastSquaresScheme::~leastSquaresScheme()
{}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
