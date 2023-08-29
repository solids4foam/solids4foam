/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     |
    \\  /    A nd           | For copyright notice see file Copyright
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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

#include "fv.H"
#include "HashTable.H"
#include "linear.H"
#include "fvMatrices.H"

#include "blockLaplacianTraceScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

// defineTypeNameAndDebug(blockLaplacianTrace, 0);
defineRunTimeSelectionTable(blockLaplacianTrace, Istream);

// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

tmp<blockLaplacianTrace> blockLaplacianTrace::New
(
    const fvMesh& mesh,
    Istream& schemeData
)
{
    if (fv::debug)
    {
        Info<< "blockLaplacianTrace::New(const fvMesh&, Istream&) : "
               "constructing blockLaplacianTrace"
            << endl;
    }

    if (schemeData.eof())
    {
        FatalIOErrorIn
        (
            "blockLaplacianTrace::New(const fvMesh&, Istream&)",
            schemeData
        )   << "fvmBlockLaplacianTrace scheme not specified" << endl << endl
            << "Valid fvmBlockLaplacianTrace schemes are :" << endl
            << IstreamConstructorTablePtr_->sortedToc()
            << exit(FatalIOError);
    }

    word schemeName(schemeData);

//     typename IstreamConstructorTable::iterator cstrIter =
    IstreamConstructorTable::iterator cstrIter =
        IstreamConstructorTablePtr_->find(schemeName);

    if (cstrIter == IstreamConstructorTablePtr_->end())
    {
        FatalIOErrorIn
        (
            "blockLaplacianTrace::New(const fvMesh&, Istream&)",
            schemeData
        )   << "unknown fvmBlockLaplacianTrace scheme " << schemeName
            << endl << endl
            << "Valid fvmBlockLaplacianTrace schemes are :" << endl
            << IstreamConstructorTablePtr_->sortedToc()
            << exit(FatalIOError);
    }

    return cstrIter()(mesh, schemeData);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

blockLaplacianTrace::~blockLaplacianTrace()
{}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
