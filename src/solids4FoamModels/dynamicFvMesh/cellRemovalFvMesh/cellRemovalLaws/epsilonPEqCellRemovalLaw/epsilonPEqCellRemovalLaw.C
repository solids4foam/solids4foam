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

#include "epsilonPEqCellRemovalLaw.H"
#include "addToRunTimeSelectionTable.H"
#include "zeroGradientFvPatchFields.H"
#include "fvc.H"
#include "removeCells.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(epsilonPEqCellRemovalLaw, 0);
    addToRunTimeSelectionTable
    (
        cellRemovalLaw, epsilonPEqCellRemovalLaw, dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from dictionary
Foam::epsilonPEqCellRemovalLaw::epsilonPEqCellRemovalLaw
(
    const word& name,
    fvMesh& mesh,
    const dictionary& dict
)
:
    cellRemovalLaw(name, mesh, dict),
    epsilonPEqCrit_(readScalar(dict.lookup("epsilonPEqCritical"))),
    epsilonPEqName_(dict.lookupOrDefault<word>("epsilonPEqName", "epsilonPEq")),
    patchID_(readInt(dict.lookup("exposedFacesPatchID")))
{}


// * * * * * * * * * * * * * * * * Destructor * * * * * * * * * * * * * * * * //

Foam::epsilonPEqCellRemovalLaw::~epsilonPEqCellRemovalLaw()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::labelField> Foam::epsilonPEqCellRemovalLaw::cellsToRemove()
{
    // Lookup the plastic equivalent strain
    if (mesh().foundObject<volScalarField>(epsilonPEqName_))
    {
        const volScalarField& epsilonPEq =
            mesh().lookupObject<volScalarField>(epsilonPEqName_);

        // Find cells with epsilonPEq greater than the critical value
        const scalarField& epsilonPEqI = epsilonPEq.internalField();

        labelHashSet cellsToRemove;

        forAll(epsilonPEqI, cellI)
        {
            if (epsilonPEqI[cellI] > epsilonPEqCrit_)
            {
                cellsToRemove.insert(cellI);
            }
        }

        return tmp<labelField>(new labelField(cellsToRemove.toc()));
    }
    else
    {
        return tmp<labelField>(new labelField(0));
    }
}


Foam::label Foam::epsilonPEqCellRemovalLaw::exposedFacesPatchID()
{
    return patchID_;
}

// ************************************************************************* //
