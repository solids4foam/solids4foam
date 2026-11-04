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

#include "fluidMeshMotionCoupler.H"

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::fluidMeshMotionCoupler> Foam::fluidMeshMotionCoupler::New
(
    dynamicFvMesh& mesh,
    const labelList& patchIDs,
    const dictionary& dict
)
{
    word couplerType(dict.lookupOrDefault<word>("type", "auto"));

    if (couplerType == "auto")
    {
        couplerType = autoSelectType(mesh);
    }

    Info<< "Selecting fluid mesh motion coupler: " << couplerType << endl;

#ifdef OPENFOAM_COM
    auto* ctorPtr = dictionaryConstructorTable(couplerType);

    if (!ctorPtr)
    {
        FatalErrorInFunction
            << "Unknown fluidMeshMotion type "
            << couplerType << nl << nl
            << "Valid fluidMeshMotion types are:" << nl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }
#else
    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(couplerType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown fluidMeshMotion type "
            << couplerType << nl << nl
            << "Valid fluidMeshMotion types are:" << nl
            << dictionaryConstructorTablePtr_->toc()
            << exit(FatalError);
    }

    auto* ctorPtr = cstrIter();
#endif

    return autoPtr<fluidMeshMotionCoupler>(ctorPtr(mesh, patchIDs, dict));
}


// ************************************************************************* //
