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

#ifdef OPENFOAM_ORG

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class SourcePatch, class TargetPatch>
Foam::autoPtr<Foam::newAMIMethod<SourcePatch, TargetPatch>>
Foam::newAMIMethod<SourcePatch, TargetPatch>::New
(
    const word& methodName,
    const SourcePatch& srcPatch,
    const TargetPatch& tgtPatch,
    const scalarField& srcMagSf,
    const scalarField& tgtMagSf,
    const faceAreaIntersect::triangulationMode& triMode,
    const bool reverseTarget,
    const bool requireMatch
)
{
    if (debug)
    {
        Info<< "Selecting newAMIMethod " << methodName << endl;
    }

#if (OPENFOAM >= 2112)
    auto* ctorPtr = componentsConstructorTable(methodName);

    if (!ctorPtr)
    {
        FatalErrorInLookup
        (
            "newAMIMethod",
             methodName,
            *componentsConstructorTablePtr_
        ) << exit(FatalError);
    }

#else
    typename componentsConstructorTable::iterator cstrIter =
        componentsConstructorTablePtr_->find(methodName);

    if (cstrIter == componentsConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown newAMIMethod type "
            << methodName << nl << nl
            << "Valid newAMIMethod types are:" << nl
            << componentsConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }
    auto* ctorPtr = cstrIter();
#endif

    return autoPtr<newAMIMethod<SourcePatch, TargetPatch>>
    (
        ctorPtr
        (
            srcPatch,
            tgtPatch,
            srcMagSf,
            tgtMagSf,
            triMode,
            reverseTarget,
            requireMatch
        )
    );
}


#endif // end of #ifdef OPENFOAM_ORG

// ************************************************************************* //
