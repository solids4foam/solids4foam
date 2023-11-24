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

#include "AMIInterpolationS4F.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::autoPtr<Foam::AMIInterpolationS4F> Foam::AMIInterpolationS4F::New
(
    const word& modelName,
    const dictionary& dict,
    const bool reverseTarget
)
{
    DebugInfo << "Selecting model " << modelName << endl;

    auto cstrIter = dictConstructorTablePtr_->cfind(modelName);

    if (!cstrIter.found())
    {
        FatalErrorInLookup
        (
            typeName,
            modelName,
            *dictConstructorTablePtr_
        ) << exit(FatalError);
    }

    return autoPtr<AMIInterpolationS4F>(cstrIter()(dict, reverseTarget));
}


Foam::autoPtr<Foam::AMIInterpolationS4F> Foam::AMIInterpolationS4F::New
(
    const word& modelName,
    const bool requireMatch,
    const bool reverseTarget,
    const scalar lowWeightCorrection
)
{
    DebugInfo << "Selecting model " << modelName << endl;

    auto cstrIter = componentConstructorTablePtr_->cfind(modelName);

    if (!cstrIter.found())
    {
        FatalErrorInLookup
        (
            typeName,
            modelName,
            *componentConstructorTablePtr_
        ) << exit(FatalError);
    }

    return autoPtr<AMIInterpolationS4F>
    (
        cstrIter()
        (
            requireMatch,
            reverseTarget,
            lowWeightCorrection
        )
    );
}

// ************************************************************************* //
