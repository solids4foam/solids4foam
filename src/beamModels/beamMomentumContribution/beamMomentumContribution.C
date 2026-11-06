/*---------------------------------------------------------------------------*\
License
    This file is part of beamFoam.

    beamFoam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    beamFoam is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with beamFoam.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "beamMomentumContribution.H"
#include "HermiteSpline.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(beamMomentumContribution, 0);
defineRunTimeSelectionTable(beamMomentumContribution, dictionary);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


autoPtr<beamMomentumContribution> beamMomentumContribution::New
(
    // const Time& runTime
    const word& beamMomentumContributionTypeName,
    const dictionary& dict
)
{
    auto* cstrIter = dictionaryConstructorTable
    (
        beamMomentumContributionTypeName
    );

    if (!cstrIter)
    {
        FatalIOErrorIn
        (
            "beamMomentumContribution::New(\n"
            "    const word& name,\n"
            "    const dictionary& dict\n"
            ")",
            dict
        )   << "Unknown beamMomentumContribution type "
            << beamMomentumContributionTypeName
            << endl << endl
            << "Valid beamMomentumContribution types are : " << endl
            << dictionaryConstructorTablePtr_->toc()
            << exit(FatalIOError);
    }

    return autoPtr<beamMomentumContribution>
    (
        cstrIter
        (
             beamMomentumContributionTypeName,
             dict
        )
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


beamMomentumContribution::beamMomentumContribution
(
    const word& name,
    const dictionary& dict
)
:
    name_(name)
{}


// * * * * * * * * * * * * * * * *  Destructors  * * * * * * * * * * * * * * //


beamMomentumContribution::~beamMomentumContribution()
{}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void beamMomentumContribution::checkSingleBeam(const beamModel& bm) const
{
    // The contributions use the spline and the radius of beam 0 alone, so a
    // region holding several beams would pair the other beams' cells with
    // out-of-range data
    if (bm.nBeams() > 1)
    {
        FatalErrorInFunction
            << "The " << type() << " beam momentum contribution supports"
            << " a single beam per region, but beam region '"
            << bm.solutionW().mesh().name() << "' holds " << bm.nBeams()
            << " beams (cellZones)." << nl
            << "See https://github.com/solids4foam/beamFoam/issues/46"
            << abort(FatalError);
    }
}


void beamMomentumContribution::checkCircularSection(const beamModel& bm) const
{
    const word sectionType(bm.crossSections()[0].type());

    if (sectionType != "circle")
    {
        FatalErrorInFunction
            << "The " << type() << " beam momentum contribution uses R()"
            << " as the radius of a circular cross-section, but beam region '"
            << bm.solutionW().mesh().name() << "' has cross-section type '"
            << sectionType << "'." << nl
            << "See https://github.com/solids4foam/beamFoam/issues/34"
            << abort(FatalError);
    }
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

tmp<vectorField> beamMomentumContribution::unitTangents
(
    const beamModel& bm
) const
{
    // Create spline using current beam points and tangents data
    // The 0 argument selects the 0th beam (beam index bI = 0).
    HermiteSpline spline
    (
        bm.currentBeamPoints(0),
        bm.currentBeamTangents(0)
    );

    // Evaluate dRdS - tangents to beam centreline at beam CV cell-centres.
    // HermiteSpline::firstDerivative already returns unit vectors; they are
    // normalised again so that the projections do not rely on it
    const vectorField& dRdScell = spline.midPointDerivatives();

    return tmp<vectorField>(new vectorField(dRdScell/(mag(dRdScell) + SMALL)));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void beamMomentumContribution::checkBeam(const beamModel& bm) const
{
    if (requiresSingleBeam())
    {
        checkSingleBeam(bm);
    }

    if (requiresCircularSection())
    {
        checkCircularSection(bm);
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
