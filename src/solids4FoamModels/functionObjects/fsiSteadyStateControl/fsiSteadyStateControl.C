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

\*----------------------------------------------------------------------------*/

#include "fsiSteadyStateControl.H"
#include "addToRunTimeSelectionTable.H"
#include "pointFields.H"
#include "OSspecific.H"
#include "PstreamReduceOps.H"

// * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(fsiSteadyStateControl, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        fsiSteadyStateControl,
        dictionary
    );
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

bool Foam::fsiSteadyStateControl::sampleDisplacement
(
    scalar& displacement
) const
{
    if (!time_.foundObject<fvMesh>(solidRegionName_))
    {
        WarningInFunction
            << "Cannot find solid mesh region " << solidRegionName_ << nl;
        return false;
    }

    const fvMesh& mesh = time_.lookupObject<fvMesh>(solidRegionName_);

    if (!mesh.foundObject<pointVectorField>("pointD"))
    {
        WarningInFunction
            << "Cannot find pointVectorField pointD in region "
            << solidRegionName_ << nl;
        return false;
    }

    const pointVectorField& pointD =
        mesh.lookupObject<pointVectorField>("pointD");

    scalar minDist = GREAT;
    scalar localDisplacement = Zero;

    forAll(mesh.points(), pointI)
    {
        const scalar dist = mag(mesh.points()[pointI] - point_);

        if (dist < minDist)
        {
            minDist = dist;
            localDisplacement = pointD[pointI].component(displacementComponent_);
        }
    }

    if (Pstream::parRun())
    {
        Tuple2<scalar, scalar> localData(minDist, localDisplacement);

        List<Tuple2<scalar, scalar>> allData(Pstream::nProcs());
        allData[Pstream::myProcNo()] = localData;
        Pstream::allGatherList(allData);

        minDist = GREAT;
        forAll(allData, proci)
        {
            if (allData[proci].first() < minDist)
            {
                minDist = allData[proci].first();
                localDisplacement = allData[proci].second();
            }
        }
    }

    if (minDist == GREAT)
    {
        return false;
    }

    displacement = localDisplacement;
    return true;
}


bool Foam::fsiSteadyStateControl::sampleForce(scalar& force) const
{
    vector pressureForce(vector::zero);
    vector viscousForce(vector::zero);
    vector internalForce(vector::zero);

    if
    (
       !getObjectResult
        (
            forcesObjectName_,
            "pressureForce",
            pressureForce
        )
    )
    {
        WarningInFunction
            << "Cannot find pressureForce result from function object "
            << forcesObjectName_ << nl;
        return false;
    }

    if
    (
        includeViscousForce_
     && !getObjectResult
        (
            forcesObjectName_,
            "viscousForce",
            viscousForce
        )
    )
    {
        WarningInFunction
            << "Cannot find viscousForce result from function object "
            << forcesObjectName_ << nl;
        return false;
    }

    if
    (
        includeInternalForce_
     && !getObjectResult
        (
            forcesObjectName_,
            "internalForce",
            internalForce
        )
    )
    {
        WarningInFunction
            << "Cannot find internalForce result from function object "
            << forcesObjectName_ << nl;
        return false;
    }

    force =
    (
        pressureForce
      + (includeViscousForce_ ? viscousForce : vector::zero)
      + (includeInternalForce_ ? internalForce : vector::zero)
    ).component(forceComponent_);

    return true;
}


void Foam::fsiSteadyStateControl::openHistoryFile()
{
    if (!historyFilePtr_.empty())
    {
        return;
    }

    if (Pstream::master())
    {
        fileName historyDir;

        const word startTimeName =
            time_.timeName(time_.startTime().value());

        if (Pstream::parRun())
        {
            historyDir = time_.path()/".."/"postProcessing"/startTimeName;
        }
        else
        {
            historyDir = time_.path()/"postProcessing"/startTimeName;
        }

        mkDir(historyDir);

        historyFilePtr_.reset
        (
            new OFstream(historyDir/"fsiSteadyStateControl.dat")
        );

        if (historyFilePtr_.valid())
        {
            historyFilePtr_()
                << "# Time displacement force displacementAbsChange"
                << " displacementRelChange forceAbsChange forceRelChange"
                << " consecutive steady" << endl;
        }
    }
}


void Foam::fsiSteadyStateControl::writeHistory
(
    const scalar displacement,
    const scalar force,
    const scalar displacementAbsChange,
    const scalar displacementRelChange,
    const scalar forceAbsChange,
    const scalar forceRelChange,
    const bool steady
)
{
    if (Pstream::master() && historyFilePtr_.valid())
    {
        historyFilePtr_()
            << time_.timeOutputValue() << " "
            << displacement << " "
            << force << " "
            << displacementAbsChange << " "
            << displacementRelChange << " "
            << forceAbsChange << " "
            << forceRelChange << " "
            << consecutive_ << " "
            << steady << endl;
    }
}


bool Foam::fsiSteadyStateControl::checkSteadyState()
{
    if (stopped_)
    {
        return true;
    }

    scalar displacement = Zero;
    scalar force = Zero;

    if (!sampleDisplacement(displacement) || !sampleForce(force))
    {
        return false;
    }

    displacementHistory_.append(displacement);
    forceHistory_.append(force);

    scalar displacementAbsChange = -1;
    scalar displacementRelChange = -1;
    scalar forceAbsChange = -1;
    scalar forceRelChange = -1;

    bool steady = false;

    if
    (
        time_.timeOutputValue() >= minTime_
     && displacementHistory_.size() > nSamples_
    )
    {
        const label oldI = displacementHistory_.size() - nSamples_ - 1;

        displacementAbsChange =
            mag(displacement - displacementHistory_[oldI]);

        forceAbsChange =
            mag(force - forceHistory_[oldI]);

        displacementRelChange =
            displacementAbsChange
           /max(mag(displacementHistory_[oldI]), displacementRefMin_);

        forceRelChange =
            forceAbsChange/max(mag(forceHistory_[oldI]), forceRefMin_);

        steady =
        (
            displacementAbsChange <= displacementAbsTol_
         || displacementRelChange <= displacementRelTol_
        )
     && (
            forceAbsChange <= forceAbsTol_
         || forceRelChange <= forceRelTol_
        );
    }

    if (steady)
    {
        ++consecutive_;
    }
    else
    {
        consecutive_ = 0;
    }

    writeHistory
    (
        displacement,
        force,
        displacementAbsChange,
        displacementRelChange,
        forceAbsChange,
        forceRelChange,
        steady
    );

    if (log_)
    {
        Info<< type() << " " << name_ << ": displacement = "
            << displacement << ", force = " << force
            << ", displacementRelChange = " << displacementRelChange
            << ", forceRelChange = " << forceRelChange
            << ", consecutive = " << consecutive_ << endl;
    }

    setResult("displacement", displacement);
    setResult("force", force);
    setResult("displacementAbsChange", displacementAbsChange);
    setResult("displacementRelChange", displacementRelChange);
    setResult("forceAbsChange", forceAbsChange);
    setResult("forceRelChange", forceRelChange);
    setResult("steady", steady);

    if (consecutive_ >= nConsecutive_)
    {
        Info<< type() << " " << name_
            << ": steady state detected at time "
            << time_.timeOutputValue() << endl;

        stopped_ = true;

        if (writeAndEnd_)
        {
            Info<< "    Writing fields and stopping calculation" << endl;

            Time& runTime = const_cast<Time&>(time_);
            runTime.writeAndEnd();
            runTime.run();
        }
    }

    return true;
}


// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fsiSteadyStateControl::fsiSteadyStateControl
(
    const word& name,
    const Time& t,
    const dictionary& dict
)
:
    functionObjects::stateFunctionObject(name, t),
    name_(name),
    solidRegionName_("solid"),
    forcesObjectName_("forces"),
    point_(vector::zero),
    displacementComponent_(vector::Y),
    forceComponent_(vector::Y),
    minTime_(0),
    nSamples_(10),
    nConsecutive_(1),
    displacementAbsTol_(SMALL),
    displacementRelTol_(1e-5),
    forceAbsTol_(SMALL),
    forceRelTol_(1e-5),
    displacementRefMin_(SMALL),
    forceRefMin_(SMALL),
    includeViscousForce_(true),
    includeInternalForce_(false),
    writeAndEnd_(true),
    log_(true),
    consecutive_(0),
    stopped_(false),
    displacementHistory_(),
    forceHistory_(),
    historyFilePtr_()
{
    Info<< "Creating " << this->name() << " function object." << endl;

    read(dict);
    openHistoryFile();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::fsiSteadyStateControl::start()
{
    return false;
}


#if FOAMEXTEND
bool Foam::fsiSteadyStateControl::execute(const bool forceWrite)
#else
bool Foam::fsiSteadyStateControl::execute()
#endif
{
    return checkSteadyState();
}


bool Foam::fsiSteadyStateControl::read(const dictionary& dict)
{
    dict.readIfPresent("solidRegion", solidRegionName_);
    dict.readIfPresent("forcesObject", forcesObjectName_);
    dict.readIfPresent("point", point_);
    dict.readIfPresent("displacementComponent", displacementComponent_);
    dict.readIfPresent("forceComponent", forceComponent_);
    dict.readIfPresent("minTime", minTime_);
    dict.readIfPresent("nSamples", nSamples_);
    dict.readIfPresent("nConsecutive", nConsecutive_);
    dict.readIfPresent("displacementAbsTol", displacementAbsTol_);
    dict.readIfPresent("displacementRelTol", displacementRelTol_);
    dict.readIfPresent("forceAbsTol", forceAbsTol_);
    dict.readIfPresent("forceRelTol", forceRelTol_);
    dict.readIfPresent("displacementRefMin", displacementRefMin_);
    dict.readIfPresent("forceRefMin", forceRefMin_);
    dict.readIfPresent("includeViscousForce", includeViscousForce_);
    dict.readIfPresent("includeInternalForce", includeInternalForce_);
    dict.readIfPresent("writeAndEnd", writeAndEnd_);
    dict.readIfPresent("log", log_);

    if (nSamples_ < 1)
    {
        FatalIOErrorInFunction(dict)
            << "nSamples must be greater than zero"
            << exit(FatalIOError);
    }

    if (nConsecutive_ < 1)
    {
        FatalIOErrorInFunction(dict)
            << "nConsecutive must be greater than zero"
            << exit(FatalIOError);
    }

    return true;
}


#ifdef OPENFOAM_NOT_EXTEND
bool Foam::fsiSteadyStateControl::write()
{
    return true;
}
#endif

// ************************************************************************* //
