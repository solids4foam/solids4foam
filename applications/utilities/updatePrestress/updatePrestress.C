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

Application
    updatePrestress

Description
    One pass of the modified updated Lagrangian formulation (Gee, Förster
    and Wall, 2010) for prestressing a solid model.

    After a prestressed run of a solid model (prestress yes; in
    mechanicalProperties), the deformation gradient F written at the latest
    time is the solver's own, from the mesh to the deformed configuration.
    Fpre for the next pass is that deformation composed with the
    Fpre the run used:

        FpreNew = F & Fpre

    F is on the left because it acts after Fpre. Copying F into Fpre
    instead discards the prestress of every earlier pass, and the passes do
    not converge.

    FpreNew is written to the start time, the startTime entry of
    system/controlDict, internal and patch values alike, at full precision,
    replacing the Fpre there, and every time directory after the start time
    is removed, so the next pass starts from the start time with zero
    displacement. Time directories before the start time, log files and
    postProcessing are left. A case written with a collated fileHandler is
    refused, since its time directories are not ones this can remove.

    The Fpre the start time held is kept, as it was, as FprePass<N>, N being
    the pass that wrote it, so that copying it back to Fpre restores it. Only
    that backup is kept, unless -keepPasses is given, in which case every
    earlier one is. Nothing is removed until the backup is made and FpreNew
    is written and read back as it is held, on every processor; if any of
    that fails, the Fpre that was there is put back and the utility stops.
    Each pass appends a line to prestressHistory in the case directory: the
    pass written, max mag(F - I) of the run that used the pass before it,
    and the range of det(Fpre).

    The manager writes the Fpre it used at every write time, at full
    precision, so Fpre in the latest time is the one the run used, and is
    the one composed with. If the start time also holds an Fpre it must be
    the same one; otherwise the run was continued from a time holding a
    different Fpre, or an earlier update did not finish, and the case is
    refused before anything is written or removed. Fpre in the start time
    but not the latest is refused, since the run did not apply it. On the
    first pass, where there is no Fpre in either, Fpre the run used is taken
    to be the identity.

    max mag(F - I) is reported as the convergence measure: it tends to zero
    as the mesh comes into equilibrium with the load under Fpre.

    A warning is given when the results are written in ascii with fewer than
    10 digits, which limits how well F, and so Fpre, is known.

    A case whose ddt or d2dt2 scheme for D is not steadyState is refused
    before anything is written or removed, unless -force is given: F is then
    one instant of a transient response rather than the equilibrium, and an
    Fpre built from it does not converge.

    A material whose law cannot be prestressed (see
    mechanicalConstitutiveLaw::admitsPrestress) is refused before anything is
    written or removed: the laws are constructed from mechanicalProperties
    and asked, as the mechanicalConstitutiveLawManager asks them, so the
    utility and the solver refuse the same laws.

Usage
    updatePrestress [-parallel] [-force] [-keepPasses]

Author
    Conor Moore, UCD.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "compatibilityFunctions.H"
#include "mechanicalConstitutiveLaw.H"
#include "Switch.H"
#include <fstream>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- True if the IOobject's file holds a volTensorField
bool foundVolTensorField(IOobject& io)
{
#ifdef OPENFOAM_NOT_EXTEND
    return io.typeHeaderOk<volTensorField>(true);
#else
    return io.headerOk();
#endif
}


//- The file a field is written to, compressed or not, or an empty name if
//  there is neither
fileName fieldFile(const fileName& path)
{
    if (isFile(path, false))
    {
        return path;
    }

    const fileName gzPath(path + ".gz");

    if (isFile(gzPath, false))
    {
        return gzPath;
    }

    return fileName();
}


//- The pass of a backup named FprePass<N>, or 0 for any other name
label backupPass(const fileName& name)
{
    const word base(name.ext() == "gz" ? name.lessExt().name() : name.name());
    const string prefix("FprePass");

    if
    (
        base.size() <= prefix.size()
     || base.compare(0, prefix.size(), prefix) != 0
    )
    {
        return 0;
    }

    const string digits(base.substr(prefix.size()));

    for (std::size_t charI = 0; charI < digits.size(); charI++)
    {
        if (!isdigit(digits[charI]))
        {
            return 0;
        }
    }

    return readLabel(IStringStream(digits)());
}


//- The backups, FprePass<N>, in the directory
fileNameList prestressBackups(const fileName& dir)
{
    const fileNameList files(readDir(dir, fileName::FILE));

    DynamicList<fileName> backups(files.size());

    forAll(files, fileI)
    {
        if (backupPass(files[fileI]) > 0)
        {
            backups.append(files[fileI]);
        }
    }

    return fileNameList(backups);
}


//- True if Fpre as read back from the instance is held: its internal values,
//  and those of every patch that writes them. A constraint patch, a symmetry
//  plane say, is evaluated from the internal values on reading rather than
//  read, so it is not compared. False, rather than fatal, if it cannot be
//  read, so that the caller can put back what was there
bool writtenAsHeld
(
    const volTensorField& held,
    const word& instance,
    const fvMesh& mesh
)
{
    FatalError.throwExceptions();
    FatalIOError.throwExceptions();

    bool same = false;

    try
    {
        const volTensorField read
        (
            IOobject
            (
                held.name(),
                instance,
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh
        );

        // At full precision or in binary the round trip is exact; the
        // tolerance only allows for the last digit
        const scalar tolerance =
            1e-14*max(1.0, max(mag(Foam::primitiveField(held))));

        same =
            max(mag(Foam::primitiveField(read) - Foam::primitiveField(held)))
         <= tolerance;

        forAll(read.boundaryField(), patchI)
        {
            const word& patchType = mesh.boundaryMesh()[patchI].type();

            if
            (
                polyPatch::constraintType(patchType)
             || read.boundaryField()[patchI].empty()
            )
            {
                continue;
            }

            same =
                same
             && max
                (
                    mag
                    (
                        read.boundaryField()[patchI]
                      - held.boundaryField()[patchI]
                    )
                )
             <= tolerance;
        }
    }
    catch (const Foam::error&)
    {
        same = false;
    }

    FatalError.dontThrowExceptions();
    FatalIOError.dontThrowExceptions();

    return same;
}


//- Set out to a & b, internal and boundary values alike
void compose
(
    volTensorField& out,
    const volTensorField& a,
    const volTensorField& b
)
{
    Foam::primitiveFieldRef(out) =
        Foam::primitiveField(a) & Foam::primitiveField(b);

    forAll(out.boundaryField(), patchI)
    {
        Foam::boundaryFieldRef(out)[patchI] ==
        (
            a.boundaryField()[patchI] & b.boundaryField()[patchI]
        );
    }
}


int main(int argc, char *argv[])
{
    argList::addBoolOption
    (
        "force",
        "Update Fpre even when the time schemes for D are not steadyState"
    );

    argList::addBoolOption
    (
        "keepPasses",
        "Keep every earlier Fpre as FprePass<N>, rather than only the one "
        "before this update"
    );

#   include "setRootCase.H"
#   include "createTime.H"

#ifdef OPENFOAM_COM
    const bool force = args.found("force");
    const bool keepPasses = args.found("keepPasses");
#else
    const bool force = args.optionFound("force");
    const bool keepPasses = args.optionFound("keepPasses");
#endif

#ifdef OPENFOAM_NOT_EXTEND
    // A collated case keeps every processor's time directories together,
    // and they are not ones rmDir on this processor's path can remove
    if (fileHandler().type() != "uncollated")
    {
        FatalErrorInFunction
            << "The fileHandler is " << fileHandler().type() << "; only an "
            << "uncollated case can be updated." << nl
            << "Write the case with '-fileHandler uncollated', or set "
            << "fileHandler uncollated in OptimisationSwitches"
            << exit(FatalError);
    }
#endif

    // The time directories, earliest first; constant is not one of them
    const instantList allTimes(runTime.times());
    DynamicList<instant> timeDirs(allTimes.size());

    forAll(allTimes, timeI)
    {
        if (allTimes[timeI].name() != runTime.constant())
        {
            timeDirs.append(allTimes[timeI]);
        }
    }

    if (timeDirs.empty())
    {
        FatalErrorInFunction
            << "There are no time directories" << exit(FatalError);
    }

    // The start time the passes start from, as the controlDict gives it,
    // rather than the earliest time directory, which may be an earlier
    // stage of the case. With startFrom firstTime a pass starts from the
    // earliest directory instead, so the two must agree; with latestTime it
    // starts from the start time, once the later directories are removed
    const word startFrom
    (
        runTime.controlDict().lookupOrDefault<word>("startFrom", "latestTime")
    );

    scalar startValue = timeDirs.first().value();

    if (runTime.controlDict().found("startTime"))
    {
        startValue = readScalar(runTime.controlDict().lookup("startTime"));

        if (startFrom == "firstTime" && !timeDirs.first().equal(startValue))
        {
            FatalErrorInFunction
                << "startFrom is firstTime, so a pass starts from time "
                << timeDirs.first().name() << ", but startTime is "
                << startValue << "." << nl
                << "Set startFrom startTime, or remove the time directories "
                << "before the start time, so that Fpre is written to the "
                << "time the next pass starts from"
                << exit(FatalError);
        }
    }
    else if (startFrom != "firstTime")
    {
        FatalErrorInFunction
            << "There is no startTime in system/controlDict to update Fpre "
            << "in, and startFrom is " << startFrom << exit(FatalError);
    }

    label startI = -1;

    forAll(timeDirs, timeI)
    {
        if (timeDirs[timeI].equal(startValue))
        {
            startI = timeI;
        }
    }

    if (startI == -1)
    {
        FatalErrorInFunction
            << "There is no time directory for the startTime " << startValue
            << " in system/controlDict" << exit(FatalError);
    }

    if (startI == timeDirs.size() - 1)
    {
        FatalErrorInFunction
            << "There are no results after the start time " << startValue
            << " to update the Fpre from." << nl
            << "Run the prestressed solid model first"
            << exit(FatalError);
    }

    const instant startTime(timeDirs[startI]);
    const instant latestTime(timeDirs.last());

    runTime.setTime(latestTime, timeDirs.size() - 1);

#   include "createMesh.H"

    Info<< nl << "Updating Fpre in time " << startTime.name()
        << " from the results at time " << latestTime.name() << nl << endl;

    {
        const IOdictionary mechanicalProperties
        (
            IOobject
            (
                "mechanicalProperties",
                runTime.constant(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                false
            )
        );

        const PtrList<entry> lawEntries
        (
            mechanicalProperties.lookup("mechanical")
        );

        // Each law is constructed as the manager constructs it, with the
        // settings the manager injects, and asked whether it can be
        // prestressed; this covers the laws a law wraps as well
        const Switch planeStress
        (
            mechanicalProperties.lookupOrDefault<Switch>("planeStress", false)
        );

        OStringStream refused;

        forAll(lawEntries, lawI)
        {
            dictionary lawDict(lawEntries[lawI].dict());

            // Given by the manager, and refused by it inside a material
            if (lawDict.found("planeStress") || lawDict.found("solutionD"))
            {
                FatalIOErrorInFunction(lawDict)
                    << "'planeStress' and 'solutionD' must not be given "
                    << "inside the '" << lawEntries[lawI].keyword()
                    << "' sub-dictionary: planeStress is set once in "
                    << "mechanicalProperties, and solutionD is supplied by "
                    << "the mesh" << exit(FatalIOError);
            }

            lawDict.add("planeStress", planeStress);
            lawDict.add("solutionD", mesh.solutionD());

            const autoPtr<mechanicalConstitutiveLaw> lawPtr
            (
                mechanicalConstitutiveLaw::New(lawDict)
            );

            if (!lawPtr->admitsPrestress())
            {
                refused
                    << "    material " << lawEntries[lawI].keyword()
                    << " uses " << lawPtr->type() << nl;
            }
        }

        if (!refused.str().empty())
        {
            FatalErrorInFunction
                << "These laws cannot be prestressed:" << nl
                << refused.str().c_str()
                << "A run with them cannot be prestressed, and an Fpre "
                << "built from one is not one they can use, so Fpre is not "
                << "updated"
                << exit(FatalError);
        }
    }

    if
    (
        runTime.writeFormat() == IOstream::ASCII
     && IOstream::defaultPrecision() < 10
    )
    {
        WarningInFunction
            << "The results are written in ascii with writePrecision "
            << IOstream::defaultPrecision() << ", so F, and Fpre "
            << "built from it, are known only to that many digits." << nl
            << "    Set writePrecision to 10 or more, or writeFormat to "
            << "binary, in system/controlDict" << nl << endl;
    }

    // Checked before anything is written or removed. A solid model takes its
    // inertia from d2dt2(D), and some its damping from ddt(D), so both have
    // to be steady for F at the latest time to be the equilibrium
    {
        const word ddtDScheme
        (
#ifdef OPENFOAM_NOT_EXTEND
            mesh.ddtScheme("ddt(D)")
#else
            mesh.schemesDict().ddtScheme("ddt(D)")
#endif
        );

        const word d2dt2DScheme(d2dt2SchemeCompat(mesh, "d2dt2(D)"));

        if (ddtDScheme != "steadyState" || d2dt2DScheme != "steadyState")
        {
            if (force)
            {
                WarningInFunction
                    << "The time schemes for D are " << ddtDScheme
                    << " (ddt) and " << d2dt2DScheme << " (d2dt2), not "
                    << "steadyState, so F at the latest time is one instant "
                    << "of a transient response. Updating Fpre anyway, as "
                    << "-force was given" << nl << endl;
            }
            else
            {
                FatalErrorInFunction
                    << "The time schemes for D are " << ddtDScheme
                    << " (ddt) and " << d2dt2DScheme << " (d2dt2), not "
                    << "steadyState, so F at the latest time is one instant "
                    << "of a transient response rather than the equilibrium "
                    << "under the load, and an Fpre built from it does not "
                    << "converge." << nl
                    << "Use steadyState for both, with predictor no, for "
                    << "prestress passes, or give -force to update Fpre "
                    << "anyway. Nothing has been written or removed"
                    << exit(FatalError);
            }
        }
    }

    const volTensorField F
    (
        IOobject
        (
            "F",
            latestTime.name(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        ),
        mesh
    );

    // Fpre the run used: in the latest time if the run applied it.
    // Unregistered, so that the updated Fpre can be registered under
    // the same name
    IOobject FpreLatestHeader
    (
        "Fpre",
        latestTime.name(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE,
        false
    );

    IOobject FpreStartHeader
    (
        "Fpre",
        startTime.name(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE,
        false
    );

    const bool FpreAtLatest = foundVolTensorField(FpreLatestHeader);
    const bool FpreAtStart = foundVolTensorField(FpreStartHeader);

    // Each processor looks at its own files, so the processors must agree on
    // which there are before any of them takes a branch with a reduction in
    // it, or those that do will wait for those that do not
    if
    (
        returnReduce(FpreAtLatest, orOp<bool>())
     != returnReduce(FpreAtLatest, andOp<bool>())
     || returnReduce(FpreAtStart, orOp<bool>())
     != returnReduce(FpreAtStart, andOp<bool>())
    )
    {
        FatalErrorInFunction
            << "Fpre is in time " << startTime.name() << " or time "
            << latestTime.name() << " on some processors but not on others."
            << nl
            << "Check the processor directories. Nothing has been written or "
            << "removed" << exit(FatalError);
    }

    autoPtr<volTensorField> FpreOldPtr;

    if (FpreAtLatest)
    {
        // The manager writes the Fpre it used at every write time, at full
        // precision, so the latest time's is the one the run used
        FpreOldPtr.reset(new volTensorField(FpreLatestHeader, mesh));

        // and the start time's, if there is one, has to be the same Fpre:
        // otherwise the run did not start from it, or an earlier update
        // wrote it without removing the results it was built from
        if (FpreAtStart)
        {
            const volTensorField FpreStart(FpreStartHeader, mesh);

            const scalar scale =
                max(1.0, gMax(mag(Foam::primitiveField(FpreStart))));

            const scalar diff = gMax
            (
                mag
                (
                    Foam::primitiveField(FpreStart)
                  - Foam::primitiveField(autoPtrRef(FpreOldPtr))
                )
            );

            if (diff > 1e-10*scale)
            {
                FatalErrorInFunction
                    << "Fpre in time " << startTime.name() << " is not the "
                    << "Fpre in time " << latestTime.name() << " (max mag "
                    << "difference " << diff << "), so the start time does "
                    << "not hold the Fpre the run used." << nl
                    << "Either the run was continued from a time holding a "
                    << "different Fpre, or an earlier update wrote Fpre to "
                    << "time " << startTime.name() << " without removing the "
                    << "results it was built from." << nl
                    << "Nothing has been written or removed"
                    << exit(FatalError);
            }
        }
    }
    else if (FpreAtStart)
    {
        FatalErrorInFunction
            << "Fpre is in time " << startTime.name() << " but not in time "
            << latestTime.name() << ", so the run did not apply it." << nl
            << "Check that mechanicalProperties has 'prestress yes;', spelt "
            << "so, and that the solver was built with prestress"
            << exit(FatalError);
    }
    else
    {
        Info<< "There is no Fpre: this is the first pass, and the run is "
            << "taken to have had none" << nl << endl;

        FpreOldPtr.reset
        (
            new volTensorField
            (
                IOobject
                (
                    "FpreOld",
                    startTime.name(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    false
                ),
                F
            )
        );

        autoPtrRef(FpreOldPtr) == dimensionedTensor("I", dimless, I);
    }

    const volTensorField& FpreOld = autoPtrRef(FpreOldPtr);

    // The updated Fpre, with the patch types of the old one
    volTensorField FpreNew
    (
        IOobject
        (
            "Fpre",
            startTime.name(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        FpreOld
    );

    compose(FpreNew, F, FpreOld);

    const volScalarField JpreNew("det(Fpre)", det(FpreNew));
    const scalar minJpre = min(JpreNew).value();

    if (minJpre <= SMALL)
    {
        FatalErrorInFunction
            << "det(Fpre) of the updated Fpre must be positive "
            << "everywhere, but its minimum is " << minJpre
            << exit(FatalError);
    }

    const volScalarField magFminusI
    (
        "mag(F - I)", mag(F - dimensionedTensor("I", dimless, I))
    );

    const scalar maxMagFminusI = max(magFminusI).value();
    const scalar maxJpre = max(JpreNew).value();

    Info<< "    max mag(F - I)    = " << maxMagFminusI << nl;

    Info<< "    det(Fpre) ranges from " << minJpre << " to " << maxJpre
        << nl << endl;

    // Writing moves an instance that is not constant or system to the
    // current time, so the current time is set back to the start time first
    runTime.setTime(startTime, 0);

    const fileName startDir(runTime.path()/startTime.name());

    // The backups are numbered by the pass that wrote them, so the Fpre in
    // the start time, if any, is from the pass after the latest backup. The
    // ranks agree on the number, so a backup has the same name on each
    label lastBackupPass = 0;
    const fileNameList backups(prestressBackups(startDir));

    forAll(backups, fileI)
    {
        lastBackupPass =
            max(lastBackupPass, backupPass(backups[fileI].name()));
    }

    reduce(lastBackupPass, maxOp<label>());

    const fileName oldFpreFile(fieldFile(startDir/"Fpre"));
    const bool backingUp = !oldFpreFile.empty();
    const label oldPass = backingUp ? lastBackupPass + 1 : lastBackupPass;
    const label newPass = oldPass + 1;

    // The Fpre in the start time is kept, as it was, as FprePass<N>, so that
    // copying it back to Fpre restores it. Nothing below removes a result
    // unless every rank has its backup made and the updated Fpre written
    // and read back as it is held here
    fileName backupFile;
    bool ok = true;

    if (backingUp)
    {
        backupFile = startDir/("FprePass" + Foam::name(oldPass));

        if (oldFpreFile.ext() == "gz")
        {
            backupFile += ".gz";
        }

        Info<< "Keeping the Fpre from pass " << oldPass << " as "
            << backupFile.name() << endl;

        ok = mv(oldFpreFile, backupFile);
    }

    if (ok)
    {
        // At full precision, so that the next pass starts from Fpre
        // computed here rather than a rounding of it
        const unsigned int precision = IOstream::defaultPrecision();
        IOstream::defaultPrecision(17);

        Info<< "Writing Fpre from pass " << newPass << " to time "
            << startTime.name() << endl;

        ok = FpreNew.write();

        IOstream::defaultPrecision(precision);
    }

    if (ok)
    {
        ok = writtenAsHeld(FpreNew, startTime.name(), mesh);
    }

    reduce(ok, andOp<bool>());

    if (!ok)
    {
        // Back as it was on every rank: the updated Fpre, if any of it was
        // written, is removed, and the backup is moved back
        const fileName newFpreFile(fieldFile(startDir/"Fpre"));

        if (!newFpreFile.empty())
        {
            rm(newFpreFile);
        }

        if (backingUp && isFile(backupFile, false))
        {
            mv(backupFile, oldFpreFile);
        }

        FatalErrorInFunction
            << "The updated Fpre could not be written to time "
            << startTime.name() << " and read back as it was computed, on "
            << "at least one processor." << nl
            << "The Fpre that was there has been put back, and no results "
            << "have been removed. The Fpre the run used is also in time "
            << latestTime.name() << exit(FatalError);
    }

    // Only the backup just made is kept, unless every pass is asked for
    if (!keepPasses)
    {
        forAll(backups, fileI)
        {
            const fileName backup(fieldFile(startDir/backups[fileI].name()));

            if (!backup.empty() && backup != backupFile)
            {
                Info<< "Removing " << backup.name() << endl;
                rm(backup);
            }
        }
    }

    // One line per pass in the case directory: the pass written, max mag(F -
    // I) of the run that used the pass before it, and the range of det(Fpre)
    if (Pstream::master())
    {
        const fileName historyFile
        (
            runTime.rootPath()/runTime.globalCaseName()/"prestressHistory"
        );

        const bool newHistory = !isFile(historyFile, false);

        std::ofstream history(historyFile.c_str(), std::ios::app);

        if (newHistory)
        {
            history
                << "# pass  maxMag(F-I)  minDet(Fpre)  maxDet(Fpre)" << nl
                << "# maxMag(F-I) is of the run that used the pass before"
                << nl;
        }

        history.precision(12);
        history
            << newPass << "  " << maxMagFminusI << "  " << minJpre << "  "
            << maxJpre << std::endl;

        if (!history.good())
        {
            WarningInFunction
                << "Could not append to " << historyFile << endl;
        }
    }

    for (label timeI = startI + 1; timeI < timeDirs.size(); timeI++)
    {
        Info<< "Removing time " << timeDirs[timeI].name() << endl;
        rmDir(runTime.path()/timeDirs[timeI].name());
    }

    Info<< nl << "End" << nl << endl;

    return(0);
}


// ************************************************************************* //
