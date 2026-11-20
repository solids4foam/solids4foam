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
    Stop a run that continues from a written time whose time state
    (<time>/uniform/time) is missing.

    A run continued from a written time reads its time index and previous
    time-step from <time>/uniform/time. Without that file the time index
    starts from zero, and the models take the branches of a new run
    (first-step time schemes, restart fields not read, quantities recomputed
    on the deformed mesh) while the fields are those of a restart: the run
    continues with a plausible but different solution (solids4foam #740).
    This happens when a time is decomposed region by region
    (decomposePar -region), which does not copy <time>/uniform, or when the
    <time>/uniform that decomposePar links to (without -copyUniform) is
    deleted.

    The old-time fields written for a restart (*_0) show that the start time
    is a restart rather than a new run; a new run from time zero may hold
    old-time fields as initial conditions, so time zero is not checked. A new
    run started from another time with old-time fields as initial conditions
    is taken for a restart and needs a <time>/uniform/time (e.g. copied from a
    written time). The time directories are listed directly, so a case written
    with the collated file handler is not checked.

\*---------------------------------------------------------------------------*/

{
    // Whether this process starts from a restart time without its time state.
    // Evaluated on every process and reduced below, as the time state can be
    // missing on some processes only
    bool missingTimeState = false;

    if
    (
        runTime.startTimeIndex() == 0
     && mag(runTime.value()) > SMALL
     && !isFile(runTime.timePath()/"uniform"/"time")
    )
    {
        // Look for old-time fields in the start time and its region
        // directories. OpenFOAM.org spells the file type as a scoped
        // enumeration where the others use a member constant
#ifdef OPENFOAM_ORG
        fileNameList dirs(readDir(runTime.timePath(), fileType::directory));
#else
        fileNameList dirs(readDir(runTime.timePath(), fileName::DIRECTORY));
#endif
        dirs.append(fileName());

        forAll(dirs, dirI)
        {
#ifdef OPENFOAM_ORG
            const fileNameList files
            (
                readDir(runTime.timePath()/dirs[dirI], fileType::file)
            );
#else
            const fileNameList files
            (
                readDir(runTime.timePath()/dirs[dirI], fileName::FILE)
            );
#endif

            forAll(files, fileI)
            {
                const word name(files[fileI].name());

                if (name.size() > 2 && name.substr(name.size() - 2) == "_0")
                {
                    missingTimeState = true;
                    break;
                }
            }

            if (missingTimeState)
            {
                break;
            }
        }
    }

    if (returnReduce(missingTimeState, orOp<bool>()))
    {
        FatalErrorIn(args.executable())
            << "The run starts from time " << runTime.timeName()
            << ", which holds old-time fields of a restart, but its time "
            << "state " << runTime.timeName() << "/uniform/time is missing"
            << (Pstream::parRun() ? " in a processor directory." : ".") << nl
            << "    The time index would restart from zero and the models "
            << "would continue as a new run (first-step time schemes, "
            << "restart fields not read), which is not the continuation of "
            << "the original run." << nl
            << "    Restore " << runTime.timeName() << "/uniform from the "
            << "original run. For a parallel run, copy it into each "
            << "processor*/" << runTime.timeName() << "/ directory"
#ifdef OPENFOAM_COM
            << ", or remove the processor directories and decompose all "
            << "regions in one call, copying it:" << nl
            << "        decomposePar -regions '(fluid solid)' -copyUniform "
            << "-time " << runTime.timeName() << nl
            << "    rather than region by region"
#endif
            << "." << exit(FatalError);
    }
}

// ************************************************************************* //
