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

        Fpre_new = F & Fpre

    F is on the left because it acts after Fpre. Copying F into Fpre
    instead discards the prestress of every earlier pass, and the passes do
    not converge.

    Fpre_new is written to the start time, the startTime entry of
    system/controlDict, internal and patch values alike, at full precision,
    replacing the Fpre there, and every time directory after the start time
    is removed, so the next pass starts from the start time with zero
    displacement. Time directories before the start time, log files and
    postProcessing are left. A case written with a collated fileHandler is
    refused, since its time directories are not ones this can remove.

    The manager writes Fpre at every write time, so Fpre in the latest time
    says the run applied it. Fpre is then taken from the start
    time, which holds it at full precision, and from the latest time only if
    the start time has none. Fpre in the start time but not the latest is
    refused, since the run did not apply it. On the first pass, where there
    is no Fpre in either, Fpre the run used is taken to be the
    identity.

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
    updatePrestress [-parallel] [-force]

Author
    Conor Moore, UCD.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "compatibilityFunctions.H"
#include "mechanicalConstitutiveLaw.H"
#include "Switch.H"

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

#   include "setRootCase.H"
#   include "createTime.H"

#ifdef OPENFOAM_COM
    const bool force = args.found("force");
#else
    const bool force = args.optionFound("force");
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

    // The start time the passes start from, as the controlDict gives it,
    // rather than the earliest time directory, which may be an earlier
    // stage of the case
    const scalar startValue
    (
        readScalar(runTime.controlDict().lookup("startTime"))
    );

    // The time directories, earliest first; constant is not one of them
    const instantList allTimes(runTime.times());
    DynamicList<instant> timeDirs(allTimes.size());
    label startI = -1;

    forAll(allTimes, timeI)
    {
        if (allTimes[timeI].name() != runTime.constant())
        {
            if (allTimes[timeI].equal(startValue))
            {
                startI = timeDirs.size();
            }

            timeDirs.append(allTimes[timeI]);
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

        const word d2dt2DScheme
        (
#ifdef OPENFOAM_NOT_EXTEND
            mesh.d2dt2Scheme("d2dt2(D)")
#else
            mesh.schemesDict().d2dt2Scheme("d2dt2(D)")
#endif
        );

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

    autoPtr<volTensorField> FpreOldPtr;

    if (FpreAtLatest && FpreAtStart)
    {
        // The same Fpre, both at full precision; the start time's is the
        // one this utility wrote, which a run started from an earlier
        // build, one that wrote Fpre at writePrecision, has not rounded
        FpreOldPtr.reset(new volTensorField(FpreStartHeader, mesh));
    }
    else if (FpreAtLatest)
    {
        FpreOldPtr.reset(new volTensorField(FpreLatestHeader, mesh));
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

    Info<< "    max mag(F - I)    = " << max(magFminusI).value() << nl;

    Info<< "    det(Fpre) ranges from " << minJpre << " to "
        << max(JpreNew).value() << nl << endl;

    // Writing moves an instance that is not constant or system to the
    // current time, so the current time is set back to the start time first
    runTime.setTime(startTime, 0);

    // At full precision, so that the next pass starts from Fpre
    // computed here rather than a rounding of it
    const unsigned int precision = IOstream::defaultPrecision();
    IOstream::defaultPrecision(17);

    Info<< "Writing Fpre to time " << startTime.name() << endl;
    FpreNew.write();

    IOstream::defaultPrecision(precision);

    // Nothing is removed unless the updated Fpre is on disk
    if (!isFile(FpreNew.objectPath()))
    {
        FatalErrorInFunction
            << "The updated Fpre was not written to "
            << FpreNew.objectPath() << ", so no results are removed"
            << exit(FatalError);
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
