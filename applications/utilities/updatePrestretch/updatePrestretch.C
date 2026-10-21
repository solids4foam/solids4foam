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
    updatePrestretch

Description
    One pass of the modified updated Lagrangian formulation (Gee, Förster
    and Wall, 2010) for prestressing a solid model.

    After a prestressed run of a solid model (prestress yes; in
    mechanicalProperties), the deformation gradient F written at the latest
    time is the solver's own, from the mesh to the deformed configuration.
    The prestretch for the next pass is that deformation composed with the
    prestretch the run used:

        Fpre_new = F & Fpre

    F is on the left because it acts after Fpre. Copying F into Fpre
    instead discards the prestretch of every earlier pass, and the passes do
    not converge.

    Fpre_new is written to the start time, internal and patch values alike,
    at full precision, replacing the Fpre there, and every time directory
    after the start time is removed, so the next pass starts from the start
    time with zero displacement. Log files and postProcessing are left.

    The manager writes Fpre at every write time, so Fpre in the latest time
    says the run applied it. The prestretch is then taken from the start
    time, which holds it at full precision, and from the latest time only if
    the start time has none. Fpre in the start time but not the latest is
    refused, since the run did not apply it. On the first pass, where there
    is no Fpre in either, the prestretch the run used is taken to be the
    identity.

    max mag(F - I) is reported as the convergence measure: it tends to zero
    as the mesh comes into equilibrium with the load under the prestretch.

    A warning is given when the results are written in ascii with fewer than
    10 digits, which limits how well F, and so the prestretch, is known, and
    when the time scheme for D is not steadyState, since F is then one
    instant of a transient response rather than the equilibrium.

    A material whose law, or a law it wraps, is one of the linearElastic laws
    is refused before anything is written or removed: they are small-strain
    laws, which prestress does not apply to.

Usage
    updatePrestretch [-parallel]

Author
    Conor Moore, UCD.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "compatibilityFunctions.H"

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


//- Append to lawTypes the type of every law in dict, including the laws
//  wrapped by another law
void collectLawTypes(const dictionary& dict, DynamicList<word>& lawTypes)
{
    forAllConstIter(dictionary, dict, iter)
    {
        if (iter().isDict())
        {
            collectLawTypes(iter().dict(), lawTypes);
        }
        else if (iter().keyword() == "type")
        {
            lawTypes.append(word(iter().stream()));
        }
    }
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
#   include "setRootCase.H"
#   include "createTime.H"

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

    if (timeDirs.size() < 2)
    {
        FatalErrorInFunction
            << "There are no results after the start time to update the "
            << "prestretch from." << nl
            << "Run the prestressed solid model first"
            << exit(FatalError);
    }

    const instant startTime(timeDirs.first());
    const instant latestTime(timeDirs.last());

    runTime.setTime(latestTime, timeDirs.size() - 1);

#   include "createMesh.H"

    Info<< nl << "Updating the prestretch in time " << startTime.name()
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

        OStringStream smallStrainLaws;

        forAll(lawEntries, lawI)
        {
            DynamicList<word> lawTypes;
            collectLawTypes(lawEntries[lawI].dict(), lawTypes);

            forAll(lawTypes, typeI)
            {
                if (lawTypes[typeI].find("linearElastic") == 0)
                {
                    smallStrainLaws
                        << "    material " << lawEntries[lawI].keyword()
                        << " uses " << lawTypes[typeI] << nl;
                }
            }
        }

        if (!smallStrainLaws.str().empty())
        {
            FatalErrorInFunction
                << "Prestress applies only to finite-strain laws and solid "
                << "models, and these are small-strain laws:" << nl
                << smallStrainLaws.str().c_str()
                << "A run with them cannot apply a prestretch, and a "
                << "prestretch built from one is not one a finite-strain law "
                << "can use, so the prestretch is not updated"
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
            << IOstream::defaultPrecision() << ", so F, and the prestretch "
            << "built from it, are known only to that many digits." << nl
            << "    Set writePrecision to 10 or more, or writeFormat to "
            << "binary, in system/controlDict" << nl << endl;
    }

    const word ddtDScheme
    (
#ifdef OPENFOAM_NOT_EXTEND
        mesh.ddtScheme("ddt(D)")
#else
        mesh.schemesDict().ddtScheme("ddt(D)")
#endif
    );

    if (ddtDScheme != "steadyState")
    {
        WarningInFunction
            << "The time scheme for D is " << ddtDScheme << ", not "
            << "steadyState, so F at the latest time is one instant of a "
            << "transient response rather than the equilibrium under the "
            << "load, and a prestretch built from it does not converge."
            << nl
            << "    Use steadyState, with predictor no, for prestress passes"
            << nl << endl;
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

    // The prestretch the run used: in the latest time if the run applied it.
    // Unregistered, so that the updated prestretch can be registered under
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
        // The same prestretch, but the start time holds it as this utility
        // wrote it, at full precision, and the latest time as the run wrote
        // it, at writePrecision
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

    const volTensorField& FpreOld = *FpreOldPtr;

    // The updated prestretch, with the patch types of the old one
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
            << "det(Fpre) of the updated prestretch must be positive "
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

    // At full precision, so that the next pass starts from the prestretch
    // computed here rather than a rounding of it
    const unsigned int precision = IOstream::defaultPrecision();
    IOstream::defaultPrecision(17);

    Info<< "Writing Fpre to time " << startTime.name() << endl;
    FpreNew.write();

    IOstream::defaultPrecision(precision);

    // Nothing is removed unless the updated prestretch is on disk
    if (!isFile(FpreNew.objectPath()))
    {
        FatalErrorInFunction
            << "The updated prestretch was not written to "
            << FpreNew.objectPath() << ", so no results are removed"
            << exit(FatalError);
    }

    for (label timeI = 1; timeI < timeDirs.size(); timeI++)
    {
        Info<< "Removing time " << timeDirs[timeI].name() << endl;
        rmDir(runTime.path()/timeDirs[timeI].name());
    }

    Info<< nl << "End" << nl << endl;

    return(0);
}


// ************************************************************************* //
