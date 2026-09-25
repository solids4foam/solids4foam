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

#include "mechanicalConstitutiveLawInputGatherer.H"
#include "mechanicalConstitutiveLawGather.H"
#include "compatibilityFunctions.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::mechanicalConstitutiveLawInputGatherer::
mechanicalConstitutiveLawInputGatherer
(
    const fvMesh& mesh,
    const PtrList<mechanicalConstitutiveLaw>& laws,
    const List<labelList>& lawCells,
    const autoPtr<mechanicalConstitutiveLawCaseInputs>& caseInputsPtr,
    HashSet<word>& reportedUnreadInputs
)
:
    mesh_(mesh),
    laws_(laws),
    lawCells_(lawCells),
    caseInputsPtr_(caseInputsPtr),
    reportedUnreadInputs_(reportedUnreadInputs),
    diskScalarInputs_(),
    diskScalarInputsTimeIndex_(-1)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::volScalarField*
Foam::mechanicalConstitutiveLawInputGatherer::scalarInputSource
(
    const label lawI,
    const word& name
) const
{
    const bool sourced =
        caseInputsPtr_.valid() && caseInputsPtr_->found(lawI, name);

    if (mesh_.foundObject<volScalarField>(name))
    {
        const volScalarField& registered =
            mesh_.lookupObject<volScalarField>(name);

        if (!sourced)
        {
            return &registered;
        }

        if (!caseInputsPtr_->owns(registered))
        {
            // Solved for, or at least supplied, by another model, which takes
            // precedence
            caseInputsPtr_->reportShadowed(lawI, name);

            return &registered;
        }
    }

    if (sourced)
    {
        return &caseInputsPtr_->field(lawI, name);
    }

    return nullptr;
}


void Foam::mechanicalConstitutiveLawInputGatherer::refreshScalarInputs() const
{
    if (caseInputsPtr_.valid())
    {
        caseInputsPtr_->refreshAll();
    }

    const label timeIndex = mesh_.time().timeIndex();

    if (diskScalarInputsTimeIndex_ == timeIndex)
    {
        return;
    }

    diskScalarInputsTimeIndex_ = timeIndex;
    diskScalarInputs_.clear();

    // Every law and every input, in the same order on every rank. Only the
    // inputs that neither another model supplies nor a case directory does
    // are read: those are the ones a law would otherwise have read from the
    // file itself, on whichever ranks hold its points
    forAll(laws_, lawI)
    {
        const wordList names(requiredScalarInputsRecursive(laws_[lawI]));

        forAll(names, i)
        {
            const word& name = names[i];

            if (diskScalarInputs_.found(name) || scalarInputSource(lawI, name))
            {
                continue;
            }

            IOobject io
            (
                name,
                mesh_.time().timeName(),
                mesh_,
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                false
            );

            // An input given once in 0 is found on a restart too
            if (!mechanicalConstitutiveLawHeaderIsA<volScalarField>(io))
            {
                io.instance() = "0";
            }

            if (mechanicalConstitutiveLawHeaderIsA<volScalarField>(io))
            {
                diskScalarInputs_.insert
                (
                    name,
                    autoPtr<volScalarField>(new volScalarField(io, mesh_))
                );
            }
        }
    }
}


const Foam::volScalarField*
Foam::mechanicalConstitutiveLawInputGatherer::diskScalarInput
(
    const word& name
) const
{
    if (diskScalarInputsTimeIndex_ != mesh_.time().timeIndex())
    {
        FatalErrorInFunction
            << "The coupling input '" << name << "' was asked for before the "
            << "inputs were refreshed for this time step." << nl
            << "Every evaluation path must call refreshScalarInputs first."
            << abort(FatalError);
    }

    if (diskScalarInputs_.found(name))
    {
        return &autoPtrRef(diskScalarInputs_[name]);
    }

    return nullptr;
}


Foam::mechanicalConstitutiveLawInputs
Foam::mechanicalConstitutiveLawInputGatherer::lawInputsPatch
(
    const label lawI,
    const label patchI,
    const labelList& faces,
    const scalar dt,
    List<HashTable<autoPtr<scalarField>>>& lawScalarInputs
) const
{
    mechanicalConstitutiveLawInputs inputs(dt, mesh_.time().value());

    const wordList names(requiredScalarInputsRecursive(laws_[lawI]));

    if (names.empty())
    {
        return inputs;
    }

    if (lawScalarInputs.empty())
    {
        lawScalarInputs.setSize(laws_.size());
    }

    HashTable<autoPtr<scalarField>>& store = lawScalarInputs[lawI];

    forAll(names, i)
    {
        const word& name = names[i];

        // The boundary values are a different set from the internal ones, so
        // they get storage of their own rather than sharing it and being
        // overwritten by whichever was gathered last
        const word key(name + "@patch");

        if (!store.found(key))
        {
            store.insert(key, autoPtr<scalarField>(new scalarField()));
        }

        scalarField& fld = store[key]();
        fld.setSize(faces.size(), 0.0);

        // A registered or case-directory field is used where it is; only a
        // missing one is built. The two are kept apart rather than put
        // through one tmp, because foam-extend's tmp refuses to be assigned
        // one that holds a reference
        const volScalarField* srcPtr = scalarInputSource(lawI, name);

        if (!srcPtr)
        {
            srcPtr = diskScalarInput(name);
        }

        if (!srcPtr)
        {
            FatalErrorInFunction
                << "Mechanical constitutive law '" << laws_[lawI].type()
                << "' reads the coupling input '" << name
                << "', which is neither registered nor present as a field in "
                << "the current or initial time directory."
                << exit(FatalError);
        }

        const volScalarField& src = *srcPtr;

        const fvPatchField<scalar>& psrc = src.boundaryField()[patchI];

        if (psrc.coupled())
        {
            // A coupled face lies between two cells, as an internal face does,
            // and the surface paths give an internal face the interpolated
            // value. So does this, with the same weights, so that a face does
            // not see a different input for having become a processor face
            const scalarField& w = mesh_.weights().boundaryField()[patchI];
            const scalarField pif(psrc.patchInternalField());
            const scalarField pnf(psrc.patchNeighbourField());

            forAll(faces, faceI)
            {
                const label f = faces[faceI];
                fld[faceI] = w[f]*pif[f] + (1.0 - w[f])*pnf[f];
            }
        }
        else
        {
            forAll(faces, faceI)
            {
                fld[faceI] = psrc[faces[faceI]];
            }
        }

        inputs.setScalar(name, fld);
    }

    return inputs;
}


Foam::mechanicalConstitutiveLawInputs
Foam::mechanicalConstitutiveLawInputGatherer::lawInputs
(
    const label lawI,
    const integrationPointTopology& topo,
    const labelList& ipIDs,
    const scalar dt,
    List<HashTable<autoPtr<scalarField>>>& lawScalarInputs
) const
{
    mechanicalConstitutiveLawInputs inputs(dt, mesh_.time().value());

    reportUnreadScalarInputs(laws_[lawI]);

    const wordList names(requiredScalarInputsRecursive(laws_[lawI]));

    if (names.empty())
    {
        // Which is every law that does not couple to another field, so this
        // costs them nothing
        return inputs;
    }

    if (lawScalarInputs.empty())
    {
        lawScalarInputs.setSize(laws_.size());
    }

    HashTable<autoPtr<scalarField>>& store = lawScalarInputs[lawI];

    forAll(names, i)
    {
        const word& name = names[i];

        if (!store.found(name))
        {
            store.insert(name, autoPtr<scalarField>(new scalarField()));
        }

        scalarField& fld = store[name]();
        fld.setSize(ipIDs.size(), 0.0);

        // The registry first, and only then the file. This is the opposite
        // order from a prescribed field, and deliberately so: a coupling input
        // is solved for, so the live field must always win and a file must
        // never shadow it.
        //
        // The file is still needed, because of when the first evaluation
        // happens. A solid model builds its implicit stiffness while it is
        // being constructed, and a derived model that solves for the coupling
        // field has not reached its own members yet - the base class runs
        // first. At that moment the field exists only as the initial condition
        // on disk, which is the right value to evaluate against anyway
        //
        // A case that reads the field from another case directory, because
        // nothing in this run solves for it, is served in between: after any
        // registered field and before the file
        const volScalarField* srcPtr = scalarInputSource(lawI, name);

        if (srcPtr)
        {
            gatherToIntegrationPoints
            (
                mesh_, lawCells_[lawI], *srcPtr, topo, ipIDs, fld
            );
        }
        else
        {
            const volScalarField* diskPtr = diskScalarInput(name);

            if (!diskPtr)
            {
                FatalErrorInFunction
                    << "Mechanical constitutive law '" << laws_[lawI].type()
                    << "' reads the coupling input '" << name
                    << "', which is neither registered nor present as a field "
                    << "in the current or initial time directory." << nl
                    << "It is produced by another model, so the solid model "
                    << "has to solve for it, or the case has to supply it."
                    << exit(FatalError);
            }

            gatherToIntegrationPoints
            (
                mesh_, lawCells_[lawI], *diskPtr, topo, ipIDs, fld
            );
        }

        inputs.setScalar(name, fld);
    }

    return inputs;
}


Foam::wordList
Foam::mechanicalConstitutiveLawInputGatherer::requiredScalarInputsRecursive
(
    const mechanicalConstitutiveLaw& law
)
{
    // Deduplicated, because the same field may be wanted at more than one
    // level of the tree and is gathered once
    HashSet<word> names(law.requiredScalarInputs());

    const wordList childNames(law.childStateNames());

    forAll(childNames, i)
    {
        const wordList childInputs
        (
            requiredScalarInputsRecursive(law.childLaw(childNames[i]))
        );

        forAll(childInputs, j)
        {
            names.insert(childInputs[j]);
        }
    }

    return names.toc();
}


void Foam::mechanicalConstitutiveLawInputGatherer::reportUnreadScalarInputs
(
    const mechanicalConstitutiveLaw& law
) const
{
    const wordList names(law.optionalScalarInputs());

    forAll(names, i)
    {
        const word& name = names[i];

        if
        (
            !reportedUnreadInputs_.found(name)
         && mesh_.foundObject<volScalarField>(name)
        )
        {
            reportedUnreadInputs_.insert(name);

            WarningInFunction
                << "Mechanical constitutive law '" << law.type()
                << "' can read the registered field '" << name
                << "' as a coupling input, but has not been asked to, so it "
                << "is using its own substitute instead. Set the law's "
                << "option to read the field if the field should drive it."
                << endl;
        }
    }

    const wordList childNames(law.childStateNames());

    forAll(childNames, i)
    {
        reportUnreadScalarInputs(law.childLaw(childNames[i]));
    }
}


// ************************************************************************* //
