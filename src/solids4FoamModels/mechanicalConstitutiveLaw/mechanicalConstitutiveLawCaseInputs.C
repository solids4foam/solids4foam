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

#include "mechanicalConstitutiveLawCaseInputs.H"
#include "compatibilityFunctions.H"

// * * * * * * * * * * * * * * Static Functions  * * * * * * * * * * * * * * //

Foam::HashTable<Foam::fileName>
Foam::mechanicalConstitutiveLawCaseInputs::readSources
(
    const dictionary& lawDict,
    const word& lawName,
    const wordList& requiredInputs
)
{
    HashTable<fileName> sources;

    const word blockName("inputCaseDirectories");

    if (lawDict.found(blockName))
    {
        if (!lawDict.isDict(blockName))
        {
            FatalIOErrorInFunction(lawDict)
                << "'" << blockName << "' in material '" << lawName
                << "' must be a dictionary of <input> <caseDirectory> "
                << "entries." << exit(FatalIOError);
        }

        const dictionary& block = lawDict.subDict(blockName);
        const wordList names(block.toc());

        forAll(names, i)
        {
            if (findIndex(requiredInputs, names[i]) == -1)
            {
                FatalIOErrorInFunction(block)
                    << "Material '" << lawName << "' gives a case directory "
                    << "for input '" << names[i] << "', which its law does "
                    << "not read." << nl
                    << "    The inputs it reads are " << requiredInputs
                    << exit(FatalIOError);
            }

            sources.insert(names[i], fileName(block.lookup(names[i])));
        }
    }

    // The legacy spelling. It is what thermoMechanicalLaw's TcaseDirectory
    // is, for its input T, and it is read for every input rather than for
    // that one, so that the legacy dictionaries keep working without the
    // framework knowing which law they came from
    forAll(requiredInputs, i)
    {
        const word& name = requiredInputs[i];
        const word key(name + "caseDirectory");

        if (lawDict.found(key))
        {
            if (sources.found(name))
            {
                FatalIOErrorInFunction(lawDict)
                    << "Material '" << lawName << "' gives the case directory "
                    << "for input '" << name << "' twice: as '" << key
                    << "' and in '" << blockName << "'. Give one."
                    << exit(FatalIOError);
            }

            sources.insert(name, fileName(lawDict.lookup(key)));
        }
    }

    return sources;
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

void Foam::mechanicalConstitutiveLawCaseInputs::makeCase(sourceCase& sc) const
{
    // Relative to this case, as the legacy law's TcaseDirectory is. In
    // parallel caseName() is the processor directory, so the second case is
    // looked for inside it, again as the legacy law does
    const fileName casePath(mesh_.time().caseName()/sc.caseDir_);

    // Checked here so that a missing case is named as such, rather than
    // reported as a controlDict the user never meant to write
    const fileName controlDict
    (
        mesh_.time().rootPath()/casePath/"system"/Time::controlDictName
    );

    bool missing = !isFile(controlDict);
    reduce(missing, orOp<bool>());

    if (missing)
    {
        FatalErrorInFunction
            << "A mechanical constitutive law input is read from the case "
            << "directory " << sc.caseDir_ << ", but there is no "
            << controlDict
            << (Pstream::parRun() ? " on at least one processor." : ".")
            << nl
            << "    The directory is relative to the case. In parallel it is "
            << "relative to each processor directory, as on the legacy "
            << "thermoMechanicalLaw: processorN/" << sc.caseDir_
            << " must hold system/controlDict and that processor's portion "
            << "of the input case, decomposed as this case is."
            << exit(FatalError);
    }

    Info<< "Creating time for the mechanical constitutive law input case "
        << casePath << endl;

    sc.runTimePtr_.reset
    (
        new Time
        (
            Time::controlDictName,
            mesh_.time().rootPath(),
            casePath,
            "system",
            "constant",
            false        // No function objects: this case is only read
#ifdef OPENFOAM_COM
          , false        // No libraries either
#endif
        )
    );

    Time& runTime = autoPtrRef(sc.runTimePtr_);
    runTime.setTime(mesh_.time());

    Info<< "Reading the mesh of the input case " << casePath << endl;

    sc.meshPtr_.reset
    (
        new fvMesh
        (
            IOobject
            (
                fvMesh::defaultRegion,
                runTime.timeName(),
                runTime,
                IOobject::MUST_READ
            )
        )
    );

    checkMesh(sc, true);
}


void Foam::mechanicalConstitutiveLawCaseInputs::checkMesh
(
    const sourceCase& sc,
    const bool checkPoints
) const
{
    const fvMesh& srcMesh = sc.meshPtr_();

    // Reduced, so that every rank takes the same branch and a mismatch on one
    // rank stops them all instead of leaving the others waiting
    bool sizesDiffer =
        mesh_.nPoints() != srcMesh.nPoints()
     || mesh_.nFaces() != srcMesh.nFaces()
     || mesh_.nCells() != srcMesh.nCells()
     || mesh_.boundaryMesh().size() != srcMesh.boundaryMesh().size();

    if (!sizesDiffer)
    {
        forAll(mesh_.boundaryMesh(), patchI)
        {
            const polyPatch& p = mesh_.boundaryMesh()[patchI];
            const polyPatch& sp = srcMesh.boundaryMesh()[patchI];

            if (p.name() != sp.name() || p.size() != sp.size())
            {
                sizesDiffer = true;
            }
        }
    }

    reduce(sizesDiffer, orOp<bool>());

    if (sizesDiffer)
    {
        FatalErrorInFunction
            << "The mesh of the input case " << sc.caseDir_ << " does not "
            << "match the mesh of this case." << nl
            << "    Here: " << mesh_.nPoints() << " points, "
            << mesh_.nFaces() << " faces, " << mesh_.nCells() << " cells, "
            << "patches " << mesh_.boundaryMesh().names() << nl
            << "    There: " << srcMesh.nPoints() << " points, "
            << srcMesh.nFaces() << " faces, " << srcMesh.nCells()
            << " cells, patches " << srcMesh.boundaryMesh().names() << nl
            << "    The input is copied cell by cell and face by face, so the "
            << "two meshes must be the same. A mesh whose topology changed "
            << "during the run no longer matches the input case either."
            << exit(FatalError);
    }

    if (checkPoints)
    {
        if (gMax(mag(mesh_.points() - srcMesh.points())()) > SMALL)
        {
            FatalErrorInFunction
                << "The points of the mesh of the input case " << sc.caseDir_
                << " differ from those of this case." << exit(FatalError);
        }

        // Equal counts and points do not make equal numbering, and the input
        // is copied by index: the face addressing must match too, or a
        // renumbered source mesh would put each value in the wrong cell
        bool addressingDiffers =
            mesh_.faceOwner() != srcMesh.faceOwner()
         || mesh_.faceNeighbour() != srcMesh.faceNeighbour();

        forAll(mesh_.boundaryMesh(), patchI)
        {
            if
            (
                mesh_.boundaryMesh()[patchI].start()
             != srcMesh.boundaryMesh()[patchI].start()
            )
            {
                addressingDiffers = true;
            }
        }

        reduce(addressingDiffers, orOp<bool>());

        if (addressingDiffers)
        {
            FatalErrorInFunction
                << "The mesh of the input case " << sc.caseDir_ << " is "
                << "numbered differently from that of this case: the face "
                << "owner, neighbour or patch addressing differs." << nl
                << "    The input is copied cell by cell and face by face, so "
                << "the two meshes must be the same mesh."
                << exit(FatalError);
        }
    }
}


void Foam::mechanicalConstitutiveLawCaseInputs::prepareCase
(
    sourceCase& sc
) const
{
    const label timeIndex = mesh_.time().timeIndex();

    if (sc.timeIndex_ == timeIndex && sc.meshPtr_.valid())
    {
        return;
    }

    if (!sc.meshPtr_.valid())
    {
        makeCase(sc);
    }
    else
    {
        checkMesh(sc, false);
    }

    autoPtrRef(sc.runTimePtr_).setTime(mesh_.time());
    sc.timeIndex_ = timeIndex;
}


void Foam::mechanicalConstitutiveLawCaseInputs::refresh
(
    sourcedField& sf
) const
{
    const label timeIndex = mesh_.time().timeIndex();

    // Once per time step, as the legacy law does: the field is the input's
    // value over the step, not something that changes between iterations
    if (sf.timeIndex_ == timeIndex && sf.fieldPtr_.valid())
    {
        return;
    }

    sf.timeIndex_ = timeIndex;

    sourceCase& sc = cases_[sf.caseI_];
    prepareCase(sc);

    Time& runTime = autoPtrRef(sc.runTimePtr_);

    IOobject io
    (
        sf.name_,
        runTime.timeName(),
        sc.meshPtr_(),
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    );

#ifdef OPENFOAM_NOT_EXTEND
    bool present = io.typeHeaderOk<volScalarField>(true);
#else
    bool present =
        io.headerOk() && io.headerClassName() == volScalarField::typeName;
#endif

    // Reading can be collective under a collated file handler, so every rank
    // must either read or retain the preceding value together
    if (Pstream::parRun())
    {
        bool allPresent = present;
        reduce(allPresent, andOp<bool>());
        present = allPresent;
    }

    if (!present)
    {
        if (sf.fieldPtr_.valid())
        {
            // The last one read stands, which is how a case gives the field
            // once at the start and keeps it
            return;
        }

        FatalErrorInFunction
            << "The mechanical constitutive law input '" << sf.name_
            << "' is read from the case directory " << sc.caseDir_
            << ", but there is no volScalarField '" << sf.name_ << "' in "
            << runTime.timePath() << " and none was read before it." << nl
            << "    Give the field at least at the start time of the input "
            << "case."
            << exit(FatalError);
    }

    Info<< "Reading " << sf.name_ << " for the mechanical constitutive laws "
        << "from " << runTime.timePath() << endl;

    const volScalarField src(io, sc.meshPtr_());

    if (!sf.fieldPtr_.valid())
    {
        // Registered under the input's name, and written, as the legacy law's
        // copy is: it lets the case's own time directories show the field the
        // stress was computed with. Only if the name is free, though. A
        // second material reading the same input from another directory must
        // not displace the first one's copy, and an object this class did not
        // make is never replaced
        const bool registerCopy = !mesh_.thisDb().found(sf.name_);

        sf.fieldPtr_.reset
        (
            new volScalarField
            (
                IOobject
                (
                    sf.name_,
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    registerCopy ? IOobject::AUTO_WRITE : IOobject::NO_WRITE,
                    registerCopy
                ),
                mesh_,
                dimensionedScalar("zero", src.dimensions(), 0.0)
            )
        );
    }

    volScalarField& fld = autoPtrRef(sf.fieldPtr_);

    // Written into the current time directory, not the one it was made in
    fld.instance() = mesh_.time().timeName();

    primitiveFieldRef(fld) = primitiveField(src);

    forAll(src.boundaryField(), patchI)
    {
        boundaryFieldRef(fld)[patchI] =
            scalarField(src.boundaryField()[patchI]);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::mechanicalConstitutiveLawCaseInputs::mechanicalConstitutiveLawCaseInputs
(
    const fvMesh& mesh,
    const label nLaws
)
:
    mesh_(mesh),
    cases_(),
    fields_(),
    lawSources_(nLaws),
    reportedShadowed_()
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::mechanicalConstitutiveLawCaseInputs::~mechanicalConstitutiveLawCaseInputs()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::mechanicalConstitutiveLawCaseInputs::addSource
(
    const label lawI,
    const word& name,
    const fileName& caseDir
)
{
    // One Time and mesh per directory, however many inputs and laws read it
    label caseI = -1;

    forAll(cases_, i)
    {
        if (cases_[i].caseDir_ == caseDir)
        {
            caseI = i;
            break;
        }
    }

    if (caseI == -1)
    {
        caseI = cases_.size();
        cases_.setSize(caseI + 1);
        cases_.set(caseI, new sourceCase());
        cases_[caseI].caseDir_ = caseDir;
        cases_[caseI].timeIndex_ = -1;
    }

    // One copy per (directory, name), shared by the laws that read it
    label fieldI = -1;

    forAll(fields_, i)
    {
        if (fields_[i].caseI_ == caseI && fields_[i].name_ == name)
        {
            fieldI = i;
            break;
        }
    }

    if (fieldI == -1)
    {
        fieldI = fields_.size();
        fields_.setSize(fieldI + 1);
        fields_.set(fieldI, new sourcedField());
        fields_[fieldI].caseI_ = caseI;
        fields_[fieldI].name_ = name;
        fields_[fieldI].timeIndex_ = -1;
    }

    lawSources_[lawI].set(name, fieldI);

    Info<< "    Input " << name << " is read from the case directory "
        << caseDir << " at each time step" << endl;
}


bool Foam::mechanicalConstitutiveLawCaseInputs::found
(
    const label lawI,
    const word& name
) const
{
    return lawSources_[lawI].found(name);
}


bool Foam::mechanicalConstitutiveLawCaseInputs::owns
(
    const regIOobject& obj
) const
{
    forAll(fields_, i)
    {
        if
        (
            fields_[i].fieldPtr_.valid()
         && static_cast<const regIOobject*>(&fields_[i].fieldPtr_()) == &obj
        )
        {
            return true;
        }
    }

    return false;
}


void Foam::mechanicalConstitutiveLawCaseInputs::refreshAll() const
{
    forAll(fields_, fieldI)
    {
        refresh(fields_[fieldI]);
    }
}


const Foam::volScalarField&
Foam::mechanicalConstitutiveLawCaseInputs::field
(
    const label lawI,
    const word& name
) const
{
    if (!found(lawI, name))
    {
        FatalErrorInFunction
            << "Input '" << name << "' of law " << lawI << " is not read "
            << "from a case directory." << abort(FatalError);
    }

    sourcedField& sf = fields_[lawSources_[lawI][name]];

    refresh(sf);

    return sf.fieldPtr_();
}


void Foam::mechanicalConstitutiveLawCaseInputs::reportShadowed
(
    const label lawI,
    const word& name
) const
{
    if (reportedShadowed_.insert(name))
    {
        Info<< "The mechanical constitutive law input " << name << " is "
            << "registered by another model, which takes precedence over its "
            << "case directory "
            << cases_[fields_[lawSources_[lawI][name]].caseI_].caseDir_
            << endl;
    }
}


// ************************************************************************* //
