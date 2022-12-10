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

#include "fluidModel.H"
#include "volFields.H"
#include "fv.H"
#include "fvc.H"
#include "elasticWallPressureFvPatchScalarField.H"
#include "elasticSlipWallVelocityFvPatchVectorField.H"
#include "elasticWallVelocityFvPatchVectorField.H"
#include "EulerDdtScheme.H"
#include "backwardDdtScheme.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(fluidModel, 0);
    defineRunTimeSelectionTable(fluidModel, dictionary);
    addToRunTimeSelectionTable(physicsModel, fluidModel, physicsModel);
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fluidModel::makePisoControl() const
{
    if (!pisoPtr_.empty())
    {
        FatalErrorIn("void Foam::fluidModel::makePisoControl() const")
            << "pointer already set" << abort(FatalError);
    }

    pisoPtr_.set
    (
        new pisoControl
        (
            const_cast<fvMesh&>
            (
                refCast<const fvMesh>(mesh())
            )
        )
    );
}


void Foam::fluidModel::makePimpleControl() const
{
    if (!pimplePtr_.empty())
    {
        FatalErrorIn("void Foam::fluidModel::makePimpleControl() const")
            << "pointer already set" << abort(FatalError);
    }

    pimplePtr_.set
    (
        new pimpleControl
        (
            const_cast<fvMesh&>
            (
                refCast<const fvMesh>(mesh())
            )
        )
    );
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::fluidModel::updateRobinFsiInterface
(
    const volScalarField& p,
    const volVectorField& U,
    surfaceScalarField& phi,
    surfaceScalarField& rAUf
)
{
    forAll(p.boundaryField(), patchI)
    {
        if
        (
            (
                isA<elasticWallPressureFvPatchScalarField>
                (
                    p.boundaryField()[patchI]
                )
             && isA<elasticSlipWallVelocityFvPatchVectorField>
                (
                    U.boundaryField()[patchI]
                )
            )
         || (
                isA<elasticWallPressureFvPatchScalarField>
                (
                    p.boundaryField()[patchI]
                )
             && isA<elasticWallVelocityFvPatchVectorField>
                (
                    U.boundaryField()[patchI]
                )
            )
        )
        {
            const word ddtScheme =
#ifdef OPENFOAMESIORFOUNDATION
                word(mesh().ddtScheme("ddt(" + U.name() +')'));
#else
                mesh().schemesDict().ddtScheme("ddt(" + U.name() +')');
#endif

            if (ddtScheme == fv::EulerDdtScheme<vector>::typeName)
            {
#ifdef OPENFOAMESIORFOUNDATION
                phi.boundaryFieldRef()[patchI] =
                    phi.oldTime().boundaryField()[patchI];
                rAUf.boundaryFieldRef()[patchI] = runTime().deltaT().value();
#else
                phi.boundaryField()[patchI] =
                    phi.oldTime().boundaryField()[patchI];
                rAUf.boundaryField()[patchI] = runTime().deltaT().value();
#endif
            }
            else if (ddtScheme == fv::backwardDdtScheme<vector>::typeName)
            {
                if (runTime().timeIndex() == 1)
                {
#ifdef OPENFOAMESIORFOUNDATION
                    phi.boundaryFieldRef()[patchI] =
                        phi.oldTime().boundaryField()[patchI];
                    rAUf.boundaryFieldRef()[patchI] = runTime().deltaT().value();
#else
                    phi.boundaryField()[patchI] =
                        phi.oldTime().boundaryField()[patchI];
                    rAUf.boundaryField()[patchI] = runTime().deltaT().value();
#endif

                    phi.oldTime().oldTime();
                }
                else
                {
                    scalar deltaT = runTime().deltaT().value();
                    scalar deltaT0 = runTime().deltaT0().value();

                    scalar Cn = 1 + deltaT/(deltaT + deltaT0);
                    scalar Coo = deltaT*deltaT/(deltaT0*(deltaT + deltaT0));
                    scalar Co = Cn + Coo;

#ifdef OPENFOAMESIORFOUNDATION
                    phi.boundaryFieldRef()[patchI] =
                        (Co/Cn)*phi.oldTime().boundaryField()[patchI]
                      - (Coo/Cn)
                       *phi.oldTime().oldTime().boundaryField()[patchI];

                    rAUf.boundaryFieldRef()[patchI] =
                        runTime().deltaT().value()/Cn;
#else
                    phi.boundaryField()[patchI] =
                        (Co/Cn)*phi.oldTime().boundaryField()[patchI]
                      - (Coo/Cn)
                       *phi.oldTime().oldTime().boundaryField()[patchI];

                    rAUf.boundaryField()[patchI] =
                        runTime().deltaT().value()/Cn;
#endif
                }
            }
        }
    }
}


void Foam::fluidModel::CourantNo
(
    scalar& CoNum,
    scalar& meanCoNum,
    scalar& velMag
) const
{
    if (mesh().nInternalFaces())
    {
        surfaceScalarField magPhi(mag(phi()));

        if (phi().dimensions() == dimVelocity*dimArea*dimDensity)
        {
            const volScalarField& rho =
                mesh().lookupObject<volScalarField>("rho");

            magPhi /= fvc::interpolate(rho);
        }

        const surfaceScalarField SfUfbyDelta
        (
            mesh().surfaceInterpolation::deltaCoeffs()*magPhi
        );

        const scalar deltaT = runTime().deltaT().value();

        CoNum = max(SfUfbyDelta/mesh().magSf()).value()*deltaT;

        meanCoNum = (sum(SfUfbyDelta)/sum(mesh().magSf())).value()*deltaT;

        velMag = max(magPhi/mesh().magSf()).value();
    }

    Info<< "Courant Number mean: " << meanCoNum
        << " max: " << CoNum
        << " velocity magnitude: " << velMag
        << endl;
}


void Foam::fluidModel::CourantNo() const
{
    scalar CoNum = 0.0;
    scalar meanCoNum = 0.0;
    scalar velMag = 0.0;
    CourantNo(CoNum, meanCoNum, velMag);
}

#if FOAMEXTEND
void Foam::fluidModel::oversetCourantNo
(
    scalar& CoNum,
    scalar& meanCoNum,
    scalar& velMag
) const
{
    if (mesh().nInternalFaces())
    {
        const surfaceScalarField magPhi = mag(osMesh().sGamma()*phi());

        const surfaceScalarField SfUfbyDelta =
            mesh().surfaceInterpolation::deltaCoeffs()*magPhi;

        const scalar deltaT = runTime().deltaT().value();

        CoNum = max(SfUfbyDelta/mesh().magSf()).value()*deltaT;

        meanCoNum = (sum(SfUfbyDelta)/sum(mesh().magSf())).value()*deltaT;

        velMag = max(magPhi/mesh().magSf()).value();
    }

    Info<< "Courant Number mean: " << meanCoNum
        << " max: " << CoNum
        << " velocity magnitude: " << velMag
        << endl;
}


void Foam::fluidModel::oversetCourantNo() const
{
    scalar CoNum = 0.0;
    scalar meanCoNum = 0.0;
    scalar velMag = 0.0;
    oversetCourantNo(CoNum, meanCoNum, velMag);
}
#endif


void Foam::fluidModel::continuityErrs()
{
    const volScalarField contErr(fvc::div(phi()));

    const scalar sumLocalContErr = runTime().deltaT().value()*
        mag(contErr)().weightedAverage(mesh().V()).value();

    const scalar globalContErr = runTime().deltaT().value()*
        contErr.weightedAverage(mesh().V()).value();

    cumulativeContErr_ += globalContErr;

    Info<< "time step continuity errors : sum local = "
        << sumLocalContErr << ", global = " << globalContErr
        << ", cumulative = " << cumulativeContErr_
        << endl;
}

#if FOAMEXTEND
void Foam::fluidModel::oversetContinuityErrs()
{
    const volScalarField contErr = osMesh().gamma()*fvc::div(phi());

    const scalar sumLocalContErr = runTime().deltaT().value()*
        mag(contErr)().weightedAverage(mesh().V()).value();

    const scalar globalContErr = runTime().deltaT().value()*
        contErr.weightedAverage(mesh().V()).value();

    cumulativeContErr_ += globalContErr;

    Info<< "time step continuity errors : sum local = "
        << sumLocalContErr << ", global = " << globalContErr
        << ", cumulative = " << cumulativeContErr_
        << endl;
}
#endif

void Foam::fluidModel::boundPU
(
    volScalarField& p,
    volVectorField& U
) const
{
    // Bound the pressure
    dimensionedScalar p1 = min(p);
    dimensionedScalar p2 = max(p);

    if (p1 < pMin_ || p2 > pMax_)
    {
        Info<< "p: " << p1.value() << " " << p2.value()
            << ".  Bounding." << endl;

        p.max(pMin_);
        p.min(pMax_);
        p.correctBoundaryConditions();
    }

    // Bound the velocity
    volScalarField magU(mag(U));
    dimensionedScalar U1(max(magU));

    if (U1 > UMax_)
    {
        Info<< "U: " << U1.value() << ".  Bounding." << endl;

        volScalarField Ulimiter
        (
            pos(magU - UMax_)*UMax_/(magU + smallU_) + neg(magU - UMax_)
        );
        Ulimiter.max(scalar(0));
        Ulimiter.min(scalar(1));

        U *= Ulimiter;
        U.correctBoundaryConditions();
    }
}

#ifdef OPENFOAMESI
Foam::meshObjects::gravity Foam::fluidModel::readG() const
#else
Foam::uniformDimensionedVectorField Foam::fluidModel::readG() const
#endif
{
    // Note: READ_IF_PRESENT is incorreclty implemented within the
    // uniformDimensionedField constructor so we will use a work-around here

    // Check if waveProperties was read from disk; if so, then g must be read
    // from disk
    IOobject wavePropertiesHeader
    (
        "waveProperties",
        runTime().caseConstant(),
        mesh(),
        IOobject::MUST_READ,
        IOobject::NO_WRITE,
        false // do not register
    );

    // Check if g exists on disk
    IOobject gHeader
    (
        "g",
        runTime().caseConstant(),
        mesh(),
        IOobject::MUST_READ,
        IOobject::NO_WRITE,
        false // do not register
    );

#ifdef OPENFOAMESIORFOUNDATION
    if
    (
        wavePropertiesHeader.typeHeaderOk<IOdictionary>(true)
     && !gHeader.typeHeaderOk<uniformDimensionedVectorField>(true)
    )
#else
    if (wavePropertiesHeader.headerOk() && !gHeader.headerOk())
#endif
    {
        FatalErrorIn(type() + "::readG() const")
            << "g field not found in the constant directory: the g field "
            << "must be specified when the waveProperties dictionary is "
            << "present!" << abort(FatalError);
    }

    // The if-else-if structure is broken to keep the compiler happy with the
    // lack of a return statement above

#ifdef OPENFOAMESIORFOUNDATION
    if
    (
        wavePropertiesHeader.typeHeaderOk<IOdictionary>(true)
     || gHeader.typeHeaderOk<uniformDimensionedVectorField>(true)
    )
#else
    if (wavePropertiesHeader.headerOk() || gHeader.headerOk())
#endif
    {
        Info<< "Reading g from constant directory" << endl;
#ifdef OPENFOAMESI
    #if OPENFOAMESI > 1812
        return meshObjects::gravity(runTime());
    #else
        return meshObjects::gravity
        (
            runTime()
#if OPENFOAM < 1912
	    ,
            IOobject
            (
                "g",
                runTime().caseConstant(),
                mesh(),
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
#endif
        );
    #endif
#else
        return uniformDimensionedVectorField
        (
            IOobject
            (
                "g",
                runTime().caseConstant(),
                mesh(),
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        );
#endif
    }
    else
    {
        Info<< "g field not found in constant directory: initialising to zero"
            << endl;
#ifdef OPENFOAMESI
    #if OPENFOAMESI > 1812
        return meshObjects::gravity(runTime());
    #else
        return meshObjects::gravity
        (
            runTime()
#if OPENFOAM < 1912
	    ,
            IOobject
            (
                "g",
                runTime().caseConstant(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            )
#endif
        );
    #endif
#else
        return uniformDimensionedVectorField
        (
            IOobject
            (
                "g",
                runTime().caseConstant(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            dimensionedVector("zero", dimAcceleration, vector::zero)
        );
#endif
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fluidModel::fluidModel
(
    const word& type,
    Time& runTime,
    const word& region,
    const bool constructNull
)
:
    physicsModel(type, runTime),
    IOdictionary
    (
        IOobject
        (
            "fluidProperties",
            // If region == "region0" then read from the main case
            // Otherwise, read from the region/sub-mesh directory e.g.
            // constant/fluid or constant/solid
            bool(region == dynamicFvMesh::defaultRegion)
          ? fileName(runTime.caseConstant())
          : fileName(runTime.caseConstant()/region),
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    meshPtr_
    (
        dynamicFvMesh::New
        (
            IOobject
            (
                region,
                runTime.timeName(),
                runTime,
                IOobject::MUST_READ
            )
        )
    ),
    fluidProperties_(subDict(type + "Coeffs")),
    pisoPtr_(),
    pimplePtr_(),
    waveProperties_
    (
        IOobject
        (
            "waveProperties",
            runTime.constant(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        )
    ),
    g_(readG()),
    Uheader_("U", runTime.timeName(), mesh(), IOobject::MUST_READ),
    pheader_("p", runTime.timeName(), mesh(), IOobject::MUST_READ),
    UPtr_
    (
        constructNull
      ? nullptr
      : new volVectorField
        (
            IOobject
            (
                "U",
                runTime.timeName(),
                mesh(),
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            mesh(),
            dimensionedVector("zero", dimVelocity, vector::zero)
        )
    ),
    pPtr_
    (
        constructNull
      ? nullptr
      : new volScalarField
        (
            IOobject
            (
                "p",
                runTime.timeName(),
                mesh(),
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            mesh(),
            dimensionedScalar("zero", dimPressure, 0.0)
        )
    ),
    gradUPtr_
    (
        constructNull
      ? nullptr
      : new volTensorField
        (
            IOobject
            (
                "grad(U)",
                runTime.timeName(),
                mesh()
            ),
            mesh(),
            dimensionedTensor("zero", dimVelocity/dimLength, tensor::zero)
        )
    ),
    gradpPtr_
    (
        constructNull
      ? nullptr
      : new volVectorField
        (
            IOobject
            (
                "grad(p)",
                runTime.timeName(),
                mesh()
            ),
            mesh(),
            dimensionedVector("zero", p().dimensions()/dimLength, vector::zero)
        )
    ),
    phiPtr_
    (
        constructNull
      ? nullptr
      : new surfaceScalarField
        (
            IOobject
            (
                "phi",
                runTime.timeName(),
                mesh(),
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            fvc::interpolate(U()) & mesh().Sf()
        )
    ),
    adjustTimeStep_
    (
        runTime.controlDict().lookupOrDefault<Switch>("adjustTimeStep", false)
    ),
    maxCo_
    (
        runTime.controlDict().lookupOrDefault<scalar>("maxCo", 1.0)
    ),
    maxDeltaT_
    (
        runTime.controlDict().lookupOrDefault<scalar>("maxDeltaT", GREAT)
    ),
    pMin_("pMin", dimPressure, 0),
    pMax_("pMax", dimPressure, 0),
    UMax_("UMax", dimVelocity, 0),
    smallU_("smallU", dimVelocity, 1e-10),
    cumulativeContErr_(0.0),
#ifdef OPENFOAMFOUNDATION
    fvModels_(fvModels::New(mesh())),
    fvConstraints_(fvConstraints::New(mesh())),
#elif OPENFOAMESI
    fvOptions_(fv::options::New(mesh())),
#endif
    fsiMeshUpdate_(false),
    fsiMeshUpdateChanged_(false),
    globalPatchesPtrList_()
{
    if (!constructNull)
    {
        gradUPtr_() = fvc::grad(UPtr_());
        gradpPtr_() = fvc::grad(pPtr_());

        pMin_.dimensions().reset(p().dimensions());
        pMax_.dimensions().reset(p().dimensions());
        if (mesh().solutionDict().found("fieldBounds"))
        {
            dictionary fieldBounds = mesh().solutionDict().subDict("fieldBounds");
            fieldBounds.lookup(p().name())
                >> pMin_.value() >> pMax_.value();
            fieldBounds.lookup(U().name())
                >> UMax_.value();
        }
    }

#ifdef OPENFOAMFOUNDATION
    // Check if any finite volume models is present
    if (!fvModels_.PtrListDictionary<fvModel>::size())
    {
        Info << "No fvModels present" << endl;
    }

    // Check if any finite volume constrains is present
    if (!fvConstraints_.PtrListDictionary<fvConstraint>::size())
    {
        Info << "No fvConstraints present" << endl;
    }
#elif OPENFOAMESI
    // Check if any finite volume option is present
    if (!fvOptions_.optionList::size())
    {
        Info << "No finite volume options present\n" << endl;
    }
#endif
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fluidModel::~fluidModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::pisoControl& Foam::fluidModel::piso()
{
    if (pisoPtr_.empty())
    {
        makePisoControl();
    }

    return pisoPtr_();
}


Foam::pimpleControl& Foam::fluidModel::pimple()
{
    if (pimplePtr_.empty())
    {
        makePimpleControl();
    }

    return pimplePtr_();
}

#if FOAMEXTEND
const Foam::oversetMesh& Foam::fluidModel::osMesh() const
{
    return oversetMesh::New(mesh());
}
#endif

Foam::tmp<Foam::vectorField> Foam::fluidModel::faceZoneViscousForce
(
    const label interfaceI
) const
{
    const vectorField patchVF
    (
        patchViscousForce(globalPatches()[interfaceI].patch().index())
    );

    return globalPatches()[interfaceI].patchFaceToGlobal(patchVF);
}


Foam::tmp<Foam::scalarField> Foam::fluidModel::faceZonePressureForce
(
    const label interfaceI
) const
{
    const scalarField patchPF
    (
        patchPressureForce(globalPatches()[interfaceI].patch().index())
    );

    return globalPatches()[interfaceI].patchFaceToGlobal(patchPF);
}


void Foam::fluidModel::UisRequired()
{
#ifdef OPENFOAMESIORFOUNDATION
    if (!Uheader_.typeHeaderOk<volVectorField>(true))
#else
    if (!Uheader_.headerOk())
#endif
    {
        FatalErrorIn(type() + "::UisRequired()")
            << "This fluidModel requires the 'U' field to be specified!"
            << abort(FatalError);
    }
}


void Foam::fluidModel::pisRequired()
{
#ifdef OPENFOAMESIORFOUNDATION
    if (!pheader_.typeHeaderOk<volScalarField>(true))
#else
    if (!pheader_.headerOk())
#endif
    {
        FatalErrorIn(type() + "::pisRequired()")
            << "This fluidModel requires the 'p' field to be specified!"
            << abort(FatalError);
    }
}


void Foam::fluidModel::makeGlobalPatches(const wordList& patchNames) const
{
    globalPatchesPtrList_.setSize(patchNames.size());

    forAll(patchNames, i)
    {
        if (globalPatchesPtrList_.set(i))
        {
            FatalErrorIn
            (
                type() + "::makeGlobalPatches(const wordList&) const"
            )
                << "Pointer already set for global patch: "
                << patchNames[i] << "!"
                << abort(FatalError);
        }

        globalPatchesPtrList_.set
        (
            i,
            new globalPolyPatch(patchNames[i], mesh())
        );
    }
}


const Foam::PtrList<Foam::globalPolyPatch>&
Foam::fluidModel::globalPatches() const
{
    if (globalPatchesPtrList_.empty())
    {
        FatalErrorIn(type() + "::globalPatches() const")
            << "makeGlobalPatches(const wordList&) must be called "
            << "before globalPatch can be called!"
            << abort(FatalError);
    }

    return globalPatchesPtrList_;
}


void Foam::fluidModel::clearGlobalPatches() const
{
    globalPatchesPtrList_.clear();
}


void Foam::fluidModel::setDeltaT(Time& runTime)
{
    if (adjustTimeStep_)
    {
        // Calculate the maximum Courant number
        // Careful to use the relative flux in the calculation
        // We have to be careful when we call makeRelative and makeAbsolute
        scalar CoNum = 0.0;
        scalar meanCoNum = 0.0;
        scalar velMag = 0.0;
        fvc::makeRelative(phi(), U());
        CourantNo(CoNum, meanCoNum, velMag);
        fvc::makeAbsolute(phi(), U());

        scalar maxDeltaTFact = maxCo_/(CoNum + SMALL);
        scalar deltaTFact =
            min(min(maxDeltaTFact, 1.0 + 0.1*maxDeltaTFact), 1.2);

        runTime.setDeltaT
        (
            min
            (
                deltaTFact*runTime.deltaT().value(),
                maxDeltaT_
            )
        );

        Info<< "deltaT = " <<  runTime.deltaT().value() << endl;
    }
}


Foam::autoPtr<Foam::fluidModel> Foam::fluidModel::New
(
    Time& runTime,
    const word& region
)
{
    word fluidModelTypeName;

    // Enclose the creation of the dictionary to ensure it is
    // deleted before the fluid is created otherwise the dictionary
    // is entered in the database twice
    {
        IOdictionary fluidProperties
        (
            IOobject
            (
                "fluidProperties",
                 bool(region == dynamicFvMesh::defaultRegion)
              ? fileName(runTime.caseConstant())
              : fileName(runTime.caseConstant()/region),
                runTime,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        );

        fluidProperties.lookup("fluidModel")
            >> fluidModelTypeName;
    }

    Info<< nl << "Selecting fluidModel " << fluidModelTypeName << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(fluidModelTypeName);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "fluidModel::New(Time&, const word&)"
        )   << "Unknown fluidModel type " << fluidModelTypeName
            << endl << endl
            << "Valid fluidModel types are :" << endl
            << dictionaryConstructorTablePtr_->toc()
            << exit(FatalError);
    }

    return autoPtr<fluidModel>(cstrIter()(runTime, region));
}


void Foam::fluidModel::writeFields(const Time& runTime)
{
    physicsModel::writeFields(runTime);
}


bool Foam::fluidModel::read()
{
    if (regIOobject::read())
    {
        fluidProperties_ = subDict(type() + "Coeffs");

        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
