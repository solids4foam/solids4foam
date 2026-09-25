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

#include "immersedBoundaryForce.H"
#include "fvMatrices.H"
#include "fvmSup.H"
#include "pimpleControl.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(immersedBoundaryForce, 0);
    addToRunTimeSelectionTable(option, immersedBoundaryForce, dictionary);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

const Foam::pimpleControl& Foam::fv::immersedBoundaryForce::pimple() const
{
    // The PIMPLE and PISO controls are registered as "solutionControl"
    const pimpleControl* pimplePtr =
        mesh_.findObject<pimpleControl>("solutionControl");

    if (!pimplePtr)
    {
        FatalErrorInFunction
            << "The " << typeName << " option " << name_
            << " requires a PIMPLE or PISO solution control on mesh "
            << mesh_.name() << abort(FatalError);
    }

    return *pimplePtr;
}


Foam::scalar Foam::fv::immersedBoundaryForce::rho()
{
    if (!rhoSet_)
    {
        rhoSet_ = true;

        const IOdictionary* transportPropertiesPtr =
            mesh_.findObject<IOdictionary>("transportProperties");

        if (transportPropertiesPtr && transportPropertiesPtr->found("rho"))
        {
            rho_ =
                dimensionedScalar
                (
                    "rho",
                    dimDensity,
                    *transportPropertiesPtr
                ).value();

            Info<< "    Immersed boundary " << name_ << ": density for the "
                << "forces " << rho_ << " from transportProperties" << endl;
        }
        else
        {
            Info<< "    Immersed boundary " << name_ << ": no density given, "
                << "the forces are per unit density" << endl;
        }
    }

    return rho_;
}


Foam::scalar Foam::fv::immersedBoundaryForce::nu()
{
    if (!nuSet_)
    {
        const IOdictionary* transportPropertiesPtr =
            mesh_.findObject<IOdictionary>("transportProperties");

        if (!transportPropertiesPtr || !transportPropertiesPtr->found("nu"))
        {
            FatalIOErrorInFunction(coeffs_)
                << "The " << typeName << " option " << name_ << " requires "
                << "the kinematic viscosity for the surface rate: give nu in "
                << "the option, or set surfaceRateCoeff 0" << exit(FatalIOError);
        }

        nu_ =
            dimensionedScalar
            (
                "nu",
                dimViscosity,
                *transportPropertiesPtr
            ).value();
        nuSet_ = true;
    }

    return nu_;
}


Foam::scalar Foam::fv::immersedBoundaryForce::cellWidth
(
    const label celli
) const
{
    const vector span
    (
        boundBox(mesh_.points(), mesh_.cellPoints()[celli], false).span()
    );

    const Vector<label>& solD = mesh_.solutionD();

    scalar w = GREAT;
    for (direction d = 0; d < vector::nComponents; ++d)
    {
        if (solD[d] == 1)
        {
            w = min(w, span[d]);
        }
    }

    return w;
}


void Foam::fv::immersedBoundaryForce::updateBodies(const bool move)
{
    lambda_.primitiveFieldRef() = 0;

    forAll(bodies_, bodyi)
    {
        immersedBody& body = bodies_[bodyi];

        if (move)
        {
            body.move(mesh_.time().value());
        }

        body.addOccupancy(lambda_, surfaceThreshold_);

        Info<< "    Immersed body " << body.name() << ": "
            << returnReduce(body.internalCells().size(), sumOp<label>())
            << " internal cells, "
            << returnReduce(body.surfaceCells().size(), sumOp<label>())
            << " surface cells" << endl;
    }

    lambda_.correctBoundaryConditions();

    // Record the mesh points for which lambda was calculated
    mesh_.setUpToDatePoints(pointsStamp_);
}


void Foam::fv::immersedBoundaryForce::calcForces(const volVectorField& U)
{
    const label timeIndex = mesh_.time().timeIndex();
    const scalar deltaT = mesh_.time().deltaTValue();
    const scalar rho = this->rho();

    // Store the fluid momentum of the previous time step, when the forces are
    // first calculated in this time step
    if (momentumTimeIndex_ != timeIndex)
    {
        momentum0Valid_ =
            momentumTimeIndex_ >= 0 && momentumTimeIndex_ == timeIndex - 1;
        momentum0_ = momentum_;
        momentumTimeIndex_ = timeIndex;
    }

    forAll(bodies_, bodyi)
    {
        const immersedBody& body = bodies_[bodyi];

        force_[bodyi] = body.force(f_, rho);
        torque_[bodyi] = body.torque(f_, rho);

        momentum_[bodyi] = body.momentum(U, lambda_);

        inertia_[bodyi] = Zero;
        if (momentum0Valid_)
        {
            inertia_[bodyi] =
                rho*(momentum_[bodyi] - momentum0_[bodyi])/deltaT;
        }

        Info<< "    Immersed body " << body.name() << ": force "
            << force_[bodyi] << ", torque " << torque_[bodyi]
            << ", inertia of the fluid inside " << inertia_[bodyi] << endl;
    }

    // In a partitioned fluid-solid interaction, the fluid may be solved
    // several times per time step: only the last solution is written
    forceTimeName_ = mesh_.time().timeName();
    forcesPending_ = true;
}


void Foam::fv::immersedBoundaryForce::writeForces()
{
    if (!forcesPending_)
    {
        return;
    }

    forcesPending_ = false;

    if (Pstream::master())
    {
        forAll(bodies_, bodyi)
        {
            OFstream& os = forceFiles_[bodyi];
            const vector& F = force_[bodyi];
            const vector& T = torque_[bodyi];
            const vector& Fi = inertia_[bodyi];

            os  << forceTimeName_ << token::TAB
                << F.x() << token::SPACE << F.y() << token::SPACE << F.z()
                << token::TAB
                << T.x() << token::SPACE << T.y() << token::SPACE << T.z()
                << token::TAB
                << Fi.x() << token::SPACE << Fi.y() << token::SPACE << Fi.z()
                << endl;
        }
    }
}


void Foam::fv::immersedBoundaryForce::setPenalty(const volVectorField& U)
{
    const scalar deltaT = mesh_.time().deltaTValue();
    const bool volumeFraction = (weighting_ == "volumeFraction");
    const scalarField& lambdaI = lambda_.primitiveField();

    Ui_ = U;

    for (const immersedBody& body : bodies_)
    {
        body.setVelocity(Ui_);
    }

    const scalar nu = (surfaceRateCoeff_ > 0 ? this->nu() : 0);

    kappa_ = 0;

    for (const immersedBody& body : bodies_)
    {
        for (const labelList* cellsPtr :
            {&body.internalCells(), &body.surfaceCells()})
        {
            for (const label celli : *cellsPtr)
            {
                const scalar lambda = lambdaI[celli];

                if (volumeFraction)
                {
                    // The penalised velocity is lambda*Ui + (1 - lambda)*U
                    kappa_[celli] =
                        min
                        (
                            penaltyCoeff_,
                            lambda/max(1 - lambda, 1/penaltyCoeff_)
                        )/deltaT;

                    // Rate independent of the time step in the partially
                    // covered cells
                    if (surfaceRateCoeff_ > 0 && lambda < 1 - surfaceThreshold_)
                    {
                        const scalar w = cellWidth(celli);

                        kappa_[celli] = min
                        (
                            kappa_[celli],
                            lambda/(1 - lambda)*surfaceRateCoeff_
                           *(nu/sqr(w) + mag(Ui_[celli])/w)
                        );
                    }
                }
                else
                {
                    kappa_[celli] = penaltyCoeff_*lambda/deltaT;
                }
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::immersedBoundaryForce::immersedBoundaryForce
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    fv::option(name, modelType, dict, mesh),
    bodies_(),
    method_("penalty"),
    penaltyCoeff_(1e3),
    weighting_("volumeFraction"),
    surfaceRateCoeff_(3),
    nu_(0),
    nuSet_(false),
    kappa_(mesh.nCells(), Zero),
    couplingCoeff_(0.8),
    surfaceThreshold_(1e-4),
    rho_(1),
    rhoSet_(false),
    writeSurfaces_(true),
    lambda_
    (
        IOobject
        (
            IOobject::scopedName(name, "lambda"),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar(dimless, Zero)
    ),
    Ui_
    (
        IOobject
        (
            IOobject::scopedName(name, "Ui"),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedVector(dimVelocity, Zero)
    ),
    f_
    (
        IOobject
        (
            IOobject::scopedName(name, "f"),
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedVector(dimVelocity/dimTime, Zero)
    ),
    f0_(mesh.nCells(), Zero),
    timeIndex_(-1),
    pointsStamp_
    (
        IOobject
        (
            IOobject::scopedName(name, "pointsStamp"),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        )
    ),
    UEqnPtr_(nullptr),
    UEqnTimeIndex_(-1),
    UEqnCorr_(-1),
    momentum_(),
    momentumTimeIndex_(-1),
    momentum0_(),
    momentum0Valid_(false),
    forceFiles_(),
    force_(),
    torque_(),
    inertia_(),
    forceTimeName_(),
    forcesPending_(false),
    warnedMatrix_(false)
{
    fieldNames_.resize(1, coeffs_.getOrDefault<word>("U", "U"));

    fv::option::resetApplied();

    Info<< "    Immersed bodies:" << endl;

    const dictionary& bodiesDict = coeffs_.subDict("bodies");

    label nBodies = 0;
    for (const entry& e : bodiesDict)
    {
        if (e.isDict())
        {
            ++nBodies;
        }
    }

    if (nBodies == 0)
    {
        FatalIOErrorInFunction(bodiesDict)
            << "No immersed bodies are given" << exit(FatalIOError);
    }

    bodies_.setSize(nBodies);

    label bodyi = 0;
    for (const entry& e : bodiesDict)
    {
        if (e.isDict())
        {
            bodies_.set
            (
                bodyi++,
                new immersedBody(e.keyword(), e.dict(), mesh)
            );
        }
    }

    momentum_.setSize(nBodies, Zero);
    momentum0_.setSize(nBodies, Zero);
    force_.setSize(nBodies, Zero);
    torque_.setSize(nBodies, Zero);
    inertia_.setSize(nBodies, Zero);

    read(dict);

    // Force files
    if (Pstream::master())
    {
        fileName dir(mesh.time().globalPath()/"postProcessing");

        if (mesh.name() != polyMesh::defaultRegion)
        {
            dir = dir/mesh.name();
        }

        dir = dir/name/mesh.time().timeName();

        mkDir(dir);

        forceFiles_.setSize(nBodies);

        forAll(bodies_, bodyi)
        {
            forceFiles_.set
            (
                bodyi,
                new OFstream(dir/(bodies_[bodyi].name() + ".dat"))
            );

            forceFiles_[bodyi].precision(10);

            forceFiles_[bodyi]
                << "# Immersed body " << bodies_[bodyi].name() << nl
                << "# Time" << token::TAB << "force (x y z)" << token::TAB
                << "torque (x y z)" << token::TAB
                << "inertia of the fluid inside (x y z)" << endl;
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fv::immersedBoundaryForce::~immersedBoundaryForce()
{
    writeForces();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fv::immersedBoundaryForce::addSup
(
    fvMatrix<vector>& eqn,
    const label fieldi
)
{
    const label timeIndex = mesh_.time().timeIndex();

    bool updated = false;

    if (timeIndex != timeIndex_)
    {
        // The previous time step is complete
        writeForces();

        // New time step: store the forcing of the previous time step (or the
        // forcing read on restart), and move the bodies to the new time
        f0_ = f_.primitiveField();

        bool moving = false;
        for (const immersedBody& body : bodies_)
        {
            moving = moving || body.moving();
        }

        if (timeIndex_ < 0 || moving || !mesh_.upToDatePoints(pointsStamp_))
        {
            updateBodies(true);
            updated = true;
        }

        timeIndex_ = timeIndex;

        if (writeSurfaces_ && mesh_.time().writeTime())
        {
            fileName dir(mesh_.time().globalPath()/"postProcessing");

            if (mesh_.name() != polyMesh::defaultRegion)
            {
                dir = dir/mesh_.name();
            }

            dir = dir/name_/"surfaces"/mesh_.time().timeName();

            for (const immersedBody& body : bodies_)
            {
                if (body.moving())
                {
                    body.writeSurface(dir/(body.name() + ".stl"));
                }
            }
        }
    }
    else if (!mesh_.upToDatePoints(pointsStamp_))
    {
        // The mesh has moved since lambda was calculated
        updateBodies(false);
        updated = true;
    }

    if (method_ == "penalty")
    {
        // Implicit volume penalisation: kappa*(Ui - U)
        const volVectorField& U = eqn.psi();

        setPenalty(U);

        volScalarField::Internal kappa
        (
            IOobject
            (
                IOobject::scopedName(name_, "kappa"),
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            dimensionedScalar(dimless/dimTime, Zero)
        );
        kappa.field() = kappa_;

        eqn += kappa*Ui_();
        eqn -= fvm::Sp(kappa, U);

        return;
    }

    if (pimple().firstIter())
    {
        // Start the time step, or a repeated solution of the time step, e.g.
        // in a partitioned fluid-solid interaction iteration, from the
        // forcing of the previous time step in the covered cells
        f_.primitiveFieldRef() = f0_*lambda_.primitiveField();
    }
    else if (updated)
    {
        // Remove the forcing from the cells no longer covered
        vectorField& fI = f_.primitiveFieldRef();
        const scalarField& lambdaI = lambda_.primitiveField();

        forAll(fI, celli)
        {
            if (lambdaI[celli] <= SMALL)
            {
                fI[celli] = Zero;
            }
        }
    }

    f_.correctBoundaryConditions();

    eqn += f_;
}


void Foam::fv::immersedBoundaryForce::addSup
(
    const volScalarField& rho,
    fvMatrix<vector>& eqn,
    const label fieldi
)
{
    FatalErrorInFunction
        << "The " << typeName << " option is only implemented for "
        << "incompressible flow" << abort(FatalError);
}


void Foam::fv::immersedBoundaryForce::constrain
(
    fvMatrix<vector>& eqn,
    const label fieldi
)
{
    UEqnPtr_ = &eqn;
    UEqnTimeIndex_ = mesh_.time().timeIndex();
    UEqnCorr_ = pimple().corr();
}


void Foam::fv::immersedBoundaryForce::correct(volVectorField& U)
{
    const pimpleControl& pimple = this->pimple();

    // Nothing to do after the momentum predictor
    if (pimple.corrPISO() < 1)
    {
        return;
    }

    if (method_ == "penalty")
    {
        // The forcing exerted by the penalisation, for the forces
        f_.primitiveFieldRef() =
            kappa_*(Ui_.primitiveField() - U.primitiveField());
        f_.correctBoundaryConditions();

        if (pimple.corrPISO() == pimple.nCorrPISO() && pimple.finalIter())
        {
            calcForces(U);
        }

        return;
    }

    // Velocity of the bodies in the covered cells, U elsewhere
    Ui_ = U;

    for (const immersedBody& body : bodies_)
    {
        body.setVelocity(Ui_);
    }

    // Forcing increment, in the cells with a nonzero occupancy
    volVectorField::Internal df
    (
        IOobject
        (
            IOobject::scopedName(name_, "df"),
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        mesh_,
        dimensionedVector(f_.dimensions(), Zero)
    );

    {
        const scalar coeff = couplingCoeff_/mesh_.time().deltaTValue();
        const scalarField& lambdaI = lambda_.primitiveField();
        const vectorField& UiI = Ui_.primitiveField();
        const vectorField& UI = U.primitiveField();
        vectorField& dfI = df.field();

        forAll(dfI, celli)
        {
            if (lambdaI[celli] > SMALL)
            {
                dfI[celli] = coeff*(UiI[celli] - UI[celli]);
            }
        }
    }

    f_.primitiveFieldRef() += df.field();
    f_.correctBoundaryConditions();

    if (pimple.corrPISO() < pimple.nCorrPISO())
    {
        // Add the increment to the source of the momentum matrix, so that the
        // next pressure corrector uses the updated forcing. The matrix is only
        // used if it is the one assembled for U in this outer corrector; it
        // is not accessed after the last corrector, as the solver may have
        // cleared it
        const bool matrixValid =
            UEqnPtr_
         && UEqnTimeIndex_ == mesh_.time().timeIndex()
         && UEqnCorr_ == pimple.corr()
         && &UEqnPtr_->psi() == &U;

        if (matrixValid)
        {
            // The forcing is on the right-hand side: eqn == f
            *UEqnPtr_ -= df;
        }
        else if (!warnedMatrix_)
        {
            warnedMatrix_ = true;

            WarningInFunction
                << "The momentum matrix of " << U.name() << " was not "
                << "constrained by the " << typeName << " option " << name_
                << " in this outer corrector, so the updated forcing is only "
                << "used from the next outer corrector" << endl;
        }
    }
    else
    {
        // Last pressure corrector
        UEqnPtr_ = nullptr;

        if (pimple.finalIter())
        {
            calcForces(U);
        }
    }
}


bool Foam::fv::immersedBoundaryForce::read(const dictionary& dict)
{
    if (fv::option::read(dict))
    {
        coeffs_.readIfPresent("method", method_);

        if (method_ != "incremental" && method_ != "penalty")
        {
            FatalIOErrorInFunction(coeffs_)
                << "Unknown method " << method_ << ": valid methods are "
                << "incremental and penalty" << exit(FatalIOError);
        }

        coeffs_.readCheckIfPresent
        (
            "penaltyCoeff",
            penaltyCoeff_,
            scalarMinMax::ge(1)
        );
        coeffs_.readIfPresent("weighting", weighting_);
        coeffs_.readIfPresent("surfaceRateCoeff", surfaceRateCoeff_);
        if (coeffs_.readIfPresent("nu", nu_))
        {
            nuSet_ = true;
        }

        if (weighting_ != "volumeFraction" && weighting_ != "occupancy")
        {
            FatalIOErrorInFunction(coeffs_)
                << "Unknown weighting " << weighting_ << ": valid weightings "
                << "are volumeFraction and occupancy" << exit(FatalIOError);
        }

        coeffs_.readCheckIfPresent
        (
            "couplingCoeff",
            couplingCoeff_,
            scalarMinMax::ge(SMALL)
        );

        // An occupancy threshold of 0.5 or more would classify every covered
        // cell as internal
        coeffs_.readCheckIfPresent
        (
            "surfaceThreshold",
            surfaceThreshold_,
            scalarMinMax(0, 0.5 - SMALL)
        );
        coeffs_.readIfPresent("writeSurfaces", writeSurfaces_);

        if (coeffs_.readIfPresent("rho", rho_))
        {
            rhoSet_ = true;
        }

        if (method_ == "penalty")
        {
            Info<< "    Penalty method: coefficient " << penaltyCoeff_
                << ", " << weighting_ << " weighting, surface rate "
                << "coefficient " << surfaceRateCoeff_;
        }
        else
        {
            Info<< "    Incremental method: coupling coefficient "
                << couplingCoeff_;
        }
        Info<< ", surface threshold " << surfaceThreshold_ << endl;

        return true;
    }

    return false;
}


// ************************************************************************* //
