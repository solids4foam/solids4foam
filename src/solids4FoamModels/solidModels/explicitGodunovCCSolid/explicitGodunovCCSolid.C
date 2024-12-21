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

#include "explicitGodunovCCSolid.H"
#include "addToRunTimeSelectionTable.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace solidModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(explicitGodunovCCSolid, 0);
addToRunTimeSelectionTable(solidModel, explicitGodunovCCSolid, dictionary);


// * * * * * * * * * * *  Private Member Functions * * * * * * * * * * * * * //

bool explicitGodunovCCSolid::converged
(
    const int iCorr,
    const dimensionedScalar pDeltaT,
    const GeometricField<vector, fvPatchField, volMesh>& vf
)
{
    // We will check three residuals:
    // - relative linear momentum residual

    bool converged = false;

    // Calculate residual based on the relative change of vf
    scalar denom = 0.0;

    // Denom is linear momentum increment
        denom = gMax
        (
#ifdef OPENFOAM_NOT_EXTEND
            DimensionedField<scalar, volMesh>
#else
            Field<scalar>
#endif
            (
                mag(vf.internalField() - vf.oldTime().internalField())
            )
        );

    if (denom < SMALL)
    {
        denom =
        max
        (
            gMax
            (
#ifdef OPENFOAM_NOT_EXTEND
                DimensionedField<scalar, volMesh>(mag(vf.internalField()))
#else
                mag(vf.internalField())
#endif
            ),
            SMALL
        );
    }

    const scalar residualvf =
        gMax
        (
#ifdef OPENFOAM_NOT_EXTEND
            DimensionedField<scalar, volMesh>
            (
                mag(vf.internalField() - vf.prevIter().internalField())
            )
#else
            mag(vf.internalField() - vf.prevIter().internalField())
#endif
        )/denom;

    // If one of the residuals has converged to an order of magnitude
    // less than the tolerance then consider the solution converged
    // force at least 1 outer iteration
   if (residualvf < solutionTol())
    {
        Info<< "    Converged" << endl;
        converged = true;
    }

    // Print residual information
    if (iCorr == 0)
    {
        Info<< "    Corr, res, pDeltaT" << endl;
    }
    else if (iCorr % infoFrequency() == 0 || converged || iCorr >= nCorr() - 1)
    {
        Info<< "    " << iCorr
            << ", " << residualvf
            << ", " << pDeltaT.value() << endl;

        if (iCorr >= nCorr())
        {
            Warning
                << "Max iterations reached within the momentum loop"
                << endl;
            converged = true;
        }
    }

    return converged;
 }


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

explicitGodunovCCSolid::explicitGodunovCCSolid
(
    Time& runTime,
    const word& region
)
:
    solidModel(typeName, runTime, region),
    runTime_(runTime),
    //reading dicts
    mechanicalProperties_
    (
        IOobject
        (
            "mechanicalProperties",
            runTime.constant(),
            mesh(),
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),

    controlDict_
    (
        IOobject
        (
            "controlDict",
            runTime.system(),
            mesh(),
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),

    incompressibilityCoefficient_
    (
        solidModelDict().lookupOrAddDefault<scalar>("incompressiblilityCoefficient", 1)
    ),

    beta_(incompressibilityCoefficient_),

    angularMomentumConservation_
    (
        solidModelDict().lookupOrAddDefault<word>("angularMomentumConservation", "no")
    ),

//-----------------------------------------------------------------

    op_(mesh()),
    magSf_(mesh().magSf()),
    Sf_(mesh().Sf()),
    h_(op_.minimumEdgeLength()),

    // Creating mesh coordinate fields
    C_(mesh().C()),

    x_
    (
        IOobject("x", mesh()),
        C_
    ),
    X_
    (
        IOobject("X", mesh()),
        C_
    ),

    xN_
    (
        IOobject("xN", mesh()),
        pMesh(),
        dimensionedVector("xN", dimensionSet(0,1,0,0,0,0,0), vector::zero)
    ),

    XN_(xN_),

    xF_(mesh().Cf()),

    // Creating mesh normal fields
    N_((Sf_ / mesh().magSf()).ref()),
    n_(N_),

    // Creating linear momentum fields
    lm_
    (
        IOobject
        (
            "lm",
            runTime.timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh(),
        dimensionedVector("lm", dimensionSet(1,-2,-1,0,0,0,0), vector::zero)
    ),
    lmN_
    (
        IOobject
        (
            "lmN",
            runTime.timeName(),
            mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        pMesh()
    ),

    F_
    (
        IOobject("F", mesh()),
        mesh(),
        tensor::I
    ),

    H_
    (
        IOobject("H", mesh()),
        det(F_)*op_.invT(F_)
    ),
    J_
    (
        IOobject("J", mesh()),
        det(F_)
    ),

    // Creating constitutive model
    model_
    (
        F_,
        mechanicalProperties_
    ),

    rho_(model_.density()),
    p_(model_.pressure()),
    P_(model_.piola()),
    Px_(op_.decomposeTensorX(P_)),
    Py_(op_.decomposeTensorY(P_)),
    Pz_(op_.decomposeTensorZ(P_)),

    mech_
    (
        F_,
        controlDict_
    ),

    Up_
    (
        IOobject("Up", mesh()),
        mesh(),
        model_.Up()/beta_
    ),

    Up_time_
    (
        IOobject("Up_time", mesh()),
        mesh(),
        model_.Up()
    ),

    Us_
    (
        IOobject("Us", mesh()),
        mesh(),
        model_.Us()*beta_
    ),

    grad_(mesh()),

    lmGrad_(grad_.gradient(lm_)),

    PxGrad_(grad_.gradient(Px_)),
    PyGrad_(grad_.gradient(Py_)),
    PzGrad_(grad_.gradient(Pz_)),

    // Reconstruction of linear momentum
    lm_M_(
        IOobject("lm_M", mesh()),
        mesh(),
        dimensionedVector("lm_M", lm_.dimensions(), vector::zero)
    ),

    lm_P_(lm_M_),

    // Reconstruction of PK1 stresses
    P_M_(
        IOobject("P_M", mesh()),
        mesh(),
        dimensionedTensor("P_M", P_.dimensions(), tensor::zero)
    ),
    P_P_(P_M_),

    // Reconstruction of traction
    t_M_
    (
        IOobject("t_M", mesh()),
        P_M_ & N_
    ),
    t_P_((P_P_ & N_).ref()),

    S_lm_(mech_.Smatrix_lm()),
    S_t_(mech_.Smatrix_t()),

    // Contact traction
    tC_(t_M_),
    // Contact linear momentum
    lmC_(lm_M_),
    t_b_
    (
        IOobject
        (
            "t_b",
            runTime.timeName(),
            mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh()
    ),

    lm_b_
    (
        IOobject
        (
            "lm_b",
            runTime.timeName(),
            mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh()
    ),

    // Constrained class
    interpolate_(mesh()),

    // Cell-averaged linear momentum
    lmR_(interpolate_.surfaceToVol(lmC_)),

    // Local gradient of cell-averaged linear momentum
    lmRgrad_(grad_.localGradient(lmR_, lmC_)),


    // Creating fields for angular momentum
    // Angular momentum class
    am_(mesh(), mechanicalProperties_),

    // RHS of linear momentum equation
    rhsLm_(
        IOobject("rhsLm", mesh()),
        mesh(),
        dimensionedVector("rhsLm", dimensionSet(1,-2,-2,0,0,0,0), vector::zero)
    ),
    rhsLm1_(rhsLm_),

    // RHS of angular momentum equation
    rhsAm_(
        IOobject("rhsAm", mesh()),
        mesh(),
        dimensionedVector("rhsAm", dimensionSet(1,-1,-2,0,0,0,0), vector::zero)
    ),
        // Nodal displacement field
    uN_(
        IOobject
        (
            "uN",
            runTime.timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        pMesh(),
        dimensionedVector("uN", dimLength, vector::zero)
    ),


    // Time increment
    deltaT_(
        "deltaT",
        dimTime,
        runTime.deltaTValue()
    ),
    // Time increment
    pDeltaT_(
        "pDeltaT",
        dimTime,
        runTime.deltaTValue()
    ),

    // Runge-Kutta stage
    RKstages_(2),

    phi_lm_
    (
        IOobject("phi_lm", mesh()),
        mesh(),
        dimensionedVector("phi_lm", dimensionSet(0,0,0,0,0,0,0), vector::zero)
    ),
    phi_P_
    (
        IOobject("phi_P", mesh()),
        mesh(),
        dimensionedTensor("phi_P", dimensionSet(0,0,0,0,0,0,0), tensor::zero)
    )


{

    Info << "Reading data from dictionaries ..." << endl;
    if
    (
        angularMomentumConservation_ != "yes" && angularMomentumConservation_ != "no"
    )
    {
        FatalErrorIn("angularMomentumConservation_")
            << "Valid type entries are 'yes' or 'no' "
            << "for angularMomentumConservation"
            << abort(FatalError);
    }
    Info << "Angular Momentum Conservation: " << angularMomentumConservation_ << endl;
    Info << "dampingCoeff: " << dampingCoeff().value() << endl;



    // Assign mesh points to the primitive field of xN_
    xN_.primitiveFieldRef() = mesh().points();
    XN_ = xN_;


    RKstages_[0] = 0;
    RKstages_[1] = 1;


    Info << "Printing data ..." << endl;

    // Print material properties
    model_.printMaterialProperties();

    // Print global linear and angular momentum
    am_.printGlobalMomentum(lm_,x_);

    // Print centroid of geometry
    mech_.printCentroid();
    #include "updateVariables.H"
    #include "riemannSolver.H"

    x_.oldTime();
    xF_.oldTime();
    lm_.oldTime();
    F_.oldTime();
    xN_.oldTime();

    lm_.oldTime().oldTime();
    F_.oldTime().oldTime();
    x_.oldTime().oldTime();
    xF_.oldTime().oldTime();
    xN_.oldTime().oldTime();

    // Set the printInfo
    physicsModel::printInfo() = bool
    (
        runTime.timeIndex() % infoFrequency() == 0
     || mag(runTime.value() - runTime.endTime().value()) < SMALL
    );

    Info<< "Frequency at which info is printed: every " << infoFrequency()
        << " time-steps" << endl;

}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


bool explicitGodunovCCSolid::evolve()
{
    Info<< "starting of evolve function" << endl;
    // Mesh update loop
    do
    {
        int iCorr = 0;

        if (physicsModel::printInfo())
        {
            Info<< "Evolving solid solver form explicitGodunovCCSolid" << endl;
        }

        // Pseudo time loop (Correction loop)
        do
        {
            if (angularMomentumConservation_ == "yes")
            {
                x_.storePrevIter();
                xF_.storePrevIter();
            }

            F_.storePrevIter();
            lm_.storePrevIter();
            xN_.storePrevIter();

            mech_.time(runTime_, pDeltaT_, max(Up_time_));

            forAll(RKstages_, stage)
            {
                #include "gEqns.H"

                if (RKstages_[stage] == 0)
                {
                    #include "updateVariables.H"
                }
            }

            if (angularMomentumConservation_ == "yes")
            {
                x_ = 0.5*(x_.prevIter() + x_);
                xF_ = 0.5*(xF_.prevIter() + xF_);
            }

            lm_ = 0.5*(lm_.prevIter() + lm_);
            F_ = 0.5*(F_.prevIter() + F_);
            xN_ = 0.5*(xN_.prevIter() + xN_);

            #include "updateVariables.H"

            pointD() = xN_ - XN_;

        }
        while
        (
            !converged
                (
                    iCorr,
                    pDeltaT_,
                    lm_
                )
         && ++iCorr < nCorr()
        );

        // Update the stress field based on the latest D field
        sigma() =  symm( (1.0 / J_) * (P_ & F_.T()));

        // Increment of point displacement
        pointDD() = pointD() - pointD().oldTime();

    }
    while (mesh().update());

    return true;
}


tmp<vectorField> explicitGodunovCCSolid::tractionBoundarySnGrad
(
    const vectorField& traction,
    const scalarField& pressure,
    const fvPatch& patch
) const
{
        Info << "tractionBoundarySnGrad is not implimented";
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solidModels

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
