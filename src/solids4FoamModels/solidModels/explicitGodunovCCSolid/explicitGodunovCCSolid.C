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
#include "fvm.H"
#include "fvc.H"
#include "fvMatrices.H"
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

void explicitGodunovCCSolid::updateStress()
{
    // Update increment of displacement
    DD() = D() - D().oldTime();

    // Update gradient of displacement
    mechanical().grad(D(), gradD());

    // Update gradient of displacement increment
    gradDD() = gradD() - gradD().oldTime();

    // Calculate the stress using run-time selectable mechanical law
    mechanical().correct(sigma());

    // Interpolate cell displacements to vertices
    mechanical().interpolate(D(), pointD());

    // Increment of displacement
    DD() = D() - D().oldTime();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

explicitGodunovCCSolid::explicitGodunovCCSolid
(
    Time& runTime,
    const word& region
)
:
    solidModel(typeName, runTime, region),
    impK_(mechanical().impK()),
    impKf_(mechanical().impKf()),
    rImpK_(1.0/impK_),
    //----------------------------------------------
    runTime_(runTime),
    //reading dicts
    mechanicalProperties_//! to be changed
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

    fvSolution_
    (
        IOobject
        (
            "fvSolution",
            runTime.system(),
            mesh(),
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),
    
    beta_
    (
        readScalar(fvSolution_.lookup("incompressiblilityCoefficient"))
    ),
   
    angularMomentumConservation_
    (
        fvSolution_.lookup("angularMomentumConservation")
    ),
    
    dampingCoeff_
    (
        "dampingCoeff",
        dimensionSet(0,0,-1, 0, 0, 0, 0),
        fvSolution_
    ),

//-----------------------------------------------------------------

    op_(mesh()),
    magSf_(mesh().magSf()),
    Sf_(mesh().Sf()),
    h_(op_.minimumEdgeLength()),
    bm_(mesh().boundaryMesh()),
    symmetricPatchID_(bm_.findPatchID("symmetric")),
    symmetricXpatchID_(bm_.findPatchID("symmetricX")),
    symmetricYpatchID_(bm_.findPatchID("symmetricY")),
    symmetricZpatchID_(bm_.findPatchID("symmetricZ")),

    // Creating mesh coordinate fields 
    C_(mesh().C()),
    
    x_
    (
        IOobject("x", mesh()),
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

    // Runge-Kutta stage
    RKstages_(2)

{

    Info << "Reading data from dictionaries ..." << endl;
    if
    (
        angularMomentumConservation_ != "yes" && angularMomentumConservation_ != "no"
    )
    {
        FatalErrorIn("readControls.H")
            << "Valid type entries are 'yes' or 'no' "
            << "for angularMomentumConservation"
            << abort(FatalError);
    }

    Info << "dampingCoeff: " << dampingCoeff_.value() << endl;


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

    DisRequired();

    // Update stress
    updateStress();
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
    

    // Mesh update loop
    do
    {
        if (physicsModel::printInfo())
        {
            Info<< "Evolving solid solver form explicitGodunovCCSolid" << endl;
        }

        lm_.oldTime();
        F_.oldTime();
        x_.oldTime();
        xF_.oldTime();
        xN_.oldTime();

        mech_.time(runTime_, deltaT_, max(Up_time_));

        // Info <<"deltaT_" << deltaT_<<endl;
        forAll(RKstages_, stage)
        {
            #include "gEqns.H"

            if (RKstages_[stage] == 0)
            {
                #include "updateVariables.H"
            }
        }

        lm_ = 0.5*(lm_.oldTime() + lm_);
        F_ = 0.5*(F_.oldTime() + F_);
        x_ = 0.5*(x_.oldTime() + x_);
        xF_ = 0.5*(xF_.oldTime() + xF_);
        xN_ = 0.5*(xN_.oldTime() + xN_);

        #include "updateVariables.H"


            if (runTime_.outputTime())
            {
                uN_ = xN_ - XN_;
                uN_.write();

                p_ = model_.pressure();
                p_.write();
                P_.write();
            }
            
            U() = lm_/rho_;
            
            // Compute displacement
            D() = D().oldTime() + deltaT_*U();

            // Enforce boundary conditions on the displacement field
            D().correctBoundaryConditions();

            // Update the stress field based on the latest D field
            updateStress();
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
    // Patch index
    const label patchID = patch.index();

    // Patch mechanical property
    const scalarField& impK = impK_.boundaryField()[patchID];

    // Patch reciprocal implicit stiffness field
    const scalarField& rImpK = rImpK_.boundaryField()[patchID];

    // Patch gradient
    const tensorField& pGradD = gradD().boundaryField()[patchID];

    // Patch stress
    const symmTensorField& pSigma = sigma().boundaryField()[patchID];

    // Patch unit normals
    const vectorField n(patch.nf());

    // Return patch snGrad
    return tmp<vectorField>
    (
        new vectorField
        (
            (
                (traction - n*pressure)
              - (n & (pSigma - impK*pGradD))
            )*rImpK
        )
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solidModels

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
