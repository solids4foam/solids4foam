/*---------------------------------------------------------------------------*\
License
    This file is part of beamFoam.

    beamFoam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    beamFoam is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with beamFoam.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "coupledTotalLagNewtonRaphsonBeam.H"
#include "fvm.H"
#include "fvc.H"
#include "fvMatrices.H"
#include "addToRunTimeSelectionTable.H"
#include "HermiteSpline.H"
#include "spinTensor.H"
#include "pseudoVector.H"
// #include "momentBeamRotationFvPatchVectorField.H"
#include "momentBeamRotationNRFvPatchVectorField.H"
// #include "forceBeamDisplacementFvPatchVectorField.H"
#include "forceBeamDisplacementNRFvPatchVectorField.H"
#include "followerForceBeamDisplacementNRFvPatchVectorField.H"
// #include "axialForceTransverseDisplacementFvPatchVectorField.H"
// #include "axialForceTransverseDisplacementNRFvPatchVectorField.H"
// #include "extrapolatedBeamRotationFvPatchVectorField.H"

#include "mergePoints.H"
#include "scalarMatrices.H"
#include "denseMatrixHelperFunctions.H"
#include "BlockEigenSolverOF.H"
#include "samplingFluid.H"

// #include "beamHelperFunctions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace beamModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(coupledTotalLagNewtonRaphsonBeam, 0);
addToRunTimeSelectionTable
(
    beamModel,
    coupledTotalLagNewtonRaphsonBeam,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

coupledTotalLagNewtonRaphsonBeam::coupledTotalLagNewtonRaphsonBeam
(
    Time& runTime,
    const word& region
)
:
    beamModel(typeName, runTime, region),
    W_
    (
        IOobject
        (
            "W",
            runTime.timeName(),
            mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh()
    ),
    U_
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
        // fvc::ddt(W_)
        dimensionedVector("0", W_.dimensions()/dimTime, vector::zero)
    ),
    almForce_
    (
        IOobject
        (
            "almForce",
            mesh().time().timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh(),
        dimensionedVector("0", dimForce/dimLength, vector::zero)
    ),
    Accl_
    (
        IOobject
        (
            "Accl",
            runTime.timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh(),
        dimensionedVector("0", U_.dimensions()/dimTime, vector::zero)
    ),
    DW_
    (
        IOobject
        (
            "DW",
            runTime.timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedVector("0", W_.dimensions(), vector::zero)
    ),
    Theta_
    (
        IOobject
        (
            "Theta",
            runTime.timeName(),
            mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh()
    ),
    Omega_
    (
        IOobject
        (
            "Omega",
            runTime.timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh(),
        //fvc::ddt(Theta_)
        dimensionedVector("0", Theta_.dimensions()/dimTime, vector::zero)
    ),
    dotOmega_
    (
        IOobject
        (
            "dotOmega",
            runTime.timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh(),
        dimensionedVector("0", Omega_.dimensions()/dimTime, vector::zero)
    ),
    DTheta_
    (
        IOobject
        (
            "DTheta",
            runTime.timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedVector("0", Theta_.dimensions(), vector::zero)
    ),
    Q_
    (
        IOobject
        (
            "Q",
            runTime.timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh(),
        dimensionedVector("0", dimForce, vector::zero)
    ),
    explicitQ_
    (
        IOobject
        (
            "explicitQ",
            runTime.timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedVector("0", dimForce, vector::zero)
    ),
    M_
    (
        IOobject
        (
            "M",
            runTime.timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh(),
        dimensionedVector("M", dimForce*dimLength, vector::zero)
    ),
    explicitM_
    (
        IOobject
        (
            "explicitM",
            runTime.timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedVector("0", dimForce*dimLength, vector::zero)
    ),
    explicitMQ_
    (
        IOobject
        (
            "explicitMQ",
            runTime.timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedVector("0", dimForce*dimLength, vector::zero)
    ),
    Lambdaf_
    (
        IOobject
        (
            "Lambdaf",
            runTime.timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh(),
        dimensionedTensor("I", dimless, tensor::I)
    ),
    Lambda_
    (
        IOobject
        (
            "Lambda",
            runTime.timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh(),
        dimensionedTensor("I", dimless, tensor::I)
    ),
    refWf_
    (
        IOobject
        (
            "refWf",
            runTime.timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh(),
        dimensionedVector("zero", dimLength, vector::zero)
    ),
    refW_
    (
        IOobject
        (
            "refW",
            runTime.timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh(),
        dimensionedVector("zero", dimLength, vector::zero)
    ),
    refLambdaf_
    (
        IOobject
        (
            "refLambdaf",
            runTime.timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh(),
        dimensionedTensor("I", dimless, tensor::I)
    ),
    refLambda_
    (
        IOobject
        (
            "refLambda",
            runTime.timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh(),
        dimensionedTensor("I", dimless, tensor::I)
    ),
    refTangent_
    (
        IOobject
        (
            "refTangent",
            runTime.timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh(),
        dimensionedVector("x-axis", dimless, vector(1, 0, 0))
    ),
    CQW_
    (
        IOobject
        (
            "CQW",
            runTime.timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedTensor("0", dimPressure*dimArea, tensor::zero)
    ),
    CQTheta_
    (
        IOobject
        (
            "CQTheta",
            runTime.timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedTensor("0", dimPressure*dimArea, tensor::zero)
    ),
    CMTheta_
    (
        IOobject
        (
            "CMTheta",
            runTime.timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedTensor("0", dimPressure*dimArea*dimArea, tensor::zero)
    ),
    CMTheta2_
    (
        IOobject
        (
            "CMTheta2",
            runTime.timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedTensor("0", M_.dimensions(), tensor::zero)
    ),
    CMQW_
    (
        IOobject
        (
            "CMQW",
            runTime.timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedTensor("0", dimForce*dimLength, tensor::zero)
    ),
    CMQTheta_
    (
        IOobject
        (
            "CMQTheta",
            runTime.timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedTensor("0", dimForce*dimLength, tensor::zero)
    ),
    dR0Ds_
    (
        IOobject
        (
            "dR0Ds",
            runTime.timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedVector("0", dimless, vector::zero)
    ),
    Gamma_
    (
        IOobject
        (
            "Gamma",
            runTime.timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh(),
        dimensionedVector("0", dimless, vector::zero)
    ),
    K_
    (
        IOobject
        (
            "K",
            runTime.timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh(),
        dimensionedVector("0", dimless/dimLength, vector::zero)
    ),
    pMesh_(pointMesh::New(mesh())),
    pointW_
    (
        IOobject
        (
            "pointW",
            runTime.timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        pMesh_,
        dimensionedVector("0", dimLength, vector::zero)
    ),
    // WTheta_
    // (
    //     IOobject
    //     (
    //         "WTheta",
    //         runTime.timeName(),
    //         mesh(),
    //         IOobject::NO_READ,
    //         IOobject::AUTO_WRITE
    //     ),
    //     mesh(),
    //     dimensionedVector6("zero", dimless, vector6(1e-6))
    // ),
    totalIter_(0),
    CQ_
    (
        IOobject
        (
            "CQ",
            runTime.timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        // mesh(),
        indicator(0)
       *dimensionedTensor
        (
            "CQ",
            EA(0).dimensions(),
            tensor
            (
                EA(0).value(), 0, 0,
                0, GA(0).value(), 0,
                0, 0, GA(0).value()
            )
        )
    ),
    CM_
    (
        IOobject
        (
            "CM",
            runTime.timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        // mesh(), //CM_
        indicator(0)
       *dimensionedTensor
        (
            "CM",
            EI().dimensions(),
            tensor
            (
                GJ().value(), 0,              0,
                0,            EIyy().value(), 0,
                0,            0,              EIzz().value()
            )
        )
    ),
    // ARho_
    // (
    //         "ARho",
    //         A().dimensions()*rho().dimensions(),
    //         scalar(ARho().value())
    // ),
    // CIRho_
    // (
    //     "CIRho",
    //     I().dimensions()*rho().dimensions(),
    //     tensor
    //     (
    //         kCI()*JRho().value(),  0,           0,
    //         0,            kCI()*IyyRho().value(), 0,
    //         0,            0,           kCI()*IzzRho().value()
    //     )
    // ),
    ARho_
    (
        IOobject
        (
            "ARho",
            runTime.timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        fvc::average
        (
            indicator(0)
           *dimensionedScalar
            (
                "ARho",
                ARho(0).dimensions(),
                scalar(ARho(0).value())
            )
        )
    ),
    CIRho_
    (
        IOobject
        (
            "CIRho",
            runTime.timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        fvc::average
        (
            indicator(0)
           *dimensionedTensor
            (
                "CIRho",
                I().dimensions()*rho().dimensions(),
                tensor
                (
                    kCI()*JRho(0).value(),  0,                      0,
                    0,                       kCI()*IyyRho(0).value(), 0,
                    0,                       0,                      kCI()*IzzRho(0).value()
                )
            )
        )
    ),
    ddtSchemeName_
    (
         mesh().ddtSchemes().found("ddt(W)")
         ?
         (mesh().ddtSchemes().lookup("ddt(W)"))
         :
         (mesh().ddtSchemes().lookup("default"))
    ),
    d2dt2SchemeName_
    (
         mesh().d2dt2Schemes().found("d2dt2(W)")
         ?
         (mesh().d2dt2Schemes().lookup("d2dt2(W)"))
         :
         (mesh().d2dt2Schemes().lookup("default"))
    ),
    // Newmark-beta time integration parameters
    betaN_(0.0),
    gammaN_(0.0),
    // Pointer to add beamMomentumContributions
    momentumContribPtr_(0),
    totalContactTime_(0),
    totalSolutionTime_()
{
    // Check to enforce user to specify the d2dt2Schemes
    // system/fvSchemes - Required for dynamic simulations

    if (!List<word>({"steadyState", "Euler", "Newmark"}).found(d2dt2SchemeName_))
    {
        FatalErrorIn
        (
            "Constructor of coupledTotalLagNewtonRaphsonBeam"
        )
            << "The beam solver consists of 2nd order PDEs " << nl
            << "Please specify the d2dt2Schemes in system/fvSchemes!" << nl
            << "The valid schemes/keywords are: steadyState, Euler, Newmark." << nl
            << "The Newmark or Newmark-beta time scheme is a second-order "
            << "accurate time scheme for solving 2nd order PDEs." << nl
            << abort(FatalError);
    }

    // Check to enforce that ddtSchemes and d2dt2Schemes are consistent
    bool consistentTimeSchemes =
        (ddtSchemeName_ == "Euler" && d2dt2SchemeName_ == "Euler")
     || (ddtSchemeName_ == "steadyState" && d2dt2SchemeName_ == "steadyState")
     || (ddtSchemeName_ == "Newmark" && d2dt2SchemeName_ == "Newmark");
        // Newmark uses only d2dt2 scheme but set both to Newmark

    if (!consistentTimeSchemes)
    {
        FatalErrorInFunction
          << "Inconsistent time integration schemes! Specified schemes are: \n"
          << "  ddtScheme:    " << ddtSchemeName_  << nl
          << "  d2dt2Scheme:  " << d2dt2SchemeName_ << nl << nl
          << "Allowed combinations of (ddtScheme, d2dt2Scheme) are:\n"
          << "  (Euler, Euler)\n"
          << "  (steadyState, steadyState)\n"
          << "  (Newmark, Newmark)\n" << nl
          << abort(FatalError);
    }

    if (d2dt2SchemeName_ == "Newmark")
    {
        Info<< "The Newmark-beta time-integration scheme is set as d2dt2Scheme.\n"
            << nl
            << "The time integration parameters of Newmark method can be "
            << "specified inside d2dt2Scheme of system/fvSchemes using "
            << "keywords 'newmarkBeta(W)' and 'newmarkGamma(W)'. \n"
            << "The default parameters set are, newmarkBeta(W) = 0.25 "
            << "& newmarkGamma(W) = 0.5\n"
            << "They represent constant acceleration method or trapezoidal rule"
            << ", and are generally good. \n"
            << endl;

        betaN_ = mesh().d2dt2Schemes().lookupOrDefault<scalar>("newmarkBeta(W)", 0.25);
        gammaN_ = mesh().d2dt2Schemes().lookupOrDefault<scalar>("newmarkGamma(W)", 0.5);

        if (betaN_ != 0.25 || gammaN_ != 0.5)
        {
            Info<< "The user specified Newmark time integration parameters are\n"
                << "newmarkBeta = " << betaN_ << " , "
                << "newmarkGamma = " << gammaN_
                << endl;
        }
    }

    W_.oldTime();
    U_.oldTime();
    Accl_.oldTime();
    Omega_.oldTime();
    dotOmega_.oldTime();

    Lambda_.oldTime();
    Gamma_.oldTime();
    K_.oldTime();
    Lambdaf_.oldTime();

    Q_.oldTime();
    M_.oldTime();

    // Correct properties for multibeam case
    label nBeams = mesh().cellZones().size();
    if (nBeams > 1)
    {
        // First beam is already initialized
        for (label i = 1; i < nBeams; i++)
        {
            CQ_ +=
                indicator(i)
               *dimensionedTensor
                (
                    "CQ",
                    EA(i).dimensions(),
                    tensor
                    (
                        EA(i).value(), 0, 0,
                        0, GA(i).value(), 0,
                        0, 0, GA(i).value()
                    )
                );

            CM_ +=
                indicator(i)
               *dimensionedTensor
                (
                    "CM",
                    EI(i).dimensions(),
                    tensor
                    (
                        GJ(i).value(), 0, 0,
                        0, EIyy(i).value(), 0,
                        0, 0, EIzz(i).value()
                    )
                );

            ARho_ +=
                fvc::average
                (
                    indicator(i)
                   *dimensionedScalar
                    (
                        "ARho",
                        ARho(i).dimensions(),
                        scalar(ARho(i).value())
                    )
                );

            CIRho_ +=
                fvc::average
                (
                    indicator(i)
                   *dimensionedTensor
                    (
                        "CIRho",
                        I().dimensions()*rho().dimensions(),
                        tensor
                        (
                            kCI()*JRho(i).value(),    0,                      0,
                            0,                        kCI()*IyyRho(i).value(), 0,
                            0,                        0,                      kCI()*IzzRho(i).value()
                        )
                    )
                );
        }
    }

    // Check whether the reference longitudinal axis of beams is
    // always set in global x-direction and throw error otherwise
    // This is a mandatory step because the system of beam equations
    // and the strains are defined in the code assuming the reference
    // beam centreline to be aligned in global x-direction
    volVectorField R0
    (
        IOobject
        (
            "R0",
            runTime.timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimLength
    );
    R0 = mesh().C();
    const vectorField refTangentError(fvc::snGrad(R0) - vector(1, 0, 0));

    if (max(mag(refTangentError)) > 1e-06)
    {
        FatalErrorIn
        (
            "Constructor of Foam::coupledTotalLagNewtonRaphsonBeam "
        )   << "The longitudinal axis of beam in the reference configuration " << nl
            << "is not aligned the global x-axis : This is a mandatory " << nl
            << "requirement before running the beam solver. Check the beam mesh " << nl
            << "in paraview. For example, " << nl
            << "if the current beam is oriented in the Z direction, run the "
            << "following command: " << nl
            << "    transformPoints -rotate-angle '((0 1 0) 90)'" << nl
            << "after creating the beam mesh to set the correct reference " << nl
            << nl
            << "You might also be using a deprecated mesh utility " << nl
            << " 'createCircularBeamMesh' or 'createMultipleBeamMesh' "
            << " along with transformPoints command! " << nl << nl
            << " Use 'createBeamMesh' utility instead! and remove  "
            << " transformPoints command! " << nl
            << abort(FatalError);
    }

    // Calculate tangents if it is not already set by the user
    // This means the initial configuration is equivalent to the
    // reference configuration and the tangents are set to vector(1,0,0).
    IOobject refTangentHeader
    (
        "refTangent",
        runTime.timeName(),
        mesh(),
        IOobject::MUST_READ
    );
    if (!refTangentHeader.typeHeaderOk<surfaceVectorField>(true))
    {
        Info<< "Calculating mean line tangents for initial configuration"
            << endl;

        //R0.boundaryFieldRef().evaluateCoupled();
        dR0Ds_ = fvc::snGrad(R0);
        dR0Ds_ /= mag(dR0Ds_);

        refTangent_ = dR0Ds_;

        if (nBeams < 2)
        {
            if (refTangent_.boundaryFieldRef()[startPatchIndex()].size())
            {
                refTangent_.boundaryFieldRef()[startPatchIndex()] *= -1;
            }
            else
            {
                // Find left processor patch
                forAll(mesh().boundary(), patchI)
                {
                    if (isA<processorFvPatch>(mesh().boundary()[patchI]))
                    {
                        const labelList& fc =
                            mesh().boundary()[patchI].faceCells();

                        if (fc[0] == 0)
                        {
                            refTangent_.boundaryFieldRef()[patchI] *= -1;
                        }
                    }
                }
            }
        }
        else
        {
            for (label bI = 0; bI < nBeams; bI++)
            {
                const cellZone& cz = mesh().cellZones()[bI];
                label firstCell = min(cz);
                // Pout << bI << ", " << firstCell << endl;

                if (refTangent_.boundaryFieldRef()[startPatchIndex(bI)].size())
                {
                    refTangent_.boundaryFieldRef()[startPatchIndex(bI)] *= -1;
                }
                else
                {
                    // Find left processor patch
                    forAll(mesh().boundary(), patchI)
                    {
                        if (isA<processorFvPatch>(mesh().boundary()[patchI]))
                        {
                            const labelList& fc =
                                mesh().boundary()[patchI].faceCells();

                            // if (Pstream::myProcNo() == 1)
                            // {
                            //     Pout << "fc: " << fc << endl;
                            // }

                            forAll(fc, cI)
                            {
                                if (cz.whichCell(fc[cI]) != -1)
                                {
                                    if (fc[cI] == firstCell)
                                    {
                                        // if (Pstream::myProcNo() == 1)
                                        // {
                                        //     Pout << fc[cI] << " == " << firstCell << endl;
                                        //     Pout << refTangent_.boundaryField()[patchI] << endl;
                                        // }

                                        refTangent_.boundaryFieldRef()[patchI][cI]
                                            *= -1;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    // Here, the values of tangent fields set by the user for an initially curved
    // or translated (and rotated) beam are assigned to dR0Ds_
    else
    {
        dR0Ds_ = refTangent_;
        for (label pI = 0; pI < nBeams; pI++)
        {
            dR0Ds_.boundaryFieldRef()[startPatchIndex(pI)] *= -1;
        }
    }

    // In the original foam-extend version, the cell-centre rotation fields
    // were called refRM and RM. For consistency they are now renamed as
    // refLambda and Lambda in OpenFOAM-v2306.
    // For backwards compatibility, we will check for refRM and RM and report
    // an error if found
    IOobject refRMheader
    (
        "refRM",
        runTime.timeName(),
        mesh(),
        IOobject::MUST_READ
    );
    IOobject RMheader
    (
        "RM",
        runTime.timeName(),
        mesh(),
        IOobject::MUST_READ
    );
    // The check for refRM and RM files, also warn user about refLambda and
    // Lambda datatypes
    if
    (
        RMheader.typeHeaderOk<volTensorField>(true)
     || refRMheader.typeHeaderOk<volTensorField>(true)
    )
    {
         FatalErrorIn
        (
            "Constructor of Foam::coupledTotalLagNewtonRaphsonBeam "
        )   << "The variables 'refRM' or 'RM' are found in the case files "
            << "which are deprecated for OpenFOAM-v23xx version. " << nl
            << "Please rename volTensorField refRM--> volTensorField refLambda "
            << "and volTensorField RM--> volTensorField Lambda.  " << nl
            << "A better solution is to clean all the files, and re-run the "
            << "simulation."
            << abort(FatalError);
    }


    // Calculate cell-centre reference rotation matrix
    // if it is not read
    IOobject refLambdaHeader
    (
        "refLambda",
        runTime.timeName(),
        mesh(),
        IOobject::MUST_READ
    );
    if (!refLambdaHeader.typeHeaderOk<volTensorField>(true))
    {
        Info<< "Calculating cell-centre reference rotation matrix" << endl;

        // Calc cell-centre reference rotation matrix
        interpolateRotationMatrix(*this, refLambdaf_, refLambda_);

        if (debug)
        {
            Info<< refLambda_ << endl;
        }
    }


    // Calc element lengths
    if (mesh().nCells() > 0)
    {
        if (nBeams < 2)
        {
            surfaceVectorField curCf(mesh().Cf() + refWf_);
            vectorField beamPoints(this->beamPointData(curCf));
            vectorField beamTangents(this->beamPointData(refTangent_));

            HermiteSpline spline
            (
                beamPoints,
                beamTangents
            );

            // Set segment lengths
            this->L().primitiveFieldRef() = spline.segLengths();
        }
        else
        {
            surfaceVectorField curCf(mesh().Cf() + refWf_);
            for (label i = 0; i < nBeams; i++)
            {
                vectorField beamPoints(this->beamPointData(curCf, i));
                vectorField beamTangents(this->beamPointData(refTangent_, i));

                if (beamPoints.size())
                {
                    HermiteSpline spline
                    (
                        beamPoints,
                        beamTangents
                    );

                    const labelList& curBeamCells = mesh().cellZones()[i];
                    scalarField curSegLengths(spline.segLengths());
                    // vectorField newCellCentres = spline.midPoints();

                    forAll(curBeamCells, cellI)
                    {
                        // Set segment lengths
                        label curCell = curBeamCells[cellI];
                        this->L().primitiveFieldRef()[curCell] =
                            curSegLengths[cellI];

                    }
                }
            }
        }
    }

    W_.storePrevIter();
    Theta_.storePrevIter();


    // Update contact for restart
    // if (contactActive())
    // {
    //     if (runTime.timeIndex()>1)
    //     {
        //     // Info<< "Inside if" << endl;
    //         contact().update();
    //         contact().finalUpdate();
    //     }
    //     else
    //     {
        //     // Info<< "else loop" << endl;
    //         contact();
    //     }
    // }


    // Set fields as oriented
    explicitQ_.setOriented(true);
    explicitM_.setOriented(true);
    Q_.setOriented(true);
    M_.setOriented(true);

    // Create momentumContributions
    // In beam-fluid interaction cases, beams are separate regions within the
    // fluid domain (e.g. constant/beam_0/). We use runTime.caseConstant()/regionName
    // to read from the main case's constant directory, not the processor's.
    // This mirrors the pattern used for beamProperties in newBeamModel.C.
    const word& regionName = mesh().name();
    IOobject momentumContribHeader
    (
        "beamMomentumContributionProperties",
        (regionName == polyMesh::defaultRegion)
          ? fileName(runTime.caseConstant())
          : fileName(runTime.caseConstant()/regionName),
        runTime,
        IOobject::MUST_READ
    );

    if (momentumContribHeader.typeHeaderOk<dictionary>(true))
    {
        Info<< "Found beamMomentumContributionProperties file: creating object"
            << endl;

        // TO CHECK: maybe we should not register this dict
        IOdictionary beamMomentContribDict(momentumContribHeader);

        const PtrList<entry> entries
        (
            beamMomentContribDict.lookup("beamMomentumContributions")
        );

        const label nContrib = entries.size();

        momentumContribPtr_.setSize(nContrib);
        // dictionary& momentumContribDicts = dict.lookup("beamContributions");

        forAll(entries, i)
        {
            momentumContribPtr_.set
            (
                i,
                beamMomentumContribution::New
                (
                    word(entries[i].dict().lookup("beamMomentumContributionType")),
                    entries[i].dict()
                )
            );

            // Stop before the first time step if the beam does not suit the
            // contribution
            momentumContribPtr_[i].checkBeam(*this);
        }
    }
    else
    {
        Info<< "beamMomentumContributionProperties file not found: skipping"
            << endl;
    }
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void coupledTotalLagNewtonRaphsonBeam::updateTotalFields()
{
    surfaceVectorField Wf(fvc::interpolate(W_) + refWf_);

    const vectorField& WfI = Wf.internalField();

    const tensorField& LambdafI = Lambdaf_.internalField();
    const tensorField& refLambdafI = refLambdaf_.internalField();
    const faceList& faces = mesh().faces();

    const vectorField& points = mesh().points();

    vectorField& pointWI(pointW_.primitiveFieldRef());

    forAll(WfI, faceI)
    {
        const face& curFace = faces[faceI];

        const vector C0 = curFace.centre(points);

        forAll(curFace, pointI)
        {
            label curPoint = curFace[pointI];

            vector oldR = points[curPoint] - C0;

            vector newR = C0 + WfI[faceI] + //(LambdafI[faceI] & oldR);
                (LambdafI[faceI] & (refLambdafI[faceI] & oldR));

            pointWI[curPoint] = newR - points[curPoint];
        }
    }

    forAll(Wf.boundaryField(), patchI)
    {
        const vectorField& pWf =
            Wf.boundaryField()[patchI];

        const tensorField& pLambdaf =
            Lambdaf_.boundaryField()[patchI];

       const tensorField& pRefLambdaf =
            refLambdaf_.boundaryField()[patchI];

        const label start =
            mesh().boundaryMesh()[patchI].start();

        forAll(pWf, faceI)
        {
            const face& curFace = faces[start + faceI];
            vector C0 = curFace.centre(points);

            forAll(curFace, pointI)
            {
                label curPoint = curFace[pointI];

                vector oldR = points[curPoint] - C0;

                vector newR = C0 + pWf[faceI] +
                 (pLambdaf[faceI] & (pRefLambdaf[faceI] & oldR));

                pointWI[curPoint] = newR - points[curPoint];
            }
        }
    }

    beamModel::updateTotalFields();
}


void coupledTotalLagNewtonRaphsonBeam::writeFields()
{
    beamModel::writeFields();
}


tmp<vectorField> coupledTotalLagNewtonRaphsonBeam::currentBeamPoints
(
    const label bI
) const
{
    const fvMesh& mesh = this->mesh();

    label nCellZones = mesh.cellZones().size();

    if (nCellZones < 2)
    {
        label nPoints = mesh.nCells() + 1;

        tmp<vectorField> tCurrentPoints
        (
            new vectorField(nPoints, vector::zero)
        );

        const surfaceVectorField Wf(fvc::interpolate(W_));
        surfaceVectorField curCf(mesh.Cf() + refWf_ + Wf);

        // Calculate cell-centre mean line tangents
        vectorField tangents(mesh.nCells(), vector::zero);
        const labelList& upperNeiFaces = upperNeiCellFaces();
        const labelList& lowerNeiFaces = lowerNeiCellFaces();
        forAll(tangents, cellI)
        {
            label upperFace = upperNeiFaces[cellI];
            label lowerFace = lowerNeiFaces[cellI];

            if (upperFace > -1)
            {
                tangents[cellI] = curCf[upperFace];
            }
            else
            {
                label patchID = -upperFace - 1;
                tangents[cellI] = curCf.boundaryField()[patchID][0];
            }

            if (lowerFace > -1)
            {
                tangents[cellI] -= curCf[lowerFace];
            }
            else
            {
                label patchID = -lowerFace - 1;
                tangents[cellI] -= curCf.boundaryField()[patchID][0];
            }

            tangents[cellI] /= mag(tangents[cellI]) + SMALL;
        }

        const volVectorField& C = mesh.C();
        const volVectorField curC(C + refW_ + W_);

        // Correct internal faces mean line position
        const labelList& own = mesh.owner();
        const labelList& nei = mesh.neighbour();
        forAll(curCf, faceI)
        {
            curCf[faceI] =
                HermiteSpline::position
                (
                    curC[own[faceI]],
                    tangents[own[faceI]],
                    curC[nei[faceI]],
                    tangents[nei[faceI]],
                    0
                );
        }

        tCurrentPoints.ref() = this->beamPointData(curCf);

        return tCurrentPoints;
    }
    else
    {
        const cellZone& cz = mesh.cellZones()[bI];
        label nPoints = cz.size() + 1;

        tmp<vectorField> tCurrentPoints
        (
            new vectorField(nPoints, vector::zero)
        );

        const surfaceVectorField Wf(fvc::interpolate(W_));
        surfaceVectorField curCf(mesh.Cf() + refWf_ + Wf);

        // Calculate cell-centre mean line tangents
        vectorField tangents(mesh.nCells(), vector::zero);
        const labelList& upperNeiFaces = upperNeiCellFaces();
        const labelList& lowerNeiFaces = lowerNeiCellFaces();
        forAll(tangents, cellI)
        {
            label upperFace = upperNeiFaces[cellI];
            label lowerFace = lowerNeiFaces[cellI];

            if (upperFace > -1)
            {
                tangents[cellI] = curCf[upperFace];
            }
            else
            {
                label patchID = -upperFace - 1;

                const labelList& fc =
                    mesh.boundary()[patchID].faceCells();

                forAll(curCf.boundaryField()[patchID], fI)
                {
                    if (fc[fI] == cellI)
                    {
                        tangents[cellI] = curCf.boundaryField()[patchID][fI];
                    }
                }
            }

            if (lowerFace > -1)
            {
                tangents[cellI] -= curCf[lowerFace];
            }
            else
            {
                label patchID = -lowerFace - 1;

                const labelList& fc =
                    mesh.boundary()[patchID].faceCells();

                forAll(curCf.boundaryField()[patchID], fI)
                {
                    if (fc[fI] == cellI)
                    {
                        tangents[cellI] -= curCf.boundaryField()[patchID][fI];
                    }
                }
            }

            tangents[cellI] /= mag(tangents[cellI]) + SMALL;
        }

        const volVectorField& C = mesh.C();
        const volVectorField curC(C + refW_ + W_);

        // Correct internal faces mean line position
        const labelList& own = mesh.owner();
        const labelList& nei = mesh.neighbour();
        forAll(curCf, faceI)
        {
            curCf[faceI] =
                HermiteSpline::position
                (
                    curC[own[faceI]],
                    tangents[own[faceI]],
                    curC[nei[faceI]],
                    tangents[nei[faceI]],
                    0
                );
        }

        tCurrentPoints.ref() = this->beamPointData(curCf, bI);

        return tCurrentPoints;
    }
}


tmp<vectorField> coupledTotalLagNewtonRaphsonBeam::currentGlobalBeamPoints
(
    const label bI
) const
{
    if (Pstream::parRun())
    {
        List<pointField> allPoints(Pstream::nProcs());
        allPoints[Pstream::myProcNo()] = currentBeamPoints(bI);

        // Number of point per processors
        labelList procNPoints(Pstream::nProcs(), 0);
        procNPoints[Pstream::myProcNo()] =
            allPoints[Pstream::myProcNo()].size();
        Pstream::gatherList(procNPoints);
        Pstream::scatterList(procNPoints);

        forAll(allPoints, procI)
        {
            if (procI != Pstream::myProcNo())
            {
                allPoints[procI] =
                    pointField(procNPoints[procI], vector::zero);
            }
        }


        // Send data
        for (label domainI = 0; domainI < Pstream::nProcs(); ++domainI)
        {
            if (domainI != Pstream::myProcNo())
            {
                // if (allPoints[Pstream::myProcNo()].size())
                {
                    // Parallel data exchange
                    OPstream::write
                    (
                        Pstream::commsTypes::blocking,
                        domainI,
                        reinterpret_cast<const char*>
                        (
                            allPoints[Pstream::myProcNo()].begin()
                        ),
                        allPoints[Pstream::myProcNo()].byteSize()
                    );
                }
            }
        }

        // Pout << "data send" << endl;
        // sleep(5);

        // Receive data
        for (label domainI = 0; domainI < Pstream::nProcs(); ++domainI)
        {
            if (domainI != Pstream::myProcNo())
            {
                // if (procNPoints[domainI])
                {
                    // Parallel data exchange
                    IPstream::read
                    (
                        Pstream::commsTypes::blocking,
                        domainI,
                        reinterpret_cast<char*>(allPoints[domainI].begin()),
                        allPoints[domainI].byteSize()
                    );
                }
            }
        }

        // Pout << "data received" << endl;
        // sleep(5);

        for (label domainI = 0; domainI < Pstream::nProcs(); ++domainI)
        {
            if (allPoints[domainI].size())
            {
                return tmp<pointField>(new pointField(allPoints[domainI]));
                break;
            }
        }
    }

    return currentBeamPoints(bI);
}


tmp<vectorField> coupledTotalLagNewtonRaphsonBeam::currentBeamTangents
(
    const label bI
) const
{
    const fvMesh& mesh = this->mesh();
    label nCellZones = mesh.cellZones().size();

    if (nCellZones < 2)
    {
        label nPoints = mesh.nCells() + 1;

        tmp<vectorField> tCurrentTangents
        (
            new vectorField(nPoints, vector::zero)
        );

        volVectorField R0
        (
            IOobject
            (
                "R0",
                runTime().timeName(),
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
             ),
            mesh,
            dimensionedVector("R0", dimLength, vector::zero)
        );
        R0 = mesh.C() + refW_ + W_;


        //R0.boundaryField().evaluateCoupled();
        surfaceVectorField curTangents(fvc::snGrad(R0));
        curTangents /= mag(curTangents);

        if (curTangents.boundaryField()[startPatchIndex()].size())
        {
            curTangents.boundaryFieldRef()[startPatchIndex()] *= -1;
        }
        else
        {
            // Find left processor patch
            forAll(mesh.boundary(), patchI)
            {
                if (isA<processorFvPatch>(mesh.boundary()[patchI]))
                {
                    const labelList& fc =
                        mesh.boundary()[patchI].faceCells();

                    if (fc[0] == 0)
                    {
                        curTangents.boundaryFieldRef()[patchI] *= -1;
                    }
                }
            }
        }

        tCurrentTangents.ref() = this->beamPointData(curTangents);

        return tCurrentTangents;
    }
    else
    {
        const cellZone& cz = mesh.cellZones()[bI];
        label nPoints = cz.size() + 1;

        tmp<vectorField> tCurrentTangents
        (
            new vectorField(nPoints, vector::zero)
        );

        volVectorField R0
        (
            IOobject
            (
                "R0",
                runTime().timeName(),
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedVector("R0", dimLength, vector::zero)
        );
        R0 = mesh.C() + refW_ + W_;
        //R0.boundaryField().evaluateCoupled();
        surfaceVectorField curTangents(fvc::snGrad(R0));
        curTangents /= mag(curTangents);

        // curTangents.boundaryField()[startPatchIndex(bI)] *= -1;

        if (curTangents.boundaryField()[startPatchIndex(bI)].size())
        {
            curTangents.boundaryFieldRef()[startPatchIndex(bI)] *= -1;
        }
        else
        {
            // Find left (start) processor patch
            forAll(mesh.boundary(), patchI)
            {
                if (isA<processorFvPatch>(mesh.boundary()[patchI]))
                {
                    const labelList& fc =
                        mesh.boundary()[patchI].faceCells();

                    forAll(fc, fI)
                    {
                        if (cz.whichCell(fc[fI]) != -1)
                        {
                            if (fc[fI] <= min(cz))
                            {
                                curTangents.boundaryFieldRef()[patchI][fI] *= -1;
                            }
                        }
                    }
                }
            }
        }

        tCurrentTangents.ref() = this->beamPointData(curTangents, bI);

        return tCurrentTangents;
    }
}



tmp<vectorField> coupledTotalLagNewtonRaphsonBeam::currentGlobalBeamTangents
(
    const label bI
) const
{
    if (Pstream::parRun())
    {
        List<pointField> allTangents(Pstream::nProcs());
        allTangents[Pstream::myProcNo()] = currentBeamTangents(bI);

        // Number of point per processors
        labelList procNTangents(Pstream::nProcs(), 0);
        procNTangents[Pstream::myProcNo()] =
            allTangents[Pstream::myProcNo()].size();
        Pstream::gatherList(procNTangents);
        Pstream::scatterList(procNTangents);

        forAll(allTangents, procI)
        {
            if (procI != Pstream::myProcNo())
            {
                allTangents[procI] =
                    pointField(procNTangents[procI], vector::zero);
            }
        }

        // PstreamBuffers pBufs(Pstream::nonBlocking);


        // Send data
        for (label domainI = 0; domainI < Pstream::nProcs(); ++domainI)
        {
            if (domainI != Pstream::myProcNo())
            {
                if (allTangents[Pstream::myProcNo()].size())
                {
                    // Parallel data exchange
                    OPstream::write
                    (
                        Pstream::commsTypes::blocking,
                        domainI,
                        reinterpret_cast<const char*>
                        (
                            allTangents[Pstream::myProcNo()].begin()
                        ),
                        allTangents[Pstream::myProcNo()].byteSize()
                    );
                }
            }
        }

        // Pout << "data send" << endl;
        // sleep(5);

        // Receive data
        for (label domainI = 0; domainI < Pstream::nProcs(); ++domainI)
        {
            if (domainI != Pstream::myProcNo())
            {
                if (procNTangents[domainI])
                {
                    // Parallel data exchange
                    IPstream::read
                    (
                        Pstream::commsTypes::blocking,
                        domainI,
                        reinterpret_cast<char*>(allTangents[domainI].begin()),
                        allTangents[domainI].byteSize()
                    );
                }
            }
        }

        // Pout << "data received" << endl;
        // sleep(5);

        for (label domainI = 0; domainI < Pstream::nProcs(); ++domainI)
        {
            if (allTangents[domainI].size())
            {
                return tmp<pointField>(new pointField(allTangents[domainI]));
                break;
            }
        }
    }

    return currentBeamTangents(bI);
}


void coupledTotalLagNewtonRaphsonBeam::currentGlobalBeamPointsAndTangents
(
    const label bI,
    vectorField& points,
    vectorField& tangents
) const
{
    if (Pstream::parRun())
    {
        List<pointField> allPoints(Pstream::nProcs());
        allPoints[Pstream::myProcNo()] = currentBeamPoints(bI);

        Pstream::gatherList(allPoints);
        Pstream::scatterList(allPoints);

        List<pointField> allTangents(Pstream::nProcs());
        allTangents[Pstream::myProcNo()] = currentBeamTangents(bI);

        Pstream::gatherList(allTangents);
        Pstream::scatterList(allTangents);

        labelListList glPtIndices(Pstream::nProcs());
        glPtIndices[Pstream::myProcNo()] = globalPointsIndices(bI);

        Pstream::gatherList(glPtIndices);
        Pstream::scatterList(glPtIndices);

        // Set size of points and tangents
        label nPoints = 0;
        forAll(glPtIndices, procI)
        {
            label curMax = max(glPtIndices[procI]);

            if (curMax > nPoints)
            {
                nPoints = curMax;
            }
        }
        nPoints++;

        points.setSize(nPoints);
        tangents.setSize(nPoints);

        // Combine data (points and tangets should be initialised befor)
        for (label domainI = 0; domainI < Pstream::nProcs(); ++domainI)
        {
            forAll(allPoints[domainI], pI)
            {
                points[glPtIndices[domainI][pI]] = allPoints[domainI][pI];
                tangents[glPtIndices[domainI][pI]] = allTangents[domainI][pI];
            }
        }
    }
    else
    {
        points = currentBeamPoints(bI);
        tangents = currentBeamTangents(bI);
    }
}


void coupledTotalLagNewtonRaphsonBeam::currentGlobalBeamPointsAndTangents
(
    const label bI,
    const labelListList& procSendPoints,
    vectorField& points,
    vectorField& tangents
) const
{
    if (Pstream::parRun())
    {
        List<pointField> allPoints(Pstream::nProcs());
        allPoints[Pstream::myProcNo()] = currentBeamPoints(bI);

        // Pstream::gatherList(allPoints);
        // Pstream::scatterList(allPoints);

        List<pointField> allTangents(Pstream::nProcs());
        allTangents[Pstream::myProcNo()] = currentBeamTangents(bI);

        // Pstream::gatherList(allTangents);
        // Pstream::scatterList(allTangents);

        // labelListList glPtIndices(Pstream::nProcs());
        // glPtIndices[Pstream::myProcNo()] =
        //     localToGlobalBeamPointsAddressing()[bI];
        // glPtIndices[Pstream::myProcNo()] = globalPointsIndices(bI);

        // Pstream::gatherList(glPtIndices);
        // Pstream::scatterList(glPtIndices);



        // if (true){

        // Comunicate number of points which will be send to
        // neighbour processors
        labelListList procNumSendPoints(Pstream::nProcs());

        procNumSendPoints[Pstream::myProcNo()] =
            labelList(Pstream::nProcs(), 0);

        forAll(procNumSendPoints[Pstream::myProcNo()], procI)
        {
            if (procI != Pstream::myProcNo())
            {
                procNumSendPoints[Pstream::myProcNo()][procI] =
                    procSendPoints[procI].size();
            }
            // else
            // {
            //     procNumSendPoints[Pstream::myProcNo()][procI] = 0;
            // }
        }

        Pstream::gatherList(procNumSendPoints);
        Pstream::scatterList(procNumSendPoints);

        // // Number of point per processors
        // labelList procNPoints(Pstream::nProcs(), 0);
        // procNPoints[Pstream::myProcNo()] =
        //     allPoints[Pstream::myProcNo()].size();

        // Pstream::gatherList(procNPoints);
        // Pstream::scatterList(procNPoints);

        forAll(allPoints, procI)
        {
            if (procI != Pstream::myProcNo())
            {
                label curProcNumReceivePoints =
                    procNumSendPoints[procI][Pstream::myProcNo()];
                    // procReceivePoints[procI].size();

                // Pout << curProcNumReceivePoints << endl;

                allPoints[procI] =
                    pointField(curProcNumReceivePoints, vector::zero);
                allTangents[procI] =
                    pointField(curProcNumReceivePoints, vector::zero);

                // allPoints[procI] =
                //     pointField(procNPoints[procI], vector::zero);
                // allTangents[procI] =
                //     pointField(procNPoints[procI], vector::zero);
            }
        }

        // Send and receive point data
        labelListList procReceivePoints_(Pstream::nProcs());
        // procReceivePoints[Pstream::myProcNo()] =
        //     procSendPoints[Pstream::myProcNo()];
        {
            // Send data
            for (label domainI = 0; domainI < Pstream::nProcs(); ++domainI)
            {
                if (domainI != Pstream::myProcNo())
                {
                    // if (allPoints[Pstream::myProcNo()].size())
                    {
                        // Parallel data exchange
                        OPstream::write
                        (
                            Pstream::commsTypes::blocking,
                            domainI,
                            reinterpret_cast<const char*>
                            (
                                procSendPoints[domainI].begin()
                            ),
                            procSendPoints[domainI].byteSize()
                        );
                    }
                }
            }

            // Receive data
            for (label domainI = 0; domainI < Pstream::nProcs(); ++domainI)
            {
                if (domainI != Pstream::myProcNo())
                {
                    procReceivePoints_[domainI].setSize
                    (
                        procNumSendPoints[domainI][Pstream::myProcNo()]
                    );

                    // if (procNPoints[domainI])
                    {
                        // Parallel data exchange
                        IPstream::read
                        (
                            Pstream::commsTypes::blocking,
                            domainI,
                            reinterpret_cast<char*>
                            (
                                procReceivePoints_[domainI].begin()
                            ),
                            procReceivePoints_[domainI].byteSize()
                        );
                    }
                }
            }
        }

        // sleep(5);

        // PstreamBuffers pBufs(Pstream::nonBlocking);

        // Send and receive points data
        {
            for (label domainI = 0; domainI < Pstream::nProcs(); ++domainI)
            {
                if (domainI != Pstream::myProcNo())
                {
                    // if (allPoints[Pstream::myProcNo()].size())
                    {
                        vectorField pointsToSend
                        (
                            allPoints[Pstream::myProcNo()],
                            procSendPoints[domainI]
                        );

                        // Parallel data exchange
                        OPstream::write
                        (
                            Pstream::commsTypes::blocking,
                            domainI,
                            reinterpret_cast<const char*>
                            (
                                pointsToSend.begin()
                                // allPoints[Pstream::myProcNo()].begin()
                            ),
                            pointsToSend.byteSize()
                            // allPoints[Pstream::myProcNo()].byteSize()
                        );
                    }
                }
            }

            // Receive data
            for (label domainI = 0; domainI < Pstream::nProcs(); ++domainI)
            {
                if (domainI != Pstream::myProcNo())
                {
                    // if (procNPoints[domainI])
                    {
                        // Parallel data exchange
                        IPstream::read
                        (
                            Pstream::commsTypes::blocking,
                            domainI,
                            reinterpret_cast<char*>(allPoints[domainI].begin()),
                            allPoints[domainI].byteSize()
                        );
                    }
                }
            }
        }

        // Send and receive tangents data
        {
            for (label domainI = 0; domainI < Pstream::nProcs(); ++domainI)
            {
                if (domainI != Pstream::myProcNo())
                {
                    // if (allPoints[Pstream::myProcNo()].size())
                    {
                        vectorField tangentsToSend
                        (
                            allTangents[Pstream::myProcNo()],
                            procSendPoints[domainI]
                        );

                        OPstream::write
                        (
                            Pstream::commsTypes::blocking,
                            domainI,
                            reinterpret_cast<const char*>
                            (
                                tangentsToSend.begin()
                                // allTangents[Pstream::myProcNo()].begin()
                            ),
                            tangentsToSend.byteSize()
                            // allTangents[Pstream::myProcNo()].byteSize()
                        );
                    }
                }
            }

            // Receive data
            for (label domainI = 0; domainI < Pstream::nProcs(); ++domainI)
            {
                if (domainI != Pstream::myProcNo())
                {
                    // if (procNPoints[domainI])
                    {
                        // Parallel data exchange
                        IPstream::read
                        (
                            Pstream::commsTypes::blocking,
                            domainI,
                            reinterpret_cast<char*>
                            (
                                allTangents[domainI].begin()
                            ),
                            allTangents[domainI].byteSize()
                        );
                    }
                }
            }
        }
        // }

        // Combine data
        for (label procI = 0; procI < Pstream::nProcs(); ++procI)
        {
            if (procI != Pstream::myProcNo())
            {
                forAll(allPoints[procI], pI)
                {
                    label locPointIndex = procReceivePoints_[procI][pI];
                    label glPointIndex =
                        localToGlobalBeamPointsAddressing()[bI][procI][locPointIndex];

                    points[glPointIndex] = allPoints[procI][pI];
                    tangents[glPointIndex] = allTangents[procI][pI];

                    // points[glPtIndices[procI][locPointIndex]] = allPoints[procI][pI];
                    // tangents[glPtIndices[procI][locPointIndex]] = allTangents[procI][pI];
                }
            }
            else
            {
                forAll(allPoints[procI], pI)
                {
                    label glPointIndex =
                        localToGlobalBeamPointsAddressing()[bI][procI][pI];

                    points[glPointIndex] = allPoints[procI][pI];
                    tangents[glPointIndex] = allTangents[procI][pI];

                    // points[glPtIndices[procI][pI]] = allPoints[procI][pI];
                    // tangents[glPtIndices[procI][pI]] = allTangents[procI][pI];
                }
            }
        }

        if (false){
        // label nPoints = sum(procNPoints);
        // pointField globalPoints(nPoints, vector::zero);
        // pointField globalTangents(nPoints, vector::zero);
        pointField globalPoints = points;
        pointField globalTangents = tangents;

        // Combine data
        label pointI = 0;
        for (label domainI = 0; domainI < Pstream::nProcs(); ++domainI)
        {
            for (label pI=0; pI<(allPoints[domainI].size()-1); pI++)
            // forAll(allPoints[domainI], pI)
            {
                globalPoints[pointI] = allPoints[domainI][pI];
                globalTangents[pointI] = allTangents[domainI][pI];
                pointI++;
            }

            // Last point
            if (domainI == Pstream::nProcs()-1)
            {
                label pI = allPoints[domainI].size()-1;
                globalPoints[pointI] = allPoints[domainI][pI];
                globalTangents[pointI] = allTangents[domainI][pI];
                pointI++;
            }

            // if (allPoints[domainI].size())
            // {
            //     points = allPoints[domainI];
            //     tangents = allTangents[domainI];
            //     // return tmp<pointField>(new pointField(allPoints[domainI]));
            //     break;
            // }
        }

        globalPoints.resize(pointI);
        globalTangents.resize(pointI);

        points = globalPoints;
        tangents = globalTangents;
        }
    }
    else
    {
        points = currentBeamPoints(bI);
        tangents = currentBeamTangents(bI);
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace beamModels

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
