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

#include "vertexCentredLinGeomSolidMMS.H"
#include "addToRunTimeSelectionTable.H"
#include "sparseMatrix.H"
#include "symmTensor4thOrder.H"
#include "vfvcCellPoint.H"
#include "vfvmCellPoint.H"
#include "fvcDiv.H"
#include "fixedValuePointPatchFields.H"
#include "solidTractionPointPatchVectorField.H"
#include "mathematicalConstants.H"
#include "boundBox.H"
#include "sparseMatrixTools.H"
#include "symmetryPointPatchFields.H"
#include "fixedDisplacementZeroShearPointPatchVectorField.H"
#ifdef USE_PETSC
    #include <petscksp.h>
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace solidModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(vertexCentredLinGeomSolidMMS, 0);
addToRunTimeSelectionTable(solidModel, vertexCentredLinGeomSolidMMS, dictionary); //adds derived contructor to hash table held by base
//addToRunTimeSelectionTable(baseType, thisType, argNames)


// * * * * * * * * * * *  Private Member Functions * * * * * * * * * * * * * //

void vertexCentredLinGeomSolidMMS::updateSource
(
    vectorField& source,
    const labelList& dualCellToPoint
)
{
    if (debug)
    {
        Info<< "void vertexCentredLinGeomSolidMMS::updateSource(...): start"
            << endl;
    }

    // Reset to zero
    source = vector::zero;

    // The source vector is -F, where:
    // F = div(sigma) + rho*g - rho*d2dt2(D)

    // Point volume field
    const scalarField& pointVolI = pointVol_.internalField();

    // Point density field
    //const scalarField& pointRhoI = pointRho_.internalField();
    
    //Coefficients for expression
    const scalar ax = 2;
    const scalar ay = 4;
    const scalar az = 6;
    
    //pi
    const scalar pi = constant::mathematical::pi; 
    
	// Young's modulus
    const scalar E = 200e9;

    // Poisson's ratio
    const scalar nu = 0.3;

    // Set the shear modulus
    const scalar mu = E/(2.0*(1.0 + nu));

    // Set lambda
    const scalar lambda = (E*nu)/((1.0 + nu)*(1.0 - 2.0*nu));

    // Calculate the tractions on the dual faces
    surfaceVectorField dualTraction
    (
        (dualMesh().Sf()/dualMesh().magSf()) & dualSigmaf_
    );

    // Enforce extract tractions on traction boundaries

    enforceTractionBoundaries
    (
        pointD(), dualTraction, mesh(), dualMeshMap().pointToDualFaces()
    );

    // Set coupled boundary (e.g. processor) traction fields to zero: this
    // ensures their global contribution is zero
    forAll(dualTraction.boundaryField(), patchI)
    {
        if (dualTraction.boundaryField()[patchI].coupled())
        {
#ifdef OPENFOAMESIORFOUNDATION
            dualTraction.boundaryFieldRef()[patchI] = vector::zero;
#else
            dualTraction.boundaryField()[patchI] = vector::zero;
#endif
        }
    }

    // Calculate divergence of stress for the dual cells
    // const vectorField dualDivSigma = fvc::div(dualMesh().Sf() & dualSigmaf_);
    
    const vectorField dualDivSigma = fvc::div(dualTraction*dualMesh().magSf());
    
	// Map dual cell field to primary mesh point field
    vectorField pointDivSigma(mesh().nPoints(), vector::zero);
    
    // Body force field
    vectorField bodyForcesI(mesh().nPoints(), vector::zero);
    
//    // List of primary cells
//    const cellList& cells = mesh().cells();
 
    // List of primary points
    const pointField& points = mesh().points();

//    // List of primary faces
//    const faceList& faces = mesh().faces();

//    // List of dual cell volumes
//    const scalarField& V = mesh().V();
    
	forAll(bodyForcesI, vertexI)
	{

	    // Take a reference to the current primary volume
//	    const scalar curV = V[vertexI];
//	    
//	    //Take a reference to the current primary cell
//	    const cell& curCell = cells[vertexI];
//	    
//	    //Create a bounding box around the current dual cell
//	    const boundBox bb(curCell.points(faces, points), false);
//	    
//	    //Get the most negative x and y coordinates of the dual cell
//	    const scalar x1 = bb.min().x();
//	    const scalar y1 = bb.min().y();
//	    const scalar z1 = bb.min().z();
//	   
//	    //Get the most negative x and y coordinates of the dual cell
//	    const scalar x2 = bb.max().x();
//	    const scalar y2 = bb.max().y();
//	    const scalar z2 = bb.max().z();
	    
	    const scalar x = points[vertexI].x();
	    const scalar y = points[vertexI].y();
	    const scalar z = points[vertexI].z();
	                                                                                                                                             
		//Set vector in dualCellI for x-equation     	
    	bodyForcesI[vertexI].x() =
    	lambda*(8*ay*pi*pi*Foam::cos(4*pi*x)*Foam::cos(2*pi*y)*Foam::sin(pi*z) 
    	+ 4*az*pi*pi*Foam::cos(4*pi*x)*Foam::cos(pi*z)*Foam::sin(2*pi*y) 
    	- 16*ax*pi*pi*Foam::sin(4*pi*x)*Foam::sin(2*pi*y)*Foam::sin(pi*z)) 
    	+ mu*(8*ax*pi*pi*Foam::cos(4*pi*x)*Foam::cos(2*pi*y)*Foam::sin(pi*z) 
    	+ 4*az*pi*pi*Foam::cos(4*pi*x)*Foam::cos(pi*z)*Foam::sin(2*pi*y) 
    	- 5*ax*pi*pi*Foam::sin(4*pi*x)*Foam::sin(2*pi*y)*Foam::sin(pi*z)) 
    	- 32*ax*mu*pi*pi*Foam::sin(4*pi*x)*Foam::sin(2*pi*y)*Foam::sin(pi*z);

		//Set vector in dualCellI for y-equation  
    	bodyForcesI[vertexI].y() =
    	lambda*(8*ax*pi*pi*Foam::cos(4*pi*x)*Foam::cos(2*pi*y)*Foam::sin(pi*z) 
    	+ 2*az*pi*pi*Foam::cos(2*pi*y)*Foam::cos(pi*z)*Foam::sin(4*pi*x) 
    	- 4*ay*pi*pi*Foam::sin(4*pi*x)*Foam::sin(2*pi*y)*Foam::sin(pi*z)) 
    	+ mu*(8*ax*pi*pi*Foam::cos(4*pi*x)*Foam::cos(2*pi*y)*Foam::sin(pi*z) 
    	+ 2*az*pi*pi*Foam::cos(2*pi*y)*Foam::cos(pi*z)*Foam::sin(4*pi*x) 
    	- 17*ay*pi*pi*Foam::sin(4*pi*x)*Foam::sin(2*pi*y)*Foam::sin(pi*z)) 
    	- 8*ax*mu*pi*pi*Foam::sin(4*pi*x)*Foam::sin(2*pi*y)*Foam::sin(pi*z);
 
    	//Set vector in dualCellI for -equation  
    	bodyForcesI[vertexI].z() =
    	lambda*(4*ax*pi*pi*Foam::cos(4*pi*x)*Foam::cos(pi*z)*Foam::sin(2*pi*y) 
    	+ 2*ay*pi*pi*Foam::cos(2*pi*y)*Foam::cos(pi*z)*Foam::sin(4*pi*x) 
    	- az*pi*pi*Foam::sin(4*pi*x)*Foam::sin(2*pi*y)*Foam::sin(pi*z)) 
    	+ mu*(4*ax*pi*pi*Foam::cos(4*pi*x)*Foam::cos(pi*z)*Foam::sin(2*pi*y) 
    	+ 2*ay*pi*pi*Foam::cos(2*pi*y)*Foam::cos(pi*z)*Foam::sin(4*pi*x) 
    	- 20*az*pi*pi*Foam::sin(4*pi*x)*Foam::sin(2*pi*y)*Foam::sin(pi*z)) 
    	- 2*az*mu*pi*pi*Foam::sin(4*pi*x)*Foam::sin(2*pi*y)*Foam::sin(pi*z);
    }
    
    // Add body force field to source
    //source -= bodyForcesI;
    
    //Info << "source: " << source;
    
    forAll(dualDivSigma, dualCellI)
    {
        const label pointID = dualCellToPoint[dualCellI];
        pointDivSigma[pointID] = dualDivSigma[dualCellI];
    }

    // Add surface forces to source
    source -= pointDivSigma*pointVolI - bodyForcesI*pointVolI;

//    // Add gravity body forces
//    source -= pointRhoI*g().value()*pointVolI;

//    // Add transient term
//    source += vfvc::d2dt2
//    (
//#ifdef OPENFOAMESIORFOUNDATION
//        mesh().d2dt2Scheme("d2dt2(pointD)"),
//#else
//        mesh().schemesDict().d2dt2Scheme("d2dt2(pointD)"),
//#endif
//        pointD(),
//        pointU_,
//        pointA_,
//        pointRho_,
//        pointVol_,
//        int(bool(debug))
//    );

    if (debug)
    {
        Info<< "void vertexCentredLinGeomSolidMMS::updateSource(...): end"
            << endl;
    }
}


void vertexCentredLinGeomSolidMMS::setFixedDofs
(
    const pointVectorField& pointD,
    boolList& fixedDofs,
    pointField& fixedDofValues,
    symmTensorField& fixedDofDirections
) const
{
    // Flag all fixed DOFs

    forAll(pointD.boundaryField(), patchI)
    {        
        if
        (
            // isA<uniformFixedValuePointPatchVectorField>
            isA<fixedValuePointPatchVectorField>
            (
                pointD.boundaryField()[patchI]
            )
        )
        {
            // const uniformFixedValuePointPatchVectorField& dispPatch =
            //     refCast<const uniformFixedValuePointPatchVectorField>
            // const fixedValuePointPatchVectorField& dispPatch =
            //     refCast<const fixedValuePointPatchVectorField>
            //     (
            //         pointD.boundaryField()[patchI]
            //     );

            // const vector& disp = dispPatch.uniformValue();

            const labelList& meshPoints =
                pointD.mesh().mesh().boundaryMesh()[patchI].meshPoints();

            forAll(meshPoints, pI)
            {
                const label pointID = meshPoints[pI];
                const vector& disp = pointD[pointID];

                // Check if this point has already been fixed
                if (fixedDofs[pointID])
                {
                    // Check if the existing prescribed displacement is
                    // consistent with the new one
                    if
                    (
                        mag
                        (
                            fixedDofDirections[pointID]
                          & (fixedDofValues[pointID] - disp)
                        ) > SMALL
                    )
                    {
                        FatalErrorIn
                        (
                            "void vertexCentredLinGeomSolidMMS::setFixedDofs(...)"
                        )   << "Inconsistent displacements prescribed at point "
                            << "= " << pointD.mesh().mesh().points()[pointID]
                            << abort(FatalError);
                    }

                    // Set all directions as fixed, just in case it was
                    // previously marked as a symmetry point
                    fixedDofDirections[pointID] = symmTensor(I);
                }
                else
                {
                    fixedDofs[pointID] = true;
                    fixedDofValues[pointID] = disp;
                    fixedDofDirections[pointID] = symmTensor(I);  
                          
                }     
            }
        }
        else if
        (
            isA<symmetryPointPatchVectorField>
            (
                pointD.boundaryField()[patchI]
            )
         || isA<fixedDisplacementZeroShearPointPatchVectorField>
            (
                pointD.boundaryField()[patchI]
            )
        )
        {
            const labelList& meshPoints =
                pointD.mesh().boundary()[patchI].meshPoints();
            const vectorField& pointNormals =
                pointD.mesh().boundary()[patchI].pointNormals();

            scalarField normalDisp(meshPoints.size(), 0.0);
            if
            (
                isA<fixedDisplacementZeroShearPointPatchVectorField>
                (
                    pointD.boundaryField()[patchI]
                )
            )
            {          
                normalDisp =
                (
                    pointNormals
                  & pointD.boundaryField()[patchI].patchInternalField()
                );
                
                if (debug)
                {
                    Info<< "normalDisp = " << normalDisp << endl;
                }
            }

            forAll(meshPoints, pI)
            {
                const label pointID = meshPoints[pI];

                // Check if this point has already been fixed
                if (fixedDofs[pointID])
                {
                    // Check if the existing prescribed displacement is
                    // consistent with the current condition
                    if
                    (
                        mag
                        (
                            (pointNormals[pI] & fixedDofValues[pointID])
                          - normalDisp[pI]
                        ) > SMALL
                    )
                    {
                        FatalErrorIn
                        (
                            "void vertexCentredLinGeomSolidMMS::setFixedDofs(...)"
                        )   << "Inconsistent displacements prescribed at point "
                            << "= " << pointD.mesh().mesh().points()[pointID]
                            << abort(FatalError);
                    }

                    // If the point is not fully fixed then make sure the normal
                    // direction is fixed
                    if (mag(fixedDofDirections[pointID] - symmTensor(I)) > 0)
                    {
                        // If the directions are orthogonal we can add them
                        const symmTensor curDir = sqr(pointNormals[pI]);
                        if (mag(fixedDofDirections[pointID] & curDir) > 0)
                        {
                            FatalError
                                << "Point " << pointID << " is fixed in two "
                                << "directions: this is only implemented for "
                                << "Cartesian axis directions" << abort(FatalError);
                        }

                        fixedDofDirections[pointID] += curDir;
                    }
                }
                else
                {
                    fixedDofs[pointID] = true;
                    fixedDofValues[pointID] = normalDisp[pI]*pointNormals[pI];
                    fixedDofDirections[pointID] = sqr(pointNormals[pI]);
                }
            }
        }
    }
}


void vertexCentredLinGeomSolidMMS::enforceTractionBoundaries
(
    const pointVectorField& pointD,
    surfaceVectorField& dualTraction,
    const fvMesh& mesh,
    const labelListList& pointToDualFaces
) const
{
    const pointMesh& pMesh = pointD.mesh();
    const fvMesh& dualMesh = dualTraction.mesh();

    forAll(pointD.boundaryField(), patchI)
    {
        if
        (
            isA<solidTractionPointPatchVectorField>
            (
                pointD.boundaryField()[patchI]
            )
        )
        {
            const solidTractionPointPatchVectorField& tracPatch =
                refCast<const solidTractionPointPatchVectorField>
                (
                    pointD.boundaryField()[patchI]
                );

            const labelList& meshPoints =
                mesh.boundaryMesh()[patchI].meshPoints();

            // Primary mesh point normals
            const vectorField& n =
                pMesh.boundary()[patchI].pointNormals();

            // Primary mesh point tractions
            const vectorField totalTraction
            (
                tracPatch.traction() - n*tracPatch.pressure()
            );

            // Create dual mesh faces traction field
            vectorField dualFaceTraction
            (
                dualMesh.boundaryMesh()[patchI].size(), vector::zero
            );

            // Multiple points map to each dual face so we will count them
            // and then divide the dualFaceTraction by this field so that it is
            // the average of all the points that map to it
            scalarField nPointsPerDualFace(dualFaceTraction.size(), 0.0);

            // Map from primary mesh point field to second mesh face field using
            // the pointToDualFaces map
            forAll(totalTraction, pI)
            {
                const label pointID = meshPoints[pI];
                const labelList& curDualFaces = pointToDualFaces[pointID];

                forAll(curDualFaces, dfI)
                {
                    const label dualFaceID = curDualFaces[dfI];

                    if (!dualMesh.isInternalFace(dualFaceID))
                    {
                        // Check which patch this dual face belongs to
                        const label dualPatchID =
                            dualMesh.boundaryMesh().whichPatch(dualFaceID);

                        if (dualPatchID == patchI)
                        {
                            // Find local face index
                            const label localDualFaceID =
                                dualFaceID
                              - dualMesh.boundaryMesh()[dualPatchID].start();

                            // Set dual face traction
                            dualFaceTraction[localDualFaceID] +=
                                totalTraction[pI];

                            // Update the count for this face
                            nPointsPerDualFace[localDualFaceID]++;
                        }
                    }
                }
            }

            if (gMin(nPointsPerDualFace) < 1)
            {
                FatalErrorIn
                (
                    "void vertexCentredLinGeomSolidMMS::"
                    "enforceTractionBoundaries(...)"
                )   << "Problem setting tractions: gMin(nPointsPerDualFace) < 1"
                    << nl << "nPointsPerDualFace = " << nPointsPerDualFace
                    << abort(FatalError);
            }

            // Take the average
            dualFaceTraction /= nPointsPerDualFace;
            // Overwrite the dual patch face traction
#ifdef OPENFOAMESIORFOUNDATION
            dualTraction.boundaryFieldRef()[patchI] = dualFaceTraction;
#else
            dualTraction.boundaryField()[patchI] = dualFaceTraction;
#endif
        }
        else if
        (
            isA<symmetryPointPatchVectorField>(pointD.boundaryField()[patchI])
         || isA<fixedDisplacementZeroShearPointPatchVectorField>
            (
                pointD.boundaryField()[patchI]
            )
        )
        {
            // Set the dual patch face shear traction to zero
            const vectorField n(dualMesh.boundary()[patchI].nf()); 
            
#ifdef OPENFOAMESIORFOUNDATION            
            dualTraction.boundaryFieldRef()[patchI] =
                (sqr(n) & dualTraction.boundaryField()[patchI]);
#else
            dualTraction.boundaryField()[patchI] =
                (sqr(n) & dualTraction.boundaryField()[patchI]);
#endif
        }
    }
}

bool vertexCentredLinGeomSolidMMS::vertexCentredLinGeomSolidMMS::converged
(
    const label iCorr,
    scalar& initResidual,
    const scalar res,
    const label nInterations,
    const pointVectorField& pointD,
    const vectorField& pointDcorr
) const
{
    // Calculate the residual as the root mean square of the correction
    const scalar residualAbs = gSum(magSqr(pointDcorr));

    // Store initial residual
    if (iCorr == 0)
    {
        initResidual = residualAbs;

        // If the initial residual is small then convergence has been achieved
        if (initResidual < SMALL)
        {
            Info<< "    Initial residual is less than 1e-15"
                << "    Converged" << endl;
            return true;
        }
        Info<< "    Initial residual = " << initResidual << endl;
    }

    // Define a normalised residual wrt the initial residual
    const scalar residualNorm = residualAbs/initResidual;

    // Calculate the maximum displacement
#ifdef OPENFOAMESIORFOUNDATION
    const scalar maxMagD = gMax(mag(pointD.primitiveField()));
#else
    const scalar maxMagD = gMax(mag(pointD.internalField()));
#endif

    // Print information
    Info<< "    Iter = " << iCorr
        << ", relRef = " << residualNorm
        << ", res = " << res
        << ", resAbs = " << residualAbs
        << ", nIters = " << nInterations
        << ", maxD = " << maxMagD << endl;

    // Check for convergence
    if (residualNorm < solutionTol())
    {
        Info<< "    Converged" << endl;
        return true;
    }
    else if (iCorr >= nCorr() - 1)
    {
        if (nCorr() > 1)
        {
            Warning
                << "Max iterations reached within the momentum Newton-Raphson "
                "loop" << endl;
        }

        return true;
    }

    // Convergence has not been reached
    return false;
}


scalar vertexCentredLinGeomSolidMMS::calculateLineSearchSlope
(
    const scalar eta,
    const vectorField& pointDcorr,
    pointVectorField& pointD,
    surfaceTensorField& dualGradDf,
    surfaceSymmTensorField& dualSigmaf,
    const scalar zeta
)
{
    // Store pointD as we will reset it after changing it
    pointD.storePrevIter();

    // Update pointD
#ifdef OPENFOAMESIORFOUNDATION
    pointD.primitiveFieldRef() += eta*pointDcorr;
#else
    pointD.internalField() += eta*pointDcorr;
#endif
    pointD.correctBoundaryConditions();

    // Calculate gradD at dual faces
    dualGradDf = vfvc::fGrad
    (
        pointD,
        mesh(),
        dualMesh(),
        dualMeshMap().dualFaceToCell(),
        dualMeshMap().dualCellToPoint(),
        zeta
    );

    // Calculate stress at dual faces
    dualMechanicalPtr_().correct(dualSigmaf);

    // Update the source vector
    vectorField source(mesh().nPoints(), vector::zero);
    pointD.correctBoundaryConditions();
    updateSource(source, dualMeshMap().dualCellToPoint());

    // Reset pointD
    pointD = pointD.prevIter();

    // Return the slope
    return gSum(pointDcorr & source);
}


scalar vertexCentredLinGeomSolidMMS::calculateLineSearchFactor
(
    const scalar rTol, // Slope reduction tolerance
    const int maxIter, // Maximum number of line search iterations
    const vectorField& pointDcorr, // Point displacement correction
    const vectorField& source, // Linear system source
    const scalar zeta // Discretisation parameter
)
{
    // Calculate initial slope
    const scalar s0 = gSum(pointDcorr & source);

    // Set initial line search parameter
    scalar eta = 1.0;
    int lineSearchIter = 0;

    // Perform backtracking to find suitable eta
    do
    {
        lineSearchIter++;

        // Calculate slope at eta
        const scalar s1 = calculateLineSearchSlope
        (
            eta, pointDcorr, pointD(), dualGradDf_, dualSigmaf_, zeta
        );

        // Calculate ratio of s1 to s0
        const scalar r = s1/s0;

        if (mag(r) < rTol)
        {
            break;
        }
        else
        {
            // Interpolate/extrapolate to find new eta
            // Limit it to be less than 10
            //eta = min(-1/(r - 1), 10);

            if (r < 0)
            {
                // Simple back tracking
                eta *= 0.5;
            }
            else
            {
                // Extrapolate
                eta = min(-1/(r - 1), 10);
            }
        }

        if (lineSearchIter == maxIter)
        {
            Warning
                << "Max line search iterations reached!" << endl;
        }
    }
    while (lineSearchIter < maxIter);

    // Update pointD and re-calculate source, then calculate s
    if (mag(eta - 1) > SMALL)
    {
        Info<< "        line search parameter = " << eta
            << ", iter = " << lineSearchIter << endl;
    }

    return eta;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

vertexCentredLinGeomSolidMMS::vertexCentredLinGeomSolidMMS
(
    Time& runTime,
    const word& region
)
:
    solidModel(typeName, runTime, region),
    dualMechanicalPtr_
    (
        new dualMechanicalModel
        (
            dualMesh(),
            nonLinGeom(),
            incremental(),
            mechanical(),
            dualMeshMap().dualFaceToCell()
        )
    ),
    fullNewton_(solidModelDict().lookup("fullNewton")),
    steadyState_(false),
    twoD_(sparseMatrixTools::checkTwoD(mesh())),
    fixedDofs_(mesh().nPoints(), false),
    fixedDofValues_(fixedDofs_.size(), vector::zero),
    fixedDofDirections_(fixedDofs_.size(), symmTensor::zero),
    fixedDofScale_
    (
        solidModelDict().lookupOrDefault<scalar>
        (
            "fixedDofScale",
            (
                average(mechanical().impK())
               *Foam::sqrt(gAverage(mesh().magSf()))
            ).value()
        )
    ),
    pointU_
    (
        IOobject
        (
            "pointU",
            runTime.timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        pMesh(),
        dimensionedVector("0", dimVelocity, vector::zero)
    ),
    pointA_
    (
        IOobject
        (
            "pointA",
            runTime.timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        pMesh(),
        dimensionedVector("0", dimVelocity/dimTime, vector::zero)
    ),
    pointRho_
    (
        IOobject
        (
            "point(rho)",
            runTime.timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        pMesh(),
        dimensionedScalar("0", dimDensity, 0.0)
    ),
    pointVol_
    (
        IOobject
        (
            "pointVolumes",
            runTime.timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        pMesh(),
        dimensionedScalar("0", dimVolume, 0.0)
    ),
    dualGradDf_
    (
        IOobject
        (
            "grad(D)f",
            runTime.timeName(),
            dualMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        dualMesh(),
        dimensionedTensor("zero", dimless, tensor::zero),
        "calculated"
    ),
    dualSigmaf_
    (
        IOobject
        (
            "sigmaf",
            runTime.timeName(),
            dualMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        dualMesh(),
        dimensionedSymmTensor("zero", dimPressure, symmTensor::zero),
        "calculated"
    ),
    globalPointIndices_(mesh())
{
    // Create dual mesh and set write option
    dualMesh().objectRegistry::writeOpt() = IOobject::NO_WRITE;

    // pointD field must be defined
    pointDisRequired();

    // Set fixed degree of freedom list
    setFixedDofs(pointD(), fixedDofs_, fixedDofValues_, fixedDofDirections_);

    // Set point density field
    mechanical().volToPoint().interpolate(rho(), pointRho_);

    // Set the pointVol field
    // Map dualMesh cell volumes to the primary mesh points
#ifdef OPENFOAMESIORFOUNDATION
    scalarField& pointVolI = pointVol_.primitiveFieldRef();
#else
    scalarField& pointVolI = pointVol_.internalField();
#endif
    const scalarField& dualCellVol = dualMesh().V();
    const labelList& dualCellToPoint = dualMeshMap().dualCellToPoint();
    forAll(dualCellToPoint, dualCellI)
    {
        // Find point which maps to this dual cell
        const label pointID = dualCellToPoint[dualCellI];

        // Map the cell volume
        pointVolI[pointID] = dualCellVol[dualCellI];
    }

    // Store old time fields
    pointD().oldTime().storeOldTime();
    pointU_.oldTime().storeOldTime();
    pointA_.storeOldTime();

    // Write fixed degree of freedom equation scale
    Info<< "fixedDofScale: " << fixedDofScale_ << endl;

    // Disable the writing of the unused fields
    D().writeOpt() = IOobject::NO_WRITE;
    D().oldTime().oldTime().writeOpt() = IOobject::NO_WRITE;
    DD().writeOpt() = IOobject::NO_WRITE;
    DD().oldTime().oldTime().writeOpt() = IOobject::NO_WRITE;
    U().writeOpt() = IOobject::NO_WRITE;
    pointDD().writeOpt() = IOobject::NO_WRITE;
}


// * * * * * * * * * * * * * * * *  Destructors  * * * * * * * * * * * * * * //

vertexCentredLinGeomSolidMMS::~vertexCentredLinGeomSolidMMS()
{
#ifdef USE_PETSC
    if (Switch(solidModelDict().lookup("usePETSc")))
    {
        PetscFinalize();
    }
#endif
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool vertexCentredLinGeomSolidMMS::evolve()
{
    Info<< "Evolving solid solver" << endl;
    
    //Update boundary conditions
    pointD().correctBoundaryConditions();

    // Initialise matrix
    sparseMatrix matrix(sum(globalPointIndices_.stencilSize())); 
    //sparseMatrix is created given the size sum(globalPointIndices_.stencilSize())

    // Store material tangent field for dual mesh faces
    Field<RectangularMatrix<scalar>> materialTangent
    (
        dualMechanicalPtr_().materialTangentFaceField()
    );
    
    // Lookup compact edge gradient factor
    const scalar zeta(solidModelDict().lookupOrDefault<scalar>("zeta", 0.2));
    if (debug)
    {
        Info<< "zeta: " << zeta << endl;
    }

    // Global point index lists
    //List showing which points are on a specific processor?
    const boolList& ownedByThisProc = globalPointIndices_.ownedByThisProc();
    //Mapping local points to global?
    const labelList& localToGlobalPointMap =
        globalPointIndices_.localToGlobalPointMap();

    if (!fullNewton_) //if not using full Newton-Raphson approach
    {
        // Assemble matrix once per time-step
        Info<< "    Assembling the matrix" << endl;

        // Add div(sigma) coefficients
        vfvm::divSigma 
        (
            matrix,
            mesh(),
            dualMesh(),
            dualMeshMap().dualFaceToCell(),
            dualMeshMap().dualCellToPoint(),
            materialTangent,
            fixedDofs_,
            fixedDofDirections_,
            fixedDofScale_,
            zeta,
            debug
        );

        // Add d2dt2 coefficients
        vfvm::d2dt2 
        (
#ifdef OPENFOAMESIORFOUNDATION
            mesh().d2dt2Scheme("d2dt2(pointD)"),
#else
            mesh().schemesDict().d2dt2Scheme("d2dt2(pointD)"),
#endif
            runTime().deltaTValue(),
            pointD().name(),
            matrix,
            pointRho_.internalField(),
            pointVol_.internalField(),
            int(bool(debug))
        );
    }

    // Solution field: point displacement correction
    vectorField pointDcorr(pointD().internalField().size(), vector::zero);

    // Newton-Raphson loop over momentum equation
    int iCorr = 0;
    scalar initResidual = 0.0;
#ifdef OPENFOAMESIORFOUNDATION
    SolverPerformance<vector> solverPerf;
#else
    BlockSolverPerformance<vector> solverPerf; //Class containing performance statistics
#endif
    do
    {
        // Calculate gradD at dual faces
        dualGradDf_ = vfvc::fGrad
        (
            pointD(),
            mesh(),
            dualMesh(),
            dualMeshMap().dualFaceToCell(),
            dualMeshMap().dualCellToPoint(),
            zeta,
            debug
        );

        // Calculate stress at dual faces
        dualMechanicalPtr_().correct(dualSigmaf_);

        // Update the source vector
        vectorField source(mesh().nPoints(), vector::zero);
        pointD().correctBoundaryConditions();
        updateSource(source, dualMeshMap().dualCellToPoint()); 

        if (fullNewton_)
        {
            // Info<< "    Assembling the matrix" << endl;

            // Assemble the matrix once per outer iteration
            matrix.clear();

            // Update material tangent
            materialTangent = dualMechanicalPtr_().materialTangentFaceField();

            // Add div(sigma) coefficients
            vfvm::divSigma
            (
                matrix,
                mesh(),
                dualMesh(),
                dualMeshMap().dualFaceToCell(),
                dualMeshMap().dualCellToPoint(),
                materialTangent,
                fixedDofs_,
                fixedDofDirections_,
                fixedDofScale_,
                zeta
            );

            // Add d2dt2 coefficients
            vfvm::d2dt2
            (
#ifdef OPENFOAMESIORFOUNDATION
                mesh().d2dt2Scheme("d2dt2(pointD)"),
#else
                mesh().schemesDict().d2dt2Scheme("d2dt2(pointD)"),
#endif
                runTime().deltaTValue(),
                pointD().name(),
                matrix,
                pointRho_.internalField(),
                pointVol_.internalField(),
                int(bool(debug))
            );
        }
        
//        Info << endl << "Before enforcing DOFs: " << endl << endl;
//        matrix.print();
//        Info << endl << "Print out the source: " << endl << endl;
//        
//        for (int i = 0; i < source.size(); i++)
//        {
//            Info << "(" << i << ", 0) : " << source[i] << endl;
//            
//        } 
//        Info << endl;

        // Enforce fixed DOF on the linear system
        sparseMatrixTools::enforceFixedDof //Loop through the matrix and overwrite coefficients for fixed DOFs
        (
            matrix,
            source,
            fixedDofs_,
            fixedDofDirections_,
            fixedDofValues_,
            fixedDofScale_
        );

//        Info << endl << "After enforcing DOFs " << endl << endl;
//        matrix.print();
//        Info << endl << "Print out the source: " << endl << endl;
//        
//        for (int i = 0; i < source.size(); i++)
//        {
//            Info << "(" << i << ", 0) : " << source[i] << endl;
//            
//        } 

        // Solve linear system for displacement correction
        if (debug)
        {
            Info<< "bool vertexCentredLinGeomSolidMMS::evolve(): "
                << " solving linear system: start" << endl;
        }
        else
        {
            Info<< "    Solving" << endl;
        }

        if (Switch(solidModelDict().lookup("usePETSc"))) //Check if PETSCs is defined
        {
#ifdef USE_PETSC
            //If defined solve linear system using PETSc
            fileName optionsFile(solidModelDict().lookup("optionsFile"));
            solverPerf = sparseMatrixTools::solveLinearSystemPETSc
            (
                matrix,
                source,
                pointDcorr,
                twoD_,
                optionsFile,
                mesh().points(),
                ownedByThisProc,
                localToGlobalPointMap,
                globalPointIndices_.stencilSizeOwned(),
                globalPointIndices_.stencilSizeNotOwned(),
                solidModelDict().lookupOrDefault<bool>("debugPETSc", false)
            );
#else
            FatalErrorIn("vertexCentredLinGeomSolidMMS::evolve()")
                << "PETSc not available. Please set the PETSC_DIR environment "
                << "variable and re-compile solids4foam" << abort(FatalError);
#endif
        }
        else
        {
            // Otherwise, use Eigen SparseLU direct solver
            sparseMatrixTools::solveLinearSystemEigen
            (
                matrix, source, pointDcorr, twoD_, false, debug
            );
        }

        if (debug)
        {
            Info<< "bool vertexCentredLinGeomSolidMMS::evolve(): "
                << " solving linear system: end" << endl;
        }

        // Update point displacement field
        if (Switch(solidModelDict().lookup("lineSearch")))
        {
            // Lookup target tolerance for slope reduction
            const scalar rTol
            (
                solidModelDict().lookupOrDefault<scalar>("lineSearchRTol", 0.8)
            );

            // Lookup the maximum number of line search iterations
            const int maxIter
            (
                solidModelDict().lookupOrDefault<scalar>
                (
                    "lineSearchMaxIter", 10
                )
            );

            // Calculate line search factor
            const scalar eta
            (
                calculateLineSearchFactor
                (
                    rTol, maxIter, pointDcorr, source, zeta
                )
            );

            // Update displacement field
#ifdef OPENFOAMESIORFOUNDATION
            pointD().primitiveFieldRef() += eta*pointDcorr;
#else
            pointD().internalField() += eta*pointDcorr;
#endif
        }
#ifdef OPENFOAMESIORFOUNDATION
        else if (mesh().relaxField(pointD().name()))
#else
        else if (mesh().solutionDict().relaxField(pointD().name()))
#endif
        {
            // Relaxing the correction can help convergence
            // This is like a simple line search
            // Of course, an actual line search would be better: to-do!

#ifdef OPENFOAMESIORFOUNDATION
            const scalar rf
            (
                mesh().fieldRelaxationFactor(pointD().name())
            );

            pointD().primitiveFieldRef() += rf*pointDcorr;
#else
            const scalar rf
            (
                mesh().solutionDict().fieldRelaxationFactor(pointD().name())
            );

            pointD().internalField() += rf*pointDcorr;
#endif
        }
        else
        {
#ifdef OPENFOAMESIORFOUNDATION
            pointD().primitiveFieldRef() += pointDcorr;
#else
            pointD().internalField() += pointDcorr;
#endif
        }
        pointD().correctBoundaryConditions();

        // Update point accelerations
        // Note: for NewmarkBeta, this needs to come before the pointU update
#ifdef OPENFOAMESIORFOUNDATION
        pointA_.primitiveFieldRef() =
            vfvc::ddt
            (
                mesh().ddtScheme("ddt(pointU)"),
                mesh().d2dt2Scheme("d2dt2(pointD)"),
                pointU_
            );

        // Update point velocities
        pointU_.primitiveFieldRef() =
            vfvc::ddt
            (
                mesh().ddtScheme("ddt(pointD)"),
                mesh().d2dt2Scheme("d2dt2(pointD)"),
                pointD()
            );
#else
        pointA_.internalField() =
            vfvc::ddt
            (
                mesh().schemesDict().ddtScheme("ddt(pointU)"),
                mesh().schemesDict().d2dt2Scheme("d2dt2(pointD)"),
                pointU_
            );

        // Update point velocities
        pointU_.internalField() =
            vfvc::ddt
            (
                mesh().schemesDict().ddtScheme("ddt(pointD)"),
                mesh().schemesDict().d2dt2Scheme("d2dt2(pointD)"),
                pointD()
            );
#endif
    }
    while
    (
        !converged
        (
            iCorr,
            initResidual,
            solverPerf.finalResidual()[vector::X],
#ifdef OPENFOAMESIORFOUNDATION
            cmptMax(solverPerf.nIterations()),
#else
            solverPerf.nIterations(),
#endif
            pointD(),
            pointDcorr
        ) && ++iCorr
    );

    // Calculate gradD at dual faces
    dualGradDf_ = vfvc::fGrad
    (
        pointD(),
        mesh(),
        dualMesh(),
        dualMeshMap().dualFaceToCell(),
        dualMeshMap().dualCellToPoint(),
        zeta,
        debug
    );

    // Update the increment of displacement
    pointDD() = pointD() - pointD().oldTime();
    
    //Info << pointD().internalField() << endl;

    // Calculate cell gradient
    // This assumes a constant gradient within each primary mesh cell
    gradD() = vfvc::grad(pointD(), mesh());

    // Map primary cell gradD field to sub-meshes for multi-material cases
    if (mechanical().PtrList<mechanicalLaw>::size() > 1) //Meaning??
    {
        mechanical().mapGradToSubMeshes(gradD()); 
    }

    // Update dual face stress field
    dualMechanicalPtr_().correct(dualSigmaf_);

    // Update primary mesh cell stress field, assuming it is constant per
    // primary mesh cell
    mechanical().correct(sigma());

    return true;
}

void vertexCentredLinGeomSolidMMS::setTraction
(
    const label interfaceI,
    const label patchID,
    const vectorField& faceZoneTraction
)
{
    // Get point field on patch
    const vectorField traction
    (
        globalPatches()[interfaceI].globalPointToPatch
        (
            globalPatches()[interfaceI].interpolator().faceToPointInterpolate
            (
                faceZoneTraction
            )()
        )
    );

    // Lookup point patch field
#ifdef OPENFOAMESIORFOUNDATION
    pointPatchVectorField& ptPatch = pointD().boundaryFieldRef()[patchID];
#else
    pointPatchVectorField& ptPatch = pointD().boundaryField()[patchID];
#endif

    if (isA<solidTractionPointPatchVectorField>(ptPatch))
    {
        solidTractionPointPatchVectorField& patchD =
            refCast<solidTractionPointPatchVectorField>(ptPatch);

        patchD.traction() = traction;
    }
    else
    {
        FatalErrorIn
        (
            "void Foam::vertexCentredLinGeomSolidMMS::setTraction\n"
            "(\n"
            "    fvPatchVectorField& tractionPatch,\n"
            "    const vectorField& traction\n"
            ")"
        )   << "Boundary condition "
            << ptPatch.type()
            << " for point patch " << ptPatch.patch().name()
            << " should instead be type "
            << solidTractionPointPatchVectorField::typeName
            << abort(FatalError);
    }
}

void vertexCentredLinGeomSolidMMS::writeFields(const Time& runTime)
{
    // Calculate gradD at the primary points using least squares: this should
    // be second-order accurate (... I think).
    const pointTensorField pGradD(vfvc::pGrad(pointD(), mesh()));

    // Calculate strain at the primary points based on pGradD
    // Note: the symm operator is not defined for pointTensorFields so we will
    // do it manually
    // const pointSymmTensorField pEpsilon("pEpsilon", symm(pGradD));
    pointSymmTensorField pEpsilon
    (
        IOobject
        (
            "pEpsilon",
            runTime.timeName(),
            runTime,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        pMesh(),
        dimensionedSymmTensor("0", dimless, symmTensor::zero)
    );

#ifdef FOAMEXTEND
    pEpsilon.internalField() = symm(pGradD.internalField());
#else
    pEpsilon.primitiveFieldRef() = symm(pGradD.internalField());
#endif
    pEpsilon.write();

    // Equivalent strain at the points
    pointScalarField pEpsilonEq
    (
        IOobject
        (
            "pEpsilonEq",
            runTime.timeName(),
            runTime,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        pMesh(),
        dimensionedScalar("0", dimless, 0.0)
    );

#ifdef FOAMEXTEND
    pEpsilonEq.internalField() =
        sqrt((2.0/3.0)*magSqr(dev(pEpsilon.internalField())));
#else
    pEpsilonEq.primitiveFieldRef() =
        sqrt((2.0/3.0)*magSqr(dev(pEpsilon.internalField())));
#endif
    pEpsilonEq.write();

    Info<< "Max pEpsilonEq = " << gMax(pEpsilonEq) << endl;

    solidModel::writeFields(runTime);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solidModels

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //Traction
