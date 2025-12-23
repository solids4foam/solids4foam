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

#ifdef USE_PETSC

#include <functional>
#include "hofvc.H"
#include "fvc.H"
#include "multiplyCoeff.H"
#include "multiplyCoeffExtended.H"
#include "sparseMatrixTools.H"
#include "cellPointLeastSquaresVectors.H"
#include "surfaceFields.H"
#include "processorPolyPatch.H"
#include "emptyPolyPatch.H"
#include "symmetryPolyPatch.H"
#include "fixedDisplacementFvPatchVectorField.H"
#include "solidTractionFvPatchVectorField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

using namespace Foam;

tmp<surfaceVectorField> hofvc::surfaceIntegrate
(
    const List<List<symmTensor>>& quadSigma,
    const CompactListList<scalar>& quadW,
    const fvMesh& mesh
)
{
    tmp<surfaceVectorField> tsf
    (
        new surfaceVectorField
        (
            IOobject
            (
                "surfaceIntegrate(sigma)",
                 mesh.time().timeName(),
                 mesh,
                 IOobject::NO_READ,
                 IOobject::NO_WRITE
            ),
            mesh,
            dimensionedVector("0", dimPressure, vector::zero)
        )
    );

    surfaceVectorField& tf = tsf.ref();

    const vectorField normal(mesh.faceAreas()/mag(mesh.faceAreas()));

    forAll (tf, faceI)
    {
        // Sigma at the quadrature points on the face
        const List<symmTensor>& faceQuadStress = quadSigma[faceI];

        const vector& faceNormal = normal[faceI];

        forAll(faceQuadStress, pI)
        {
            // Add traction contribution of this quadrature point
            tf[faceI] += faceNormal & (faceQuadStress[pI] * quadW[faceI][pI]);
        }
    }

    forAll (tf.boundaryField(), patchI)
    {
        vectorField& tfPatch = tf.boundaryFieldRef()[patchI];
        forAll(tfPatch, faceI)
        {
            const label globalFaceID =
                mesh.boundaryMesh()[patchI].start() + faceI;

            // Sigma at the Gauss quadrature points on the face
            const List<symmTensor>& faceQuadStress = quadSigma[globalFaceID];

            const vector& faceNormal =  normal[globalFaceID];

            forAll(faceQuadStress, pI)
            {
                // Add traction contribution of this quadrature point
                tfPatch[faceI] +=
                    faceNormal & (faceQuadStress[pI] * quadW[faceI][pI]);
            }
        }
    }

    return tsf;
};


tmp<surfaceVectorField> hofvc::surfaceIntegrate
(
    const List<List<symmTensor>>& quadSigma,
    const List<List<tensor>>& gradDQuad,
    const CompactListList<scalar>& quadW,
    const fvMesh& mesh
)
{
    tmp<surfaceVectorField> tsf
    (
        new surfaceVectorField
        (
            IOobject
            (
                "surfaceIntegrate(sigma)",
                 mesh.time().timeName(),
                 mesh,
                 IOobject::NO_READ,
                 IOobject::NO_WRITE
            ),
            mesh,
            dimensionedVector("0", dimPressure, vector::zero)
        )
    );

    surfaceVectorField& tf = tsf.ref();

    const vectorField normal(mesh.faceAreas()/mag(mesh.faceAreas()));

    forAll (tf, faceI)
    {
        // Sigma at the quadrature points on the face
        const List<symmTensor>& faceQuadSigma = quadSigma[faceI];

        const vector& faceNormal = normal[faceI];

        forAll(faceQuadSigma, pI)
        {
            const tensor& gradD = gradDQuad[faceI][pI];
            const tensor F = I + gradD.T();
            const tensor invFT = inv(F.T());
            const scalar J = det(F);

            // Add traction contribution of this quadrature point
            tf[faceI] +=
                quadW[faceI][pI] * J * ((faceQuadSigma[pI] & invFT) & faceNormal);
        }
    }

    forAll (tf.boundaryField(), patchI)
    {
        vectorField& tfPatch = tf.boundaryFieldRef()[patchI];
        forAll(tfPatch, faceI)
        {
            const label globalFaceID =
                mesh.boundaryMesh()[patchI].start() + faceI;

            // Sigma at the Gauss quadrature points on the face
            const List<symmTensor>& faceQuadSigma = quadSigma[globalFaceID];

            // Grad of displacement for the quadrature points on the face
            const List<tensor>& faceQuadGradD = gradDQuad[globalFaceID];

            const vector& faceNormal = normal[globalFaceID];

            forAll(faceQuadSigma, pI)
            {
                const tensor& gradD = faceQuadGradD[pI];
                const tensor F = I + gradD.T();
                const tensor invFT = inv(F.T());
                const scalar J = det(F);

                // Add traction contribution of this quadrature point
                tfPatch[faceI] +=
                    quadW[faceI][pI] * J * ((faceQuadSigma[pI] & invFT) & faceNormal);
            }
        }
    }

    return tsf;
};


autoPtr<CompactListList<vector>> hofvc::cellToQuad
(
    const volVectorField& D,
    const LRE& LREInterp
)
{
    const fvMesh& mesh = LREInterp.mesh();
    const bool twoD = (mesh.nGeometricD() == 2);

    // Cell quadrature points
    const CompactListList<point>& quadPts = LREInterp.cellQuadPoints();

    labelList nQpPerCell(mesh.nCells(), 0);
    forAll(nQpPerCell, cellI)
    {
	nQpPerCell[cellI] = quadPts[cellI].size();
    }

    autoPtr<CompactListList<vector>> DQuadPtr;
    DQuadPtr.set(new CompactListList<vector>(nQpPerCell));

    CompactListList<vector>& DQuad  = DQuadPtr();

    const vectorField& C = mesh.C();

    tmp<volTensorField> tGradD = LREInterp.grad(D);
    const volTensorField& gradD = tGradD();
    const tensorField& gradDI = gradD.internalField();

    // Decompose the vector field D into its scalar components
    volScalarField Dx(D.component(vector::X));
    volScalarField Dy(D.component(vector::Y));
    volScalarField Dz(D.component(vector::Z));

    tmp<volSymmTensorField> tHessDx = LREInterp.hessian(Dx);
    tmp<volSymmTensorField> tHessDy = LREInterp.hessian(Dy);
    tmp<volSymmTensorField> tHessDz = LREInterp.hessian(Dz);

    const symmTensorField& hessDxI = tHessDx->internalField();
    const symmTensorField& hessDyI = tHessDy->internalField();
    const symmTensorField& hessDzI = tHessDz->internalField();

    autoPtr<List<LRE::symmTensor3Order>> tThirdDerDx;
    autoPtr<List<LRE::symmTensor3Order>> tThirdDerDy;
    autoPtr<List<LRE::symmTensor3Order>> tThirdDerDz;
    if (LREInterp.order() >= 3)
    {
	tThirdDerDx = LREInterp.thirdDeriv(Dx);
	tThirdDerDy = LREInterp.thirdDeriv(Dy);
	tThirdDerDz = LREInterp.thirdDeriv(Dz);
    }

    // Displacement at quadrature points
    forAll(DQuad, cellI)
    {
	const vector& cellCentreD = D[cellI];
	const vector& centre = C[cellI];

	// Loop over quadrature points
	forAll(DQuad[cellI], ptI)
	{
	    // Distance from quadrature point to the cell centre
	    const vector d = quadPts[cellI][ptI] - centre;

	    // Set to zero
	    DQuad[cellI][ptI] = vector::zero;

	    // Linear part
	    DQuad[cellI][ptI] += cellCentreD + (d & gradDI[cellI]);

	    // Quadratic part
	    if (LREInterp.order() >= 2)
	    {
		DQuad[cellI][ptI].x() += 0.5 * (d & (hessDxI[cellI] & d));
		DQuad[cellI][ptI].y() += 0.5 * (d & (hessDyI[cellI] & d));
		DQuad[cellI][ptI].z() += 0.5 * (d & (hessDzI[cellI] & d));
	    }

	    // Cubic part
	    if (LREInterp.order() >= 3)
	    {
		DQuad[cellI][ptI].x() += (1.0/6.0) * LRE::cubicForm((*tThirdDerDx)[cellI], d, twoD);
		DQuad[cellI][ptI].y() += (1.0/6.0) * LRE::cubicForm((*tThirdDerDy)[cellI], d, twoD);
		DQuad[cellI][ptI].z() += (1.0/6.0) * LRE::cubicForm((*tThirdDerDz)[cellI], d, twoD);
	    }
	}
    }

    return DQuadPtr;
}


autoPtr<CompactListList<vector>> hofvc::ddt
(
    const volVectorField& D,
    const LRE& LREInterp
)
{
    const fvMesh& mesh = LREInterp.mesh();

    //Check what scheme is prescribed, continue only in the case of backward
    word schemeName;
    (mesh.ddtSchemes().found(D.name())
      ? mesh.ddtSchemes().lookup(D.name())
      : mesh.ddtSchemes().lookup("default")) >> schemeName;

    if (schemeName != "backward")
    {
	FatalErrorInFunction
            << "Unsupported ddt scheme" << abort(FatalError);
    }

    // Cell quadrature points
    const CompactListList<scalar>& quadW = LREInterp.cellQuadWeight();

    // Construct acceleration field
    labelList nQpPerCell(mesh.nCells(), 0);
    forAll(nQpPerCell, cellI)
    {
	nQpPerCell[cellI] = quadW[cellI].size();
    }

    autoPtr<CompactListList<vector>> velocityPtr;
    velocityPtr.set(new CompactListList<vector>(nQpPerCell));

    CompactListList<vector>& velocity  = velocityPtr();

    // Extrapolate D, D0, D00 to quad points
    autoPtr<CompactListList<vector>> DQuadPtr = hofvc::cellToQuad(D, LREInterp);
    autoPtr<CompactListList<vector>> DQuad0Ptr = hofvc::cellToQuad(D.oldTime(), LREInterp);
    autoPtr<CompactListList<vector>> DQuad00Ptr = hofvc::cellToQuad(D.oldTime().oldTime(), LREInterp);

    const CompactListList<vector>& DQuad = DQuadPtr();
    const CompactListList<vector>& DQuad0 = DQuad0Ptr();
    const CompactListList<vector>& DQuad00 = DQuad00Ptr();

    const scalar rDeltaT = 1.0/mesh.time().deltaT().value();

    const bool oldOldTime =
	(D.oldTime().timeIndex() == D.oldTime().oldTime().timeIndex());
    const scalar deltaT = mesh.time().deltaT().value();
    const scalar deltaT0 = oldOldTime ? GREAT : mesh.time().deltaT0().value();

    const scalar coefft = 1 + deltaT/(deltaT + deltaT0);
    const scalar coefft00 = deltaT*deltaT/(deltaT0*(deltaT + deltaT0));
    const scalar coefft0 = coefft + coefft00;

    // Velocity at quadrature points, BDF2 scheme
    forAll(velocity, cellI)
    {
	// Loop over quadrature points
	forAll(velocity[cellI], ptI)
	{
	    velocity[cellI][ptI] =
   		rDeltaT*
		(
		    coefft*DQuad[cellI][ptI]
		  - coefft0*DQuad0[cellI][ptI]
		  + coefft00*DQuad00[cellI][ptI]
	        );
	}
    }

    return velocityPtr;
}

tmp<volVectorField> hofvc::d2dt2
(
    const volVectorField& D,
    const LRE& LREInterp
)
{
    // Hard-coded BDF2 scheme

    const fvMesh& mesh = LREInterp.mesh();

    tmp<volVectorField> tvf
    (
        new volVectorField
        (
            IOobject
            (
	         "d2dt2(" + D.name() + ")",
                 mesh.time().timeName(),
                 mesh,
                 IOobject::NO_READ,
                 IOobject::NO_WRITE
            ),
            mesh,
	    dimensionedVector
	    (
	        "zero",
		D.dimensions()/(dimTime*dimTime),
		vector::zero
	    )
	 )
    );

    volVectorField& tf = tvf.ref();

    //Check what scheme is prescribed, continue only in the case of backward
    word schemeName;
    (mesh.d2dt2Schemes().found(D.name())
      ? mesh.d2dt2Schemes().lookup(D.name())
      : mesh.d2dt2Schemes().lookup("default")) >> schemeName;

    if (schemeName != "backward")
    {
	return tvf;
    }

    // Cell quadrature points weights
    const CompactListList<scalar>& quadW = LREInterp.cellQuadWeight();

    // Construct acceleration field
    labelList nQpPerCell(mesh.nCells(), 0);
    forAll(nQpPerCell, cellI)
    {
	nQpPerCell[cellI] = quadW[cellI].size();
    }
    CompactListList<vector> a(nQpPerCell);

    autoPtr<CompactListList<vector>> vPtr = hofvc::ddt(D, LREInterp);
    autoPtr<CompactListList<vector>> v0Ptr = hofvc::ddt(D.oldTime(), LREInterp);
    autoPtr<CompactListList<vector>> v00Ptr = hofvc::ddt(D.oldTime().oldTime(), LREInterp);

    const CompactListList<vector>& v = vPtr();
    const CompactListList<vector>& v0 = v0Ptr();
    const CompactListList<vector>& v00 = v00Ptr();

    const scalar rDeltaT = 1.0/mesh.time().deltaT().value();

    const bool oldOldTime =
	(D.oldTime().timeIndex() == D.oldTime().oldTime().timeIndex());
    const scalar deltaT = mesh.time().deltaT().value();
    const scalar deltaT0 = oldOldTime ? GREAT : mesh.time().deltaT0().value();

    const scalar coefft = 1 + deltaT/(deltaT + deltaT0);
    const scalar coefft00 = deltaT*deltaT/(deltaT0*(deltaT + deltaT0));
    const scalar coefft0 = coefft + coefft00;

    // Acceleration at quadrature points, BDF2 scheme
    forAll(a, cellI)
    {
	// Loop over quadrature points
	forAll(a[cellI], ptI)
	{
	    a[cellI][ptI] =
		rDeltaT*
		(
		    coefft*v[cellI][ptI]
		  - coefft0*v0[cellI][ptI]
		  + coefft00*v00[cellI][ptI]
	        );
	}
    }

    // Integrate inertia term
    forAll(tf, cellI)
    {
	forAll(a[cellI], ptI)
	{
	    tf[cellI] += a[cellI][ptI] * quadW[cellI][ptI];
	}
    }

    return tvf;
}

// ************************************************************************* //

#endif // #ifdef USE_PETSC
