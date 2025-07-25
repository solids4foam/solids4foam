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
#include "hofvm.H"
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


Foam::tensor Foam::hofvm::laplacianCoeff
(
    const scalar& gammaMagSf,
    const scalar& quadWeight,
    const vector& gradInterpCoeff,
    const vector& faceNormal
)
{
    return I * gammaMagSf * quadWeight * (gradInterpCoeff & faceNormal);
}


Foam::tensor Foam::hofvm::laplacianTransposeCoeff
(
    const scalar& gammaMagSf,
    const scalar& quadWeight,
    const vector& gradInterpCoeff,
    const vector& faceNormal
)
{
    const vector& c = gradInterpCoeff;
    const vector& n = faceNormal;

    const tensor t = tensor
	(
	    c.x()*n.x(), c.x()*n.y(), c.x()*n.z(),
	    c.y()*n.x(), c.y()*n.y(), c.y()*n.z(),
	    c.z()*n.x(), c.z()*n.y(), c.z()*n.z()
	);

    return gammaMagSf * quadWeight * t;
}


Foam::tensor Foam::hofvm::laplacianTraceCoeff
(
    const scalar& gammaMagSf,
    const scalar& quadWeight,
    const vector& gradInterpCoeff,
    const vector& faceNormal
)
{
    const vector& c = gradInterpCoeff;
    const vector& n = faceNormal;

    const tensor t = tensor
	(
	    c.x()*n.x(), c.y()*n.x(), c.z()*n.x(),
	    c.x()*n.y(), c.y()*n.y(), c.z()*n.y(),
	    c.x()*n.z(), c.y()*n.z(), c.z()*n.z()
	);

    return gammaMagSf * quadWeight * t;
}


void Foam::hofvm::laplacianIntoPETScMatrix
(
    Mat matrix,
    const foamPetscSnesHelper& petscSnesHelper,
    const fvMesh& mesh,
    const volVectorField& D,
    const volScalarField& diffusivity,
    const LRE& lre,
    const label nScalarEqns,
    const label rowOffset,
    const label colOffset
)
{
    hofvmLaplacianPETSc
    (
        matrix,
	petscSnesHelper,
        mesh,
	D,
        diffusivity,
        lre,
        hofvm::laplacianCoeff,
        rowOffset,
	colOffset,
	nScalarEqns
    );
}


void Foam::hofvm::laplacianIntoSparseMatrix
(
    sparseMatrix& matrix,
    vectorField& source,
    const fvMesh& mesh,
    const volVectorField& D,
    const volScalarField& diffusivity,
    const LRE& lre
)
{
    hofvmLaplacianSparseMatrix
    (
        matrix,
        source,
        mesh,
	D,
        diffusivity,
        lre,
        hofvm::laplacianCoeff
    );
}


void Foam::hofvm::laplacianTraceIntoPETScMatrix
(
    Mat matrix,
    const foamPetscSnesHelper& petscSnesHelper,
    const fvMesh& mesh,
    const volVectorField& D,
    const volScalarField& diffusivity,
    const LRE& lre,
    const label nScalarEqns,
    const label rowOffset,
    const label colOffset
)
{
    hofvmLaplacianPETSc
    (
        matrix,
	petscSnesHelper,
        mesh,
	D,
	diffusivity,
        lre,
        hofvm::laplacianTraceCoeff,
        rowOffset,
	colOffset,
	nScalarEqns
    );
}


void Foam::hofvm::laplacianTraceIntoSparseMatrix
(
    sparseMatrix& matrix,
    vectorField& source,
    const fvMesh& mesh,
    const volVectorField& D,
    const volScalarField& diffusivity,
    const LRE& lre
)
{
    hofvmLaplacianSparseMatrix
    (
        matrix,
        source,
        mesh,
	D,
        diffusivity,
        lre,
        hofvm::laplacianTraceCoeff
    );
}


void Foam::hofvm::laplacianTransposeIntoPETScMatrix
(
    Mat matrix,
    const foamPetscSnesHelper& petscSnesHelper,
    const fvMesh& mesh,
    const volVectorField& D,
    const volScalarField& diffusivity,
    const LRE& lre,
    const label nScalarEqns,
    const label rowOffset,
    const label colOffset
)
{
    hofvmLaplacianPETSc
    (
        matrix,
	petscSnesHelper,
        mesh,
	D,
        diffusivity,
        lre,
        hofvm::laplacianTransposeCoeff,
        rowOffset,
	colOffset,
	nScalarEqns
    );
}


void Foam::hofvm::laplacianTransposeIntoSparseMatrix
(
    sparseMatrix& matrix,
    vectorField& source,
    const fvMesh& mesh,
    const volVectorField& D,
    const volScalarField& diffusivity,
    const LRE& lre
)
{
    hofvmLaplacianSparseMatrix
    (
        matrix,
        source,
        mesh,
	D,
        diffusivity,
        lre,
        hofvm::laplacianTransposeCoeff
    );
}


Foam::label Foam::hofvm::hofvmLaplacianPETSc
(
    Mat matrix,
    const foamPetscSnesHelper& petscSnesHelper,
    const fvMesh& mesh,
    const volVectorField& D,
    const volScalarField& diffusivity,
    const LRE& lre,
    tensor (*calcCoeff)
    (
        const scalar& gammaMagSf,
	const scalar& quadWeight,
	const vector& gradInterpCoeff,
	const vector& faceNormal
    ),
    const label rowOffset,
    const label colOffset,
    const label nScalarEqns
)
{
    const labelUList& owner = mesh.owner();
    const labelUList& neighbour = mesh.neighbour();
    const scalarField& magSfI = mesh.magSf().internalField();
    const surfaceVectorField n(mesh.Sf()/mesh.magSf());

    // Diffusion coefficient linearly interpolated to face centres.
    // TO_DO: interpolate diffusivity to quad points using hofvc::interpolate
    const surfaceScalarField gamma(fvc::interpolate(diffusivity));
    const scalarField& gammaI = gamma.internalField();

    // Face quadrature points weights
    const List<List<scalar>>& facesQuadWeights = lre.faceGaussPointsWeight();

    // Faces stencil
    const List<labelList>& stencils = lre.globalFaceStencils();

    // Gradient interpolation coefficients
    const List<List<DynamicList<vector>>>& gradCoeffs =
	lre.QRGradFaceGPCoeffs();

    // Get the blockSize
    label blockSize;
    MatGetBlockSize(matrix, &blockSize);

    // Initialise the block coefficient
    const label nCoeffCmpts = blockSize*blockSize;
    List<PetscScalar> values(nCoeffCmpts, 0.0);

    // Loop over internal faces
    forAll(owner, faceI)
    {
	// Preliminaries
	const vector& faceNormal = n[faceI];

	const scalar gammaMagSf = magSfI[faceI] * gammaI[faceI];

	// Face interpolation molecule
	const labelList& faceStencil = stencils[faceI];

	// Face quadrature points weights
	const List<scalar>& faceQuadWeight = facesQuadWeights[faceI];

	// Loop over face quadrature points
	forAll(faceQuadWeight, pointI)
	{
	    // Quad point weight
	    const scalar& quadPointW = faceQuadWeight[pointI];

	    // Loop over interpolation stencil
	    for(label cI = 0; cI < faceStencil.size(); cI++)
	    {
		const label globalCellID = faceStencil[cI];
		const vector& cellGradCoeff = gradCoeffs[faceI][pointI][cI];

		const tensor coeff =
		    calcCoeff(gammaMagSf, quadPointW, cellGradCoeff, faceNormal);

		// Local block row ID
		const label ownCellID = owner[faceI];

		// Local block column ID
		const label neiCellID = neighbour[faceI];

		// Global block row ID
		const label globalBlockRowI =
		    petscSnesHelper.globalCells().toGlobal(ownCellID);

		// Global block column ID
		const label globalBlockColI =
		    petscSnesHelper.globalCells().toGlobal(neiCellID);

		for (label i = 0; i < (blockSize - 1); ++i)
		{
                    for (label j = 0; j < (blockSize - 1); ++j)
		    {
			// Copy 3x3 (or 2x2 in 2-D) coeff into the top left of
			// the 4x4 (or 3x3 in 2-D) values matrix
			values[(i + rowOffset)*blockSize + j + colOffset] =
			    coeff[i*3 + j];
		    }
		}

		CHKERRQ
		(
		    MatSetValuesBlocked
		    (
		        matrix, 1, &globalBlockRowI, 1, &globalCellID,
			values.cdata(),
			ADD_VALUES
		    )
		);

		for (label i = 0; i < nCoeffCmpts; ++i)
		{
		    values[i] = -values[i];
		}

		CHKERRQ
		(
		    MatSetValuesBlocked
		    (
		        matrix, 1, &globalBlockColI, 1, &globalCellID,
			values.cdata(),
			ADD_VALUES
		    )
		);
	    }
	}
    }

    // Loop over boundary faces
    forAll(mesh.boundaryMesh(), patchI)
    {
	const word& patchType = mesh.boundaryMesh()[patchI].type();
	const scalarField& pMagSf = mesh.magSf().boundaryField()[patchI];
	const scalarField& pGamma = gamma.boundaryField()[patchI];
	const vectorField pNormal(mesh.boundary()[patchI].nf());

	const label start = mesh.boundaryMesh()[patchI].start();

	if (patchType == emptyPolyPatch::typeName)
	{
	    // Skip empty patches
	}
	else if (patchType == processorPolyPatch::typeName)
	{
	    NotImplemented;
	}
	else if (patchType == symmetryPolyPatch::typeName)
	{
	    NotImplemented;
	}
	else if
	(
	    isA<fixedGradientFvPatchVectorField>(D.boundaryField()[patchI])
	)
        {
	    // For now, skip faces with prescribed gradient.
	    // Contribution is added directly to the source vector in the case
	    // of traction
	}
	else if
	(
	    isA<fixedValueFvPatchVectorField>(D.boundaryField()[patchI])
	)
	{
            forAll(mesh.boundaryMesh()[patchI], faceI)
            {
		// Preliminaries
		const vector& faceNormal = pNormal[faceI];
		const scalar gammaMagSf = pMagSf[faceI] * pGamma[faceI];

		// Get global face index, needed for lists from LRE class
		const label faceID = faceI + start;

		// Face interpolation molecule
		const labelList& faceStencil = stencils[faceID];

		// Face quadrature points weights
		const List<scalar>& faceQuadWeight = facesQuadWeights[faceID];

		forAll(faceQuadWeight, pointI)
		{
		    // Quad point weight
		    const scalar& quadPointW = faceQuadWeight[pointI];

		    // Loop over interpolation stencil.
		    for(label cI = 0; cI < faceStencil.size(); cI++)
		    {
			const label globalCellID = faceStencil[cI];
			const vector& cellGradCoeff =
			    gradCoeffs[faceID][pointI][cI];

			const tensor coeff =
			    calcCoeff
			    (
			        gammaMagSf,
				quadPointW,
				cellGradCoeff,
				faceNormal
			    );

			for (label i = 0; i < (blockSize - 1); ++i)
			{
			    for (label j = 0; j < (blockSize - 1); ++j)
			    {
				// Copy 3x3 (or 2x2 in 2-D) coeff into the top left of
				// the 4x4 (or 3x3 in 2-D) values matrix
				values[(i + rowOffset)*blockSize + j + colOffset] =
				    coeff[i*3 + j];
			    }
			}

			// Local block row ID
			const label ownCellID = owner[faceI];

			// Global block row ID
			const label globalBlockRowI =
			    petscSnesHelper.globalCells().toGlobal(ownCellID);

			CHKERRQ
			(
			    MatSetValuesBlocked
			    (
			        matrix, 1, &globalBlockRowI, 1, &globalCellID,
				values.cdata(),
				ADD_VALUES
			     )
			 );
		    }

		    // Include boundary value to the source
		    // Last item in gradCoeff refers to boundary face.
		    {
		        const label size = faceStencil.size();

		        const vector& cellGradCoeff =
		            gradCoeffs[faceID][pointI][size];

		        const tensor coeff =
		            calcCoeff
		            (
		                gammaMagSf,
		        	quadPointW,
		        	cellGradCoeff,
		        	faceNormal
		            );

                        //source[owner[faceID]] -=
			//coeff & D.boundaryField()[patchI][faceI];
		    }
 	        }
	    }
	}
	else
	{
	    // Currently only displacement and traction boundary are implemented
	    NotImplemented;
	}
    }

    return 0;
}

void Foam::hofvm::hofvmLaplacianSparseMatrix
(
    sparseMatrix& matrix,
    vectorField& source,
    const fvMesh& mesh,
    const volVectorField& D,
    const volScalarField& diffusivity,
    const LRE& lre,
    tensor (*calcCoeff)
    (
        const scalar& gammaMagSf,
	const scalar& quadWeight,
	const vector& gradInterpCoeff,
	const vector& faceNormal
    )
)
{
    const labelUList& owner = mesh.owner();
    const labelUList& neighbour = mesh.neighbour();
    const scalarField& magSfI = mesh.magSf().internalField();
    const surfaceVectorField n(mesh.Sf()/mesh.magSf());

    // Diffusion coefficient linearly interpolated to face centres.
    // TO_DO: interpolate diffusivity to quad points using hofvc::interpolate
    const surfaceScalarField gamma(fvc::interpolate(diffusivity));
    const scalarField& gammaI = gamma.internalField();

    // Face quadrature points weights
    const List<List<scalar>>& facesQuadWeights = lre.faceGaussPointsWeight();

    // Face quadrature points
    const List<List<point>>& faceQuadPoints = lre.faceGaussPoints();

    // Faces stencil
    const List<labelList>& stencils = lre.globalFaceStencils();

    // Gradient interpolation coefficients
    const List<List<DynamicList<vector>>>& gradCoeffs =
	lre.QRGradFaceGPCoeffs();

    // Loop over internal faces
    forAll(owner, faceI)
    {
	// Preliminaries
	const vector& faceNormal = n[faceI];

	const scalar gammaMagSf = magSfI[faceI] * gammaI[faceI];

	// Face interpolation molecule
	const labelList& faceStencil = stencils[faceI];

	// Face quadrature points weights
	const List<scalar>& faceQuadWeight = facesQuadWeights[faceI];

	// Loop over face quadrature points
	forAll(faceQuadWeight, pointI)
	{
	    // Quad point weight
	    const scalar& quadPointW = faceQuadWeight[pointI];

	    // Loop over interpolation stencil
	    for(label cI = 0; cI < faceStencil.size(); cI++)
	    {
		const label globalCellID = faceStencil[cI];
		const vector& cellGradCoeff = gradCoeffs[faceI][pointI][cI];

		const tensor coeff =
		    calcCoeff(gammaMagSf, quadPointW, cellGradCoeff, faceNormal);

		matrix(owner[faceI], globalCellID) += coeff;
		matrix(neighbour[faceI], globalCellID) -= coeff;
	    }
	}
    }

    // Loop over boundary faces
    forAll(mesh.boundaryMesh(), patchI)
    {
	const word& patchType = mesh.boundaryMesh()[patchI].type();
	const scalarField& pMagSf = mesh.magSf().boundaryField()[patchI];
	const scalarField& pGamma = gamma.boundaryField()[patchI];
	const vectorField pNormal(mesh.boundary()[patchI].nf());

	const label start = mesh.boundaryMesh()[patchI].start();

	if (patchType == emptyPolyPatch::typeName)
	{
	    // Skip empty patches
	}
	else if (patchType == processorPolyPatch::typeName)
	{
	    NotImplemented;
	}
	else if (patchType == symmetryPolyPatch::typeName)
	{
	    NotImplemented;
	}
	else if
	(
	    isA<fixedGradientFvPatchVectorField>(D.boundaryField()[patchI])
	)
        {
	    // For now, skip faces with prescribed gradient.
	    // Contribution is added directly to the source vector in the case
	    // of traction
	}
	else if
	(
	    isA<fixedValueFvPatchVectorField>(D.boundaryField()[patchI])
	)
	{
            forAll(mesh.boundaryMesh()[patchI], faceI)
            {
		// Preliminaries
		const vector& faceNormal = pNormal[faceI];
		const scalar gammaMagSf = pMagSf[faceI] * pGamma[faceI];

		// Get global face index, needed for lists from LRE class
		const label faceID = faceI + start;

		// Face interpolation molecule
		const labelList& faceStencil = stencils[faceID];

		// Face quadrature points weights
		const List<scalar>& faceQuadWeight = facesQuadWeights[faceID];

		forAll(faceQuadWeight, pointI)
		{
		    // Quad point weight
		    const scalar& quadPointW = faceQuadWeight[pointI];

		    // Loop over interpolation stencil.
		    for(label cI = 0; cI < faceStencil.size(); cI++)
		    {
			const label globalCellID = faceStencil[cI];
			const vector& cellGradCoeff =
			    gradCoeffs[faceID][pointI][cI];

			const tensor coeff =
			    calcCoeff
			    (
			        gammaMagSf,
				quadPointW,
				cellGradCoeff,
				faceNormal
			    );

			matrix(owner[faceID], globalCellID) += coeff;
		    }

		    // Include boundary value to the source
		    // Last item in gradCoeff refers to boundary face.
		    {
		        const label size = faceStencil.size();

		        const vector& cellGradCoeff =
		            gradCoeffs[faceID][pointI][size];

		        const tensor coeff =
		            calcCoeff
		            (
		                gammaMagSf,
		        	quadPointW,
		        	cellGradCoeff,
		        	faceNormal
			     );
			{
			    Info<<"Modified for MMS case"<<__FILE__<<__LINE__<<endl;
			    const label start = mesh.boundaryMesh()[patchI].start();
			    const point qp = faceQuadPoints[faceI+start][pointI];
			    const vector qpDisp = vector
				(
				    Foam::exp(Foam::sqr(qp.x()))*Foam::sin(qp.y()),
				    Foam::log(3+qp.y())*Foam::cos(qp.x()) + Foam::sin(qp.y()),
				    0.0
				);

			    source[owner[faceID]] -= coeff & qpDisp;

			}
                        //source[owner[faceID]] -=
		        //    coeff & D.boundaryField()[patchI][faceI];
		    }
		}
 	    }
	}
	else
	{
	    // Currently only displacement and traction boundary are implemented
	    NotImplemented;
	}
    }
}


void Foam::hofvm::addTractionBoundaries
(
    vectorField& source,
    const fvMesh& mesh,
    const volVectorField& D
)
{
    const labelUList& owner = mesh.owner();

    forAll(mesh.boundaryMesh(), patchI)
    {
	if (isA<solidTractionFvPatchVectorField>(D.boundaryField()[patchI]))
	{
	    const solidTractionFvPatchVectorField& tracPatch =
		refCast<const solidTractionFvPatchVectorField>
		(
		    D.boundaryField()[patchI]
		);

            const vectorField pNormal(mesh.boundary()[patchI].nf());

	    const vectorField traction
		(
		    tracPatch.traction() - pNormal*tracPatch.pressure()
	        );

            const label start = mesh.boundaryMesh()[patchI].start();
            const scalarField& pMagSf = mesh.magSf().boundaryField()[patchI];

	    forAll(mesh.boundaryMesh()[patchI], faceI)
	    {
		// Get global face index
		const label faceID = faceI + start;

		// Face force
		const vector force = pMagSf[faceI] * traction[faceI];

		source[owner[faceID]] -= force;
	    }
	}
    }
}

// ************************************************************************* //

#endif // #ifdef USE_PETSC
