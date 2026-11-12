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

Description
    Processor (coupled) patch contribution to the coupled W-Theta block
    matrix of coupledTotalLagNewtonRaphsonBeam, i.e. the parallel
    counterpart of the internal-face loop in assembleMatrixCoefficients().
    For each face of a coupled patch it adds, to the 6x6 blocks:
      - W equation: CQW (Laplacian) and CQTheta (weight-interpolated)
        coefficients and the explicitQ source
      - Theta equation: CMTheta, CMTheta2, the (dr x Q) terms CMQW and
        CMQTheta, and the explicitM and explicitMQ sources; the axial
        (local x) component of explicitMQ is removed via Lambda/refLambda
    The across-processor coefficients go to the interface (coupleUpper,
    coupleLower), the owner-cell part to the diagonal d and source.

    STATUS: REFERENCE ONLY - NOT COMPILED.
    This is foam-extend code (fvBlockMatrix<vector6>, coupleUpper(),
    coupleLower(), asSquare(), tensor6Field, unallocLabelList) and does not
    compile against OpenFOAM.com. It was #include-d from the matrix
    assembly in the foam-extend version (commit 8c38453,
    coupledWThetaEqn_TLNR.H) and disabled in 1e34034 (2023-09-01) when the
    solver was ported. It is kept as the derivation for parallel support.

    Where it belongs:
        coupledTotalLagNewtonRaphsonBeamBoundaryConditions.C,
        assembleBoundaryConditions(), in the "if (patch.coupled())" branch
        that currently calls notImplemented("Not ported for parallel yet!")
    The variables patchI, d and source are expected from that scope.

    Before re-enabling, port it to the current assembly:
      - d is Field<scalarSquareMatrix> (unchanged indexing) and source is
        Field<scalarRectangularMatrix>: use source[c](i, 0)
      - there is no storage for processor-interface coefficients (cu, cl):
        d, l, u only cover internal faces and BlockEigenSolverOF is a
        serial solver on local owner/neighbour addressing, so the linear
        system needs global cell numbering and a parallel solver
      - assembleMatrixCoefficients() no longer removes the axial component
        of explicitMQ; reconcile this treatment with the internal faces
      - the final sign change follows the foam-extend convention for
        interface coefficients; check it against the new solver

\*---------------------------------------------------------------------------*/

{
    // Pout << "Update coupling coeffs for patch: " << patch.name() << endl;

    tensor6Field& cu = WThetaEqn.coupleUpper()[patchI].asSquare();
    cu = tensor6::zero;

    tensor6Field& cl = WThetaEqn.coupleLower()[patchI].asSquare();
    cl = tensor6::zero;

    scalarField deltaf = 1.0/mesh().deltaCoeffs().boundaryField()[patchI];
    const scalarField& wf = mesh().weights().boundaryField()[patchI];

    const tensorField& CQWI = CQW_.boundaryField()[patchI];
    const tensorField& CQThetaI = CQTheta_.boundaryField()[patchI];
    const tensorField& CMThetaI = CMTheta_.boundaryField()[patchI];
    const tensorField& CMTheta2I = CMTheta2_.boundaryField()[patchI];
    const tensorField& CMQWI = CMQW_.boundaryField()[patchI];
    const tensorField& CMQThetaI = CMQTheta_.boundaryField()[patchI];

    const vectorField& explicitQI = explicitQ_.boundaryField()[patchI];
    const vectorField& explicitMI = explicitM_.boundaryField()[patchI];
    const vectorField& explicitMQI = explicitMQ_.boundaryField()[patchI];

    const tensorField& LambdaI = Lambda_.internalField();
    const tensorField& refLambdaI = refLambda_.internalField();
    const unallocLabelList& fc = mesh().boundary()[patchI].faceCells();

    forAll(cu, faceI)
    {
        ////// W equation

        // W part (Laplacian)
        // upper
        cu[faceI](0,0) += CQWI[faceI].xx()/deltaf[faceI];
        cu[faceI](0,1) += CQWI[faceI].xy()/deltaf[faceI];
        cu[faceI](0,2) += CQWI[faceI].xz()/deltaf[faceI];

        cu[faceI](1,0) += CQWI[faceI].yx()/deltaf[faceI];
        cu[faceI](1,1) += CQWI[faceI].yy()/deltaf[faceI];
        cu[faceI](1,2) += CQWI[faceI].yz()/deltaf[faceI];

        cu[faceI](2,0) += CQWI[faceI].zx()/deltaf[faceI];
        cu[faceI](2,1) += CQWI[faceI].zy()/deltaf[faceI];
        cu[faceI](2,2) += CQWI[faceI].zz()/deltaf[faceI];

        // lower
        cl[faceI](0,0) += -CQWI[faceI].xx()/deltaf[faceI];
        cl[faceI](0,1) += -CQWI[faceI].xy()/deltaf[faceI];
        cl[faceI](0,2) += -CQWI[faceI].xz()/deltaf[faceI];

        cl[faceI](1,0) += -CQWI[faceI].yx()/deltaf[faceI];
        cl[faceI](1,1) += -CQWI[faceI].yy()/deltaf[faceI];
        cl[faceI](1,2) += -CQWI[faceI].yz()/deltaf[faceI];

        cl[faceI](2,0) += -CQWI[faceI].zx()/deltaf[faceI];
        cl[faceI](2,1) += -CQWI[faceI].zy()/deltaf[faceI];
        cl[faceI](2,2) += -CQWI[faceI].zz()/deltaf[faceI];

        // diag
        d[fc[faceI]](0,0) += -CQWI[faceI].xx()/deltaf[faceI];
        d[fc[faceI]](0,1) += -CQWI[faceI].xy()/deltaf[faceI];
        d[fc[faceI]](0,2) += -CQWI[faceI].xz()/deltaf[faceI];

        d[fc[faceI]](1,0) += -CQWI[faceI].yx()/deltaf[faceI];
        d[fc[faceI]](1,1) += -CQWI[faceI].yy()/deltaf[faceI];
        d[fc[faceI]](1,2) += -CQWI[faceI].yz()/deltaf[faceI];

        d[fc[faceI]](2,0) += -CQWI[faceI].zx()/deltaf[faceI];
        d[fc[faceI]](2,1) += -CQWI[faceI].zy()/deltaf[faceI];
        d[fc[faceI]](2,2) += -CQWI[faceI].zz()/deltaf[faceI];


        // Theta part

        // upper
        cu[faceI](0,3) += (1-wf[faceI])*CQThetaI[faceI].xx();
        cu[faceI](0,4) += (1-wf[faceI])*CQThetaI[faceI].xy();
        cu[faceI](0,5) += (1-wf[faceI])*CQThetaI[faceI].xz();

        cu[faceI](1,3) += (1-wf[faceI])*CQThetaI[faceI].yx();
        cu[faceI](1,4) += (1-wf[faceI])*CQThetaI[faceI].yy();
        cu[faceI](1,5) += (1-wf[faceI])*CQThetaI[faceI].yz();

        cu[faceI](2,3) += (1-wf[faceI])*CQThetaI[faceI].zx();
        cu[faceI](2,4) += (1-wf[faceI])*CQThetaI[faceI].zy();
        cu[faceI](2,5) += (1-wf[faceI])*CQThetaI[faceI].zz();

        // lower
        cl[faceI](0,3) += wf[faceI]*CQThetaI[faceI].xx();
        cl[faceI](0,4) += wf[faceI]*CQThetaI[faceI].xy();
        cl[faceI](0,5) += wf[faceI]*CQThetaI[faceI].xz();

        cl[faceI](1,3) += wf[faceI]*CQThetaI[faceI].yx();
        cl[faceI](1,4) += wf[faceI]*CQThetaI[faceI].yy();
        cl[faceI](1,5) += wf[faceI]*CQThetaI[faceI].yz();

        cl[faceI](2,3) += wf[faceI]*CQThetaI[faceI].zx();
        cl[faceI](2,4) += wf[faceI]*CQThetaI[faceI].zy();
        cl[faceI](2,5) += wf[faceI]*CQThetaI[faceI].zz();

        // diag
        d[fc[faceI]](0,3) += wf[faceI]*CQThetaI[faceI].xx();
        d[fc[faceI]](0,4) += wf[faceI]*CQThetaI[faceI].xy();
        d[fc[faceI]](0,5) += wf[faceI]*CQThetaI[faceI].xz();

        d[fc[faceI]](1,3) += wf[faceI]*CQThetaI[faceI].yx();
        d[fc[faceI]](1,4) += wf[faceI]*CQThetaI[faceI].yy();
        d[fc[faceI]](1,5) += wf[faceI]*CQThetaI[faceI].yz();

        d[fc[faceI]](2,3) += wf[faceI]*CQThetaI[faceI].zx();
        d[fc[faceI]](2,4) += wf[faceI]*CQThetaI[faceI].zy();
        d[fc[faceI]](2,5) += wf[faceI]*CQThetaI[faceI].zz();

        // source
        source[fc[faceI]](0) -= explicitQI[faceI].x();
        source[fc[faceI]](1) -= explicitQI[faceI].y();
        source[fc[faceI]](2) -= explicitQI[faceI].z();


        ////// Theta equation

        // Laplacian part

        // upper
        cu[faceI](3,3) += CMThetaI[faceI].xx()/deltaf[faceI];
        cu[faceI](3,4) += CMThetaI[faceI].xy()/deltaf[faceI];
        cu[faceI](3,5) += CMThetaI[faceI].xz()/deltaf[faceI];

        cu[faceI](4,3) += CMThetaI[faceI].yx()/deltaf[faceI];
        cu[faceI](4,4) += CMThetaI[faceI].yy()/deltaf[faceI];
        cu[faceI](4,5) += CMThetaI[faceI].yz()/deltaf[faceI];

        cu[faceI](5,3) += CMThetaI[faceI].zx()/deltaf[faceI];
        cu[faceI](5,4) += CMThetaI[faceI].zy()/deltaf[faceI];
        cu[faceI](5,5) += CMThetaI[faceI].zz()/deltaf[faceI];

        // lower
        cl[faceI](3,3) += -CMThetaI[faceI].xx()/deltaf[faceI];
        cl[faceI](3,4) += -CMThetaI[faceI].xy()/deltaf[faceI];
        cl[faceI](3,5) += -CMThetaI[faceI].xz()/deltaf[faceI];

        cl[faceI](4,3) += -CMThetaI[faceI].yx()/deltaf[faceI];
        cl[faceI](4,4) += -CMThetaI[faceI].yy()/deltaf[faceI];
        cl[faceI](4,5) += -CMThetaI[faceI].yz()/deltaf[faceI];

        cl[faceI](5,3) += -CMThetaI[faceI].zx()/deltaf[faceI];
        cl[faceI](5,4) += -CMThetaI[faceI].zy()/deltaf[faceI];
        cl[faceI](5,5) += -CMThetaI[faceI].zz()/deltaf[faceI];

        // diag
        d[fc[faceI]](3,3) += -CMThetaI[faceI].xx()/deltaf[faceI];
        d[fc[faceI]](3,4) += -CMThetaI[faceI].xy()/deltaf[faceI];
        d[fc[faceI]](3,5) += -CMThetaI[faceI].xz()/deltaf[faceI];

        d[fc[faceI]](4,3) += -CMThetaI[faceI].yx()/deltaf[faceI];
        d[fc[faceI]](4,4) += -CMThetaI[faceI].yy()/deltaf[faceI];
        d[fc[faceI]](4,5) += -CMThetaI[faceI].yz()/deltaf[faceI];

        d[fc[faceI]](5,3) += -CMThetaI[faceI].zx()/deltaf[faceI];
        d[fc[faceI]](5,4) += -CMThetaI[faceI].zy()/deltaf[faceI];
        d[fc[faceI]](5,5) += -CMThetaI[faceI].zz()/deltaf[faceI];

        // Theta part

        // upper
        cu[faceI](3,3) += (1-wf[faceI])*CMTheta2I[faceI].xx();
        cu[faceI](3,4) += (1-wf[faceI])*CMTheta2I[faceI].xy();
        cu[faceI](3,5) += (1-wf[faceI])*CMTheta2I[faceI].xz();

        cu[faceI](4,3) += (1-wf[faceI])*CMTheta2I[faceI].yx();
        cu[faceI](4,4) += (1-wf[faceI])*CMTheta2I[faceI].yy();
        cu[faceI](4,5) += (1-wf[faceI])*CMTheta2I[faceI].yz();

        cu[faceI](5,3) += (1-wf[faceI])*CMTheta2I[faceI].zx();
        cu[faceI](5,4) += (1-wf[faceI])*CMTheta2I[faceI].zy();
        cu[faceI](5,5) += (1-wf[faceI])*CMTheta2I[faceI].zz();

        // lower
        cl[faceI](3,3) += wf[faceI]*CMTheta2I[faceI].xx();
        cl[faceI](3,4) += wf[faceI]*CMTheta2I[faceI].xy();
        cl[faceI](3,5) += wf[faceI]*CMTheta2I[faceI].xz();

        cl[faceI](4,3) += wf[faceI]*CMTheta2I[faceI].yx();
        cl[faceI](4,4) += wf[faceI]*CMTheta2I[faceI].yy();
        cl[faceI](4,5) += wf[faceI]*CMTheta2I[faceI].yz();

        cl[faceI](5,3) += wf[faceI]*CMTheta2I[faceI].zx();
        cl[faceI](5,4) += wf[faceI]*CMTheta2I[faceI].zy();
        cl[faceI](5,5) += wf[faceI]*CMTheta2I[faceI].zz();

        // diag
        d[fc[faceI]](3,3) += wf[faceI]*CMTheta2I[faceI].xx();
        d[fc[faceI]](3,4) += wf[faceI]*CMTheta2I[faceI].xy();
        d[fc[faceI]](3,5) += wf[faceI]*CMTheta2I[faceI].xz();

        d[fc[faceI]](4,3) += wf[faceI]*CMTheta2I[faceI].yx();
        d[fc[faceI]](4,4) += wf[faceI]*CMTheta2I[faceI].yy();
        d[fc[faceI]](4,5) += wf[faceI]*CMTheta2I[faceI].yz();

        d[fc[faceI]](5,3) += wf[faceI]*CMTheta2I[faceI].zx();
        d[fc[faceI]](5,4) += wf[faceI]*CMTheta2I[faceI].zy();
        d[fc[faceI]](5,5) += wf[faceI]*CMTheta2I[faceI].zz();


        // Explicit part
        source[fc[faceI]](3) -= explicitMI[faceI].x();
        source[fc[faceI]](4) -= explicitMI[faceI].y();
        source[fc[faceI]](5) -= explicitMI[faceI].z();


        //---- (dr x  Q) term

        // W part

        // upper
        cu[faceI](3,0) += CMQWI[faceI].xx()/deltaf[faceI];
        cu[faceI](3,1) += CMQWI[faceI].xy()/deltaf[faceI];
        cu[faceI](3,2) += CMQWI[faceI].xz()/deltaf[faceI];

        cu[faceI](4,0) += CMQWI[faceI].yx()/deltaf[faceI];
        cu[faceI](4,1) += CMQWI[faceI].yy()/deltaf[faceI];
        cu[faceI](4,2) += CMQWI[faceI].yz()/deltaf[faceI];

        cu[faceI](5,0) += CMQWI[faceI].zx()/deltaf[faceI];
        cu[faceI](5,1) += CMQWI[faceI].zy()/deltaf[faceI];
        cu[faceI](5,2) += CMQWI[faceI].zz()/deltaf[faceI];

        // lower
        cl[faceI](3,0) += -CMQWI[faceI].xx()/deltaf[faceI];
        cl[faceI](3,1) += -CMQWI[faceI].xy()/deltaf[faceI];
        cl[faceI](3,2) += -CMQWI[faceI].xz()/deltaf[faceI];

        cl[faceI](4,0) += -CMQWI[faceI].yx()/deltaf[faceI];
        cl[faceI](4,1) += -CMQWI[faceI].yy()/deltaf[faceI];
        cl[faceI](4,2) += -CMQWI[faceI].yz()/deltaf[faceI];

        cl[faceI](5,0) += -CMQWI[faceI].zx()/deltaf[faceI];
        cl[faceI](5,1) += -CMQWI[faceI].zy()/deltaf[faceI];
        cl[faceI](5,2) += -CMQWI[faceI].zz()/deltaf[faceI];

        // diag
        d[fc[faceI]](3,0) += -CMQWI[faceI].xx()/deltaf[faceI];
        d[fc[faceI]](3,1) += -CMQWI[faceI].xy()/deltaf[faceI];
        d[fc[faceI]](3,2) += -CMQWI[faceI].xz()/deltaf[faceI];

        d[fc[faceI]](4,0) += -CMQWI[faceI].yx()/deltaf[faceI];
        d[fc[faceI]](4,1) += -CMQWI[faceI].yy()/deltaf[faceI];
        d[fc[faceI]](4,2) += -CMQWI[faceI].yz()/deltaf[faceI];

        d[fc[faceI]](5,0) += -CMQWI[faceI].zx()/deltaf[faceI];
        d[fc[faceI]](5,1) += -CMQWI[faceI].zy()/deltaf[faceI];
        d[fc[faceI]](5,2) += -CMQWI[faceI].zz()/deltaf[faceI];


        // Theta part

        // upper
        cu[faceI](3,3) += (1-wf[faceI])*CMQThetaI[faceI].xx();
        cu[faceI](3,4) += (1-wf[faceI])*CMQThetaI[faceI].xy();
        cu[faceI](3,5) += (1-wf[faceI])*CMQThetaI[faceI].xz();

        cu[faceI](4,3) += (1-wf[faceI])*CMQThetaI[faceI].yx();
        cu[faceI](4,4) += (1-wf[faceI])*CMQThetaI[faceI].yy();
        cu[faceI](4,5) += (1-wf[faceI])*CMQThetaI[faceI].yz();

        cu[faceI](5,3) += (1-wf[faceI])*CMQThetaI[faceI].zx();
        cu[faceI](5,4) += (1-wf[faceI])*CMQThetaI[faceI].zy();
        cu[faceI](5,5) += (1-wf[faceI])*CMQThetaI[faceI].zz();

        // lower
        cl[faceI](3,3) += wf[faceI]*CMQThetaI[faceI].xx();
        cl[faceI](3,4) += wf[faceI]*CMQThetaI[faceI].xy();
        cl[faceI](3,5) += wf[faceI]*CMQThetaI[faceI].xz();

        cl[faceI](4,3) += wf[faceI]*CMQThetaI[faceI].yx();
        cl[faceI](4,4) += wf[faceI]*CMQThetaI[faceI].yy();
        cl[faceI](4,5) += wf[faceI]*CMQThetaI[faceI].yz();

        cl[faceI](5,3) += wf[faceI]*CMQThetaI[faceI].zx();
        cl[faceI](5,4) += wf[faceI]*CMQThetaI[faceI].zy();
        cl[faceI](5,5) += wf[faceI]*CMQThetaI[faceI].zz();

        // diag
        d[fc[faceI]](3,3) += wf[faceI]*CMQThetaI[faceI].xx();
        d[fc[faceI]](3,4) += wf[faceI]*CMQThetaI[faceI].xy();
        d[fc[faceI]](3,5) += wf[faceI]*CMQThetaI[faceI].xz();

        d[fc[faceI]](4,3) += wf[faceI]*CMQThetaI[faceI].yx();
        d[fc[faceI]](4,4) += wf[faceI]*CMQThetaI[faceI].yy();
        d[fc[faceI]](4,5) += wf[faceI]*CMQThetaI[faceI].yz();

        d[fc[faceI]](5,3) += wf[faceI]*CMQThetaI[faceI].zx();
        d[fc[faceI]](5,4) += wf[faceI]*CMQThetaI[faceI].zy();
        d[fc[faceI]](5,5) += wf[faceI]*CMQThetaI[faceI].zz();

        // Explicit part
        vector correctedOwnExplicitMQ = explicitMQI[faceI];

        correctedOwnExplicitMQ =
        (
            refLambdaI[fc[faceI]].T()
          & (LambdaI[fc[faceI]].T() & correctedOwnExplicitMQ)
        );
        correctedOwnExplicitMQ.x() = 0;
        correctedOwnExplicitMQ =
            (LambdaI[fc[faceI]] & (refLambdaI[fc[faceI]] & correctedOwnExplicitMQ));

        source[fc[faceI]](3) -= correctedOwnExplicitMQ.x();
        source[fc[faceI]](4) -= correctedOwnExplicitMQ.y();
        source[fc[faceI]](5) -= correctedOwnExplicitMQ.z();

        // Change sign
        cu[faceI] *= -1;
        cl[faceI] *= -1;
    }

    // Info <<  WThetaEqn.coupleLower()[patchI].asSquare() << endl;
    // Info <<  WThetaEqn.coupleUpper()[patchI].asSquare() << endl;

    // Pout << "cu: " << cu << endl;
}
