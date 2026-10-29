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

#include "mergePoints.H"
#include "scalarMatrices.H"
#include "denseMatrixHelperFunctions.H"
#include "BlockEigenSolverOF.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace beamModels
{


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void coupledTotalLagNewtonRaphsonBeam::assembleMatrixCoefficients
(
    Field<scalarSquareMatrix>& d,
    Field<scalarSquareMatrix>& l,
    Field<scalarSquareMatrix>& u,
    Field<scalarRectangularMatrix>& source
)
{
    const tensorField& CQWI = CQW_.internalField();
    const tensorField& CQThetaI = CQTheta_.internalField();
    const tensorField& CMThetaI = CMTheta_.internalField();
    const tensorField& CMTheta2I = CMTheta2_.internalField();
    const tensorField& CMQWI = CMQW_.internalField();
    const tensorField& CMQThetaI = CMQTheta_.internalField();
    const vectorField& explicitQI = explicitQ_.internalField();
    const vectorField& explicitMI = explicitM_.internalField();
    const vectorField& explicitMQI = explicitMQ_.internalField();


    const scalarField deltaf = 1.0/mesh().deltaCoeffs().internalField();
    const scalarField& wf = mesh().weights().internalField();
    const labelList& own = mesh().owner(); // unallocLabelList => labelList (ESI)
    const labelList& nei = mesh().neighbour();

    // Internal faces
    forAll(u, faceI)
    {
        //-W part (Laplacian)
        u[faceI](0,0) += CQWI[faceI].xx()/deltaf[faceI];
        u[faceI](0,1) += CQWI[faceI].xy()/deltaf[faceI];
        u[faceI](0,2) += CQWI[faceI].xz()/deltaf[faceI];

        u[faceI](1,0) += CQWI[faceI].yx()/deltaf[faceI];
        u[faceI](1,1) += CQWI[faceI].yy()/deltaf[faceI];
        u[faceI](1,2) += CQWI[faceI].yz()/deltaf[faceI];

        u[faceI](2,0) += CQWI[faceI].zx()/deltaf[faceI];
        u[faceI](2,1) += CQWI[faceI].zy()/deltaf[faceI];
        u[faceI](2,2) += CQWI[faceI].zz()/deltaf[faceI];

        d[own[faceI]](0,0) += -CQWI[faceI].xx()/deltaf[faceI];
        d[own[faceI]](0,1) += -CQWI[faceI].xy()/deltaf[faceI];
        d[own[faceI]](0,2) += -CQWI[faceI].xz()/deltaf[faceI];

        d[own[faceI]](1,0) += -CQWI[faceI].yx()/deltaf[faceI];
        d[own[faceI]](1,1) += -CQWI[faceI].yy()/deltaf[faceI];
        d[own[faceI]](1,2) += -CQWI[faceI].yz()/deltaf[faceI];

        d[own[faceI]](2,0) += -CQWI[faceI].zx()/deltaf[faceI];
        d[own[faceI]](2,1) += -CQWI[faceI].zy()/deltaf[faceI];
        d[own[faceI]](2,2) += -CQWI[faceI].zz()/deltaf[faceI];

        l[faceI](0,0) += CQWI[faceI].xx()/deltaf[faceI];
        l[faceI](0,1) += CQWI[faceI].xy()/deltaf[faceI];
        l[faceI](0,2) += CQWI[faceI].xz()/deltaf[faceI];

        l[faceI](1,0) += CQWI[faceI].yx()/deltaf[faceI];
        l[faceI](1,1) += CQWI[faceI].yy()/deltaf[faceI];
        l[faceI](1,2) += CQWI[faceI].yz()/deltaf[faceI];

        l[faceI](2,0) += CQWI[faceI].zx()/deltaf[faceI];
        l[faceI](2,1) += CQWI[faceI].zy()/deltaf[faceI];
        l[faceI](2,2) += CQWI[faceI].zz()/deltaf[faceI];

        d[nei[faceI]](0,0) += -CQWI[faceI].xx()/deltaf[faceI];
        d[nei[faceI]](0,1) += -CQWI[faceI].xy()/deltaf[faceI];
        d[nei[faceI]](0,2) += -CQWI[faceI].xz()/deltaf[faceI];

        d[nei[faceI]](1,0) += -CQWI[faceI].yx()/deltaf[faceI];
        d[nei[faceI]](1,1) += -CQWI[faceI].yy()/deltaf[faceI];
        d[nei[faceI]](1,2) += -CQWI[faceI].yz()/deltaf[faceI];

        d[nei[faceI]](2,0) += -CQWI[faceI].zx()/deltaf[faceI];
        d[nei[faceI]](2,1) += -CQWI[faceI].zy()/deltaf[faceI];
        d[nei[faceI]](2,2) += -CQWI[faceI].zz()/deltaf[faceI];

        //- Theta part
        u[faceI](0,3) += (1 - wf[faceI])*CQThetaI[faceI].xx();
        u[faceI](0,4) += (1 - wf[faceI])*CQThetaI[faceI].xy();
        u[faceI](0,5) += (1 - wf[faceI])*CQThetaI[faceI].xz();

        u[faceI](1,3) += (1 - wf[faceI])*CQThetaI[faceI].yx();
        u[faceI](1,4) += (1 - wf[faceI])*CQThetaI[faceI].yy();
        u[faceI](1,5) += (1 - wf[faceI])*CQThetaI[faceI].yz();

        u[faceI](2,3) += (1 - wf[faceI])*CQThetaI[faceI].zx();
        u[faceI](2,4) += (1 - wf[faceI])*CQThetaI[faceI].zy();
        u[faceI](2,5) += (1 - wf[faceI])*CQThetaI[faceI].zz();

        d[own[faceI]](0,3) += wf[faceI]*CQThetaI[faceI].xx();
        d[own[faceI]](0,4) += wf[faceI]*CQThetaI[faceI].xy();
        d[own[faceI]](0,5) += wf[faceI]*CQThetaI[faceI].xz();

        d[own[faceI]](1,3) += wf[faceI]*CQThetaI[faceI].yx();
        d[own[faceI]](1,4) += wf[faceI]*CQThetaI[faceI].yy();
        d[own[faceI]](1,5) += wf[faceI]*CQThetaI[faceI].yz();

        d[own[faceI]](2,3) += wf[faceI]*CQThetaI[faceI].zx();
        d[own[faceI]](2,4) += wf[faceI]*CQThetaI[faceI].zy();
        d[own[faceI]](2,5) += wf[faceI]*CQThetaI[faceI].zz();

        source[own[faceI]](0,0) -= explicitQI[faceI].x();
        source[own[faceI]](1,0) -= explicitQI[faceI].y();
        source[own[faceI]](2,0) -= explicitQI[faceI].z();

        l[faceI](0,3) += -wf[faceI]*CQThetaI[faceI].xx();
        l[faceI](0,4) += -wf[faceI]*CQThetaI[faceI].xy();
        l[faceI](0,5) += -wf[faceI]*CQThetaI[faceI].xz();

        l[faceI](1,3) += -wf[faceI]*CQThetaI[faceI].yx();
        l[faceI](1,4) += -wf[faceI]*CQThetaI[faceI].yy();
        l[faceI](1,5) += -wf[faceI]*CQThetaI[faceI].yz();

        l[faceI](2,3) += -wf[faceI]*CQThetaI[faceI].zx();
        l[faceI](2,4) += -wf[faceI]*CQThetaI[faceI].zy();
        l[faceI](2,5) += -wf[faceI]*CQThetaI[faceI].zz();

        d[nei[faceI]](0,3) += -(1 - wf[faceI])*CQThetaI[faceI].xx();
        d[nei[faceI]](0,4) += -(1 - wf[faceI])*CQThetaI[faceI].xy();
        d[nei[faceI]](0,5) += -(1 - wf[faceI])*CQThetaI[faceI].xz();

        d[nei[faceI]](1,3) += -(1 - wf[faceI])*CQThetaI[faceI].yx();
        d[nei[faceI]](1,4) += -(1 - wf[faceI])*CQThetaI[faceI].yy();
        d[nei[faceI]](1,5) += -(1 - wf[faceI])*CQThetaI[faceI].yz();

        d[nei[faceI]](2,3) += -(1 - wf[faceI])*CQThetaI[faceI].zx();
        d[nei[faceI]](2,4) += -(1 - wf[faceI])*CQThetaI[faceI].zy();
        d[nei[faceI]](2,5) += -(1 - wf[faceI])*CQThetaI[faceI].zz();

        source[nei[faceI]](0,0) -= -explicitQI[faceI].x();
        source[nei[faceI]](1,0) -= -explicitQI[faceI].y();
        source[nei[faceI]](2,0) -= -explicitQI[faceI].z();


        ////// Theta equation

        //- Laplacian part

        u[faceI](3,3) += CMThetaI[faceI].xx()/deltaf[faceI];
        u[faceI](3,4) += CMThetaI[faceI].xy()/deltaf[faceI];
        u[faceI](3,5) += CMThetaI[faceI].xz()/deltaf[faceI];

        u[faceI](4,3) += CMThetaI[faceI].yx()/deltaf[faceI];
        u[faceI](4,4) += CMThetaI[faceI].yy()/deltaf[faceI];
        u[faceI](4,5) += CMThetaI[faceI].yz()/deltaf[faceI];

        u[faceI](5,3) += CMThetaI[faceI].zx()/deltaf[faceI];
        u[faceI](5,4) += CMThetaI[faceI].zy()/deltaf[faceI];
        u[faceI](5,5) += CMThetaI[faceI].zz()/deltaf[faceI];

        d[own[faceI]](3,3) += -CMThetaI[faceI].xx()/deltaf[faceI];
        d[own[faceI]](3,4) += -CMThetaI[faceI].xy()/deltaf[faceI];
        d[own[faceI]](3,5) += -CMThetaI[faceI].xz()/deltaf[faceI];

        d[own[faceI]](4,3) += -CMThetaI[faceI].yx()/deltaf[faceI];
        d[own[faceI]](4,4) += -CMThetaI[faceI].yy()/deltaf[faceI];
        d[own[faceI]](4,5) += -CMThetaI[faceI].yz()/deltaf[faceI];

        d[own[faceI]](5,3) += -CMThetaI[faceI].zx()/deltaf[faceI];
        d[own[faceI]](5,4) += -CMThetaI[faceI].zy()/deltaf[faceI];
        d[own[faceI]](5,5) += -CMThetaI[faceI].zz()/deltaf[faceI];

        l[faceI](3,3) += CMThetaI[faceI].xx()/deltaf[faceI];
        l[faceI](3,4) += CMThetaI[faceI].xy()/deltaf[faceI];
        l[faceI](3,5) += CMThetaI[faceI].xz()/deltaf[faceI];

        l[faceI](4,3) += CMThetaI[faceI].yx()/deltaf[faceI];
        l[faceI](4,4) += CMThetaI[faceI].yy()/deltaf[faceI];
        l[faceI](4,5) += CMThetaI[faceI].yz()/deltaf[faceI];

        l[faceI](5,3) += CMThetaI[faceI].zx()/deltaf[faceI];
        l[faceI](5,4) += CMThetaI[faceI].zy()/deltaf[faceI];
        l[faceI](5,5) += CMThetaI[faceI].zz()/deltaf[faceI];

        d[nei[faceI]](3,3) += -CMThetaI[faceI].xx()/deltaf[faceI];
        d[nei[faceI]](3,4) += -CMThetaI[faceI].xy()/deltaf[faceI];
        d[nei[faceI]](3,5) += -CMThetaI[faceI].xz()/deltaf[faceI];

        d[nei[faceI]](4,3) += -CMThetaI[faceI].yx()/deltaf[faceI];
        d[nei[faceI]](4,4) += -CMThetaI[faceI].yy()/deltaf[faceI];
        d[nei[faceI]](4,5) += -CMThetaI[faceI].yz()/deltaf[faceI];

        d[nei[faceI]](5,3) += -CMThetaI[faceI].zx()/deltaf[faceI];
        d[nei[faceI]](5,4) += -CMThetaI[faceI].zy()/deltaf[faceI];
        d[nei[faceI]](5,5) += -CMThetaI[faceI].zz()/deltaf[faceI];


        //- Theta part

        u[faceI](3,3) += (1 - wf[faceI])*CMTheta2I[faceI].xx();
        u[faceI](3,4) += (1 - wf[faceI])*CMTheta2I[faceI].xy();
        u[faceI](3,5) += (1 - wf[faceI])*CMTheta2I[faceI].xz();

        u[faceI](4,3) += (1 - wf[faceI])*CMTheta2I[faceI].yx();
        u[faceI](4,4) += (1 - wf[faceI])*CMTheta2I[faceI].yy();
        u[faceI](4,5) += (1 - wf[faceI])*CMTheta2I[faceI].yz();

        u[faceI](5,3) += (1 - wf[faceI])*CMTheta2I[faceI].zx();
        u[faceI](5,4) += (1 - wf[faceI])*CMTheta2I[faceI].zy();
        u[faceI](5,5) += (1 - wf[faceI])*CMTheta2I[faceI].zz();

        d[own[faceI]](3,3) += wf[faceI]*CMTheta2I[faceI].xx();
        d[own[faceI]](3,4) += wf[faceI]*CMTheta2I[faceI].xy();
        d[own[faceI]](3,5) += wf[faceI]*CMTheta2I[faceI].xz();

        d[own[faceI]](4,3) += wf[faceI]*CMTheta2I[faceI].yx();
        d[own[faceI]](4,4) += wf[faceI]*CMTheta2I[faceI].yy();
        d[own[faceI]](4,5) += wf[faceI]*CMTheta2I[faceI].yz();

        d[own[faceI]](5,3) += wf[faceI]*CMTheta2I[faceI].zx();
        d[own[faceI]](5,4) += wf[faceI]*CMTheta2I[faceI].zy();
        d[own[faceI]](5,5) += wf[faceI]*CMTheta2I[faceI].zz();


        l[faceI](3,3) += -wf[faceI]*CMTheta2I[faceI].xx();
        l[faceI](3,4) += -wf[faceI]*CMTheta2I[faceI].xy();
        l[faceI](3,5) += -wf[faceI]*CMTheta2I[faceI].xz();

        l[faceI](4,3) += -wf[faceI]*CMTheta2I[faceI].yx();
        l[faceI](4,4) += -wf[faceI]*CMTheta2I[faceI].yy();
        l[faceI](4,5) += -wf[faceI]*CMTheta2I[faceI].yz();

        l[faceI](5,3) += -wf[faceI]*CMTheta2I[faceI].zx();
        l[faceI](5,4) += -wf[faceI]*CMTheta2I[faceI].zy();
        l[faceI](5,5) += -wf[faceI]*CMTheta2I[faceI].zz();

        d[nei[faceI]](3,3) += -(1 - wf[faceI])*CMTheta2I[faceI].xx();
        d[nei[faceI]](3,4) += -(1 - wf[faceI])*CMTheta2I[faceI].xy();
        d[nei[faceI]](3,5) += -(1 - wf[faceI])*CMTheta2I[faceI].xz();

        d[nei[faceI]](4,3) += -(1 - wf[faceI])*CMTheta2I[faceI].yx();
        d[nei[faceI]](4,4) += -(1 - wf[faceI])*CMTheta2I[faceI].yy();
        d[nei[faceI]](4,5) += -(1 - wf[faceI])*CMTheta2I[faceI].yz();

        d[nei[faceI]](5,3) += -(1 - wf[faceI])*CMTheta2I[faceI].zx();
        d[nei[faceI]](5,4) += -(1 - wf[faceI])*CMTheta2I[faceI].zy();
        d[nei[faceI]](5,5) += -(1 - wf[faceI])*CMTheta2I[faceI].zz();

        // Explicit part
        source[own[faceI]](3,0) -= explicitMI[faceI].x();
        source[own[faceI]](4,0) -= explicitMI[faceI].y();
        source[own[faceI]](5,0) -= explicitMI[faceI].z();

        source[nei[faceI]](3,0) -= -explicitMI[faceI].x();
        source[nei[faceI]](4,0) -= -explicitMI[faceI].y();
        source[nei[faceI]](5,0) -= -explicitMI[faceI].z();



        //---- (dr x  Q) term

        // W part
        u[faceI](3,0) += CMQWI[faceI].xx()/deltaf[faceI];
        u[faceI](3,1) += CMQWI[faceI].xy()/deltaf[faceI];
        u[faceI](3,2) += CMQWI[faceI].xz()/deltaf[faceI];

        u[faceI](4,0) += CMQWI[faceI].yx()/deltaf[faceI];
        u[faceI](4,1) += CMQWI[faceI].yy()/deltaf[faceI];
        u[faceI](4,2) += CMQWI[faceI].yz()/deltaf[faceI];

        u[faceI](5,0) += CMQWI[faceI].zx()/deltaf[faceI];
        u[faceI](5,1) += CMQWI[faceI].zy()/deltaf[faceI];
        u[faceI](5,2) += CMQWI[faceI].zz()/deltaf[faceI];

        d[own[faceI]](3,0) += -CMQWI[faceI].xx()/deltaf[faceI];
        d[own[faceI]](3,1) += -CMQWI[faceI].xy()/deltaf[faceI];
        d[own[faceI]](3,2) += -CMQWI[faceI].xz()/deltaf[faceI];

        d[own[faceI]](4,0) += -CMQWI[faceI].yx()/deltaf[faceI];
        d[own[faceI]](4,1) += -CMQWI[faceI].yy()/deltaf[faceI];
        d[own[faceI]](4,2) += -CMQWI[faceI].yz()/deltaf[faceI];

        d[own[faceI]](5,0) += -CMQWI[faceI].zx()/deltaf[faceI];
        d[own[faceI]](5,1) += -CMQWI[faceI].zy()/deltaf[faceI];
        d[own[faceI]](5,2) += -CMQWI[faceI].zz()/deltaf[faceI];

        l[faceI](3,0) += -CMQWI[faceI].xx()/deltaf[faceI];
        l[faceI](3,1) += -CMQWI[faceI].xy()/deltaf[faceI];
        l[faceI](3,2) += -CMQWI[faceI].xz()/deltaf[faceI];

        l[faceI](4,0) += -CMQWI[faceI].yx()/deltaf[faceI];
        l[faceI](4,1) += -CMQWI[faceI].yy()/deltaf[faceI];
        l[faceI](4,2) += -CMQWI[faceI].yz()/deltaf[faceI];

        l[faceI](5,0) += -CMQWI[faceI].zx()/deltaf[faceI];
        l[faceI](5,1) += -CMQWI[faceI].zy()/deltaf[faceI];
        l[faceI](5,2) += -CMQWI[faceI].zz()/deltaf[faceI];

        d[nei[faceI]](3,0) += CMQWI[faceI].xx()/deltaf[faceI];
        d[nei[faceI]](3,1) += CMQWI[faceI].xy()/deltaf[faceI];
        d[nei[faceI]](3,2) += CMQWI[faceI].xz()/deltaf[faceI];

        d[nei[faceI]](4,0) += CMQWI[faceI].yx()/deltaf[faceI];
        d[nei[faceI]](4,1) += CMQWI[faceI].yy()/deltaf[faceI];
        d[nei[faceI]](4,2) += CMQWI[faceI].yz()/deltaf[faceI];

        d[nei[faceI]](5,0) += CMQWI[faceI].zx()/deltaf[faceI];
        d[nei[faceI]](5,1) += CMQWI[faceI].zy()/deltaf[faceI];
        d[nei[faceI]](5,2) += CMQWI[faceI].zz()/deltaf[faceI];

        // Theta part
        u[faceI](3,3) += (1 - wf[faceI])*CMQThetaI[faceI].xx();
        u[faceI](3,4) += (1 - wf[faceI])*CMQThetaI[faceI].xy();
        u[faceI](3,5) += (1 - wf[faceI])*CMQThetaI[faceI].xz();

        u[faceI](4,3) += (1 - wf[faceI])*CMQThetaI[faceI].yx();
        u[faceI](4,4) += (1 - wf[faceI])*CMQThetaI[faceI].yy();
        u[faceI](4,5) += (1 - wf[faceI])*CMQThetaI[faceI].yz();

        u[faceI](5,3) += (1 - wf[faceI])*CMQThetaI[faceI].zx();
        u[faceI](5,4) += (1 - wf[faceI])*CMQThetaI[faceI].zy();
        u[faceI](5,5) += (1 - wf[faceI])*CMQThetaI[faceI].zz();

        d[own[faceI]](3,3) += wf[faceI]*CMQThetaI[faceI].xx();
        d[own[faceI]](3,4) += wf[faceI]*CMQThetaI[faceI].xy();
        d[own[faceI]](3,5) += wf[faceI]*CMQThetaI[faceI].xz();

        d[own[faceI]](4,3) += wf[faceI]*CMQThetaI[faceI].yx();
        d[own[faceI]](4,4) += wf[faceI]*CMQThetaI[faceI].yy();
        d[own[faceI]](4,5) += wf[faceI]*CMQThetaI[faceI].yz();

        d[own[faceI]](5,3) += wf[faceI]*CMQThetaI[faceI].zx();
        d[own[faceI]](5,4) += wf[faceI]*CMQThetaI[faceI].zy();
        d[own[faceI]](5,5) += wf[faceI]*CMQThetaI[faceI].zz();


        // Shouldn't the lower and diag[nei] coefficients have a negative sign?
        l[faceI](3,3) += wf[faceI]*CMQThetaI[faceI].xx();
        l[faceI](3,4) += wf[faceI]*CMQThetaI[faceI].xy();
        l[faceI](3,5) += wf[faceI]*CMQThetaI[faceI].xz();

        l[faceI](4,3) += wf[faceI]*CMQThetaI[faceI].yx();
        l[faceI](4,4) += wf[faceI]*CMQThetaI[faceI].yy();
        l[faceI](4,5) += wf[faceI]*CMQThetaI[faceI].yz();

        l[faceI](5,3) += wf[faceI]*CMQThetaI[faceI].zx();
        l[faceI](5,4) += wf[faceI]*CMQThetaI[faceI].zy();
        l[faceI](5,5) += wf[faceI]*CMQThetaI[faceI].zz();

        d[nei[faceI]](3,3) += (1 - wf[faceI])*CMQThetaI[faceI].xx();
        d[nei[faceI]](3,4) += (1 - wf[faceI])*CMQThetaI[faceI].xy();
        d[nei[faceI]](3,5) += (1 - wf[faceI])*CMQThetaI[faceI].xz();

        d[nei[faceI]](4,3) += (1 - wf[faceI])*CMQThetaI[faceI].yx();
        d[nei[faceI]](4,4) += (1 - wf[faceI])*CMQThetaI[faceI].yy();
        d[nei[faceI]](4,5) += (1 - wf[faceI])*CMQThetaI[faceI].yz();

        d[nei[faceI]](5,3) += (1 - wf[faceI])*CMQThetaI[faceI].zx();
        d[nei[faceI]](5,4) += (1 - wf[faceI])*CMQThetaI[faceI].zy();
        d[nei[faceI]](5,5) += (1 - wf[faceI])*CMQThetaI[faceI].zz();

        // Explicit part
        vector correctedOwnExplicitMQ = explicitMQI[faceI];

        source[own[faceI]](3,0) -= correctedOwnExplicitMQ.x();
        source[own[faceI]](4,0) -= correctedOwnExplicitMQ.y();
        source[own[faceI]](5,0) -= correctedOwnExplicitMQ.z();

        vector correctedNeiExplicitMQ = explicitMQI[faceI];

        // Similarly here, the explicit term is not in the negative format? // understood this
        source[nei[faceI]](3,0) -= correctedNeiExplicitMQ.x();
        source[nei[faceI]](4,0) -= correctedNeiExplicitMQ.y();
        source[nei[faceI]](5,0) -= correctedNeiExplicitMQ.z();
    }

    // Include the boundary contributions to l,d,u and source fields
    assembleBoundaryConditions(d, l, u, source);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace beamModels

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
