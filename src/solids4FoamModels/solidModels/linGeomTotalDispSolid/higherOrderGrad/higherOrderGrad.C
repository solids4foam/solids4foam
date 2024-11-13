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

#include "higherOrderGrad.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "emptyPolyPatch.H"
#include <Eigen/Dense>


namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(higherOrderGrad, 0);


// * * * * * * * * * * *  Private Member Functions * * * * * * * * * * * * * //


void higherOrderGrad::makeStencils() const
{
    Info<< "higherOrderGrad::makeStencils()" << endl;

    if (stencilsPtr_)
    {
        FatalErrorInFunction
            << "Pointer already set!" << abort(FatalError);
    }

    const fvMesh& mesh = mesh_;
    const labelListList& cellCells = mesh.cellCells();

    stencilsPtr_.set(new List<DynamicList<label>>(mesh.nCells()));
    List<DynamicList<label>>& stencils = *stencilsPtr_;

    forAll(stencils, cellI)
    {
       DynamicList<label>& curStencil = stencils[cellI];
       curStencil.setCapacity(maxStencilSize_);
       const labelList& curCellCells = cellCells[cellI];

       labelHashSet stencilCells;
       labelHashSet prevLayer;

       // Add first layer of cells
       forAll(curCellCells, cI)
       {
           stencilCells.insert(curCellCells[cI]);
           prevLayer.insert(curCellCells[cI]);
       }

       // Remaining layers of cells
       for (int layerI = 1; layerI < nLayers_; layerI++)
       {
           labelList prevLayerCells(prevLayer.toc());
           labelHashSet curLayer;

           // Loop over previous layer and add one level of
           // layer neighbours
           for (const label cI : prevLayerCells)
           {
               const labelList& cellINei= mesh.cellCells()[cI];
               forAll (cellINei, nei)
               {
                   if (!stencilCells.found(cellINei[nei]))
                   {
                       curLayer.insert(cellINei[nei]);
                   }
               }
           }

           // Now we have curent layer which we need to add to the stencil
           // and the current layer will now be the previous layer for next
           // loop
           prevLayer.clear();
           prevLayer = curLayer;
           stencilCells.merge(curLayer);
       }

       if (Pstream::parRun())
       {
           notImplemented("not implemented for parallel run");
       }

       // Neighbours of first layer will store cellI in stencil so we
       // need to remove it
       stencilCells.erase(cellI);

       curStencil.append(stencilCells.toc());
    }

    forAll(stencils, cellI)
    {
       // Note: stencils are sorted but I do not see problem in that
       stencils[cellI].shrink();
    }

    Info<< "higherOrderGrad::makeStencils(): end" << endl;
}


const List<DynamicList<label>> higherOrderGrad::stencils() const
{
    if (!stencilsPtr_)
    {
        makeStencils();
    }

    return *stencilsPtr_;
}


void higherOrderGrad::calcCoeffs() const
{
    Info<< "higherOrderGrad::calcCoeffs()" << endl;

    if (interpCoeffsPtr_ || interpGradCoeffsPtr_ || cellBoundaryFacesPtr_)
    {
        FatalErrorInFunction
            << "Pointers already set!" << abort(FatalError);
    }

    const fvMesh& mesh = mesh_;

    interpCoeffsPtr_.set(new List<DynamicList<scalar>>(mesh.nCells()));
    List<DynamicList<scalar>>& interpCoeffs = *interpCoeffsPtr_;

    interpGradCoeffsPtr_.set(new List<DynamicList<vector>>(mesh.nCells()));
    List<DynamicList<vector>>& interpGradCoeffs = *interpGradCoeffsPtr_;

    cellBoundaryFacesPtr_.set(new List<DynamicList<label>>(mesh.nCells()));
    List<DynamicList<label>>& cellBoundaryFaces = *cellBoundaryFacesPtr_;

    // Refernces for brevity and efficiency
    const vectorField& CI = mesh.C();
    const vectorField& CfI = mesh.Cf();

    // Number of terms in Taylor expression
    // 1 for zero order, 4 for 1 order, 10 for second order
    const label Np = factorial(N_ + 3)/(factorial(N_)*factorial(3));

    if (calcConditionNumber_)
    {
        // Cells condition numbers
        condNumberPtr_.set
        (
            new volScalarField
            (
               IOobject
               (
                   "condNumber",
                   mesh.time().timeName(),
                   mesh,
                   IOobject::NO_READ,
                   IOobject::AUTO_WRITE
               ),
               mesh,
               dimensionedScalar("0", dimless, Zero)
            )
        );
    }

    List<DynamicList<scalar>> c(mesh.nCells());
    List<DynamicList<scalar>> cx(mesh.nCells());
    List<DynamicList<scalar>> cy(mesh.nCells());
    List<DynamicList<scalar>> cz(mesh.nCells());

    const polyBoundaryMesh& boundaryMesh = mesh.boundaryMesh();

    forAll(boundaryMesh, patchI)
    {
        if (boundaryMesh[patchI].type() != emptyPolyPatch::typeName)
        {
            const labelUList& faceCells = boundaryMesh[patchI].faceCells();

            forAll(faceCells, faceI)
            {
                const label cellID = faceCells[faceI];
                const label gI = faceI + boundaryMesh[patchI].start();
                cellBoundaryFaces[cellID].append(gI);
            }
        }
    }

    forAll(cellBoundaryFaces, cellI)
    {
        cellBoundaryFaces[cellI].shrink();
    }

    const List<DynamicList<label>> stencils = this->stencils();

    forAll(stencils, cellI)
    {
        const DynamicList<label>& curStencil = stencils[cellI];

        // Find max distance in this stencil
        scalar maxDist = 0.0;
        forAll(curStencil, cI)
        {
            const label neiCellID = curStencil[cI];
            const scalar d = mag(CI[neiCellID] - CI[cellI]);
            if (d > maxDist)
            {
                maxDist = d;
            }
        }

        // Loop over neighbours and construct matrix Q
        const label Nn = curStencil.size() + cellBoundaryFaces[cellI].size();

        // For now I will use matrix format from Eigen/Dense library
        Eigen::MatrixXd Q = Eigen::MatrixXd::Zero(Np, Nn);

        // Check to avoid Eigen error
        if (Nn < Np)
        {
            FatalErrorInFunction
                << "Interpolation stencil needs to be bigger than the "
               "number of elements in Taylor order!"
                << exit(FatalError);
        }

        // Loop over cells in stencil, each cell have its corresponding
        // row in Q matrix
        for (label cI = 0; cI < Nn; cI++)
        {
            vector neiC = vector::zero;

            if (cI < curStencil.size())
            {
                const label neiCellID = curStencil[cI];
                neiC = CI[neiCellID];
            }
            else
            {
                // For boundary cells we need to add boundary face as
                // neigbour
                const label i = cI - curStencil.size();
                const label globalFaceID = cellBoundaryFaces[cellI][i];
                neiC = CfI[globalFaceID];
            }

            const vector& C = CI[cellI];

            // Add linear interpolation part N = 1
            if (N_ > 0)
            {
                Q(0, cI) = 1;
                Q(1, cI) = neiC.x() - C.x();
                Q(2, cI) = neiC.y() - C.y();
                Q(3, cI) = neiC.z() - C.z();
            }

            // Add quadratic interpolation part N = 2
            if (N_ > 1)
            {
                Q(4, cI) = (1.0/2.0)*pow(neiC.x() - C.x(), 2);
                Q(5, cI) = (1.0/2.0)*pow(neiC.y() - C.y(), 2);
                Q(6, cI) = (1.0/2.0)*pow(neiC.z() - C.z(), 2);
                Q(7, cI) = (neiC.x() - C.x())*(neiC.y() - C.y());
                Q(8, cI) = (neiC.x() - C.x())*(neiC.z() - C.z());
                Q(9, cI) = (neiC.y() - C.y())*(neiC.z() - C.z());
            }

            // Todo: generalise to higher orders
            if ( N_ > 2)
            {
                // This has 20 terms
                notImplemented("Orders higher that quadratic not implemented");
            }
        }

        Eigen::MatrixXd W = Eigen::MatrixXd::Zero(Nn, Nn);

        for (label cI = 0; cI < Nn; cI++)
        {
            vector neiC = vector::zero;

            if (cI < curStencil.size())
            {
                const label neiCellID = curStencil[cI];
                neiC = CI[neiCellID];
            }
            else
            {
                // For boundary cells we need to add boundary face as
                // neigbour
                const label i = cI - curStencil.size();
                const label globalFaceID = cellBoundaryFaces[cellI][i];
                neiC = CfI[globalFaceID];
            }

            const vector& C = CI[cellI];
            const scalar d = mag(neiC - C);

            // Smoothing length
            const scalar dm = 2*maxDist;

            // Weight using radially symmetric exponential function
            const scalar sqrK = -pow(k_,2);
            const scalar w =
                (
                    Foam::exp(pow(d/dm, 2)*sqrK) - Foam::exp(sqrK)
                )/(1 - exp(sqrK));

            W(cI, cI) = w;
        }

        // Now when we have W and Q, next step is QR decomposition
        const Eigen::MatrixXd sqrtW = W.cwiseSqrt();
        const Eigen::MatrixXd Qhat = Q*sqrtW;
        Eigen::HouseholderQR<Eigen::MatrixXd> qr(Qhat.transpose());

        // Q and R matrices
        const Eigen::MatrixXd O = qr.householderQ();
        const Eigen::MatrixXd R = qr.matrixQR().triangularView<Eigen::Upper>();

        // B hat
        const Eigen::MatrixXd Bhat =
            sqrtW*Eigen::MatrixXd::Identity(W.rows(), W.cols());

        // Slice Rbar and Qbar, as we do not need full matrix
        const Eigen::MatrixXd Rbar = R.topLeftCorner(Np, Np);
        const Eigen::MatrixXd Qbar = O.leftCols(Np);

        // Solve to get A
        const Eigen::MatrixXd A =
            Rbar.colPivHouseholderQr().solve(Qbar.transpose()*Bhat);

        // To be aware of interpolation accuracy we need to control the
        // condition number
        if (calcConditionNumber_)
        {
            Eigen::JacobiSVD<Eigen::MatrixXd> svd
            (
                Rbar, Eigen::ComputeFullU | Eigen::ComputeFullV
            );
            Eigen::VectorXd singularValues = svd.singularValues();

            volScalarField& condNumber = *condNumberPtr_;
            condNumber[cellI] =
                singularValues(0)
               /(singularValues(singularValues.size() - 1) + VSMALL);
        }

        c[cellI].setCapacity(A.cols());
        cx[cellI].setCapacity(A.cols());
        cy[cellI].setCapacity(A.cols());
        cz[cellI].setCapacity(A.cols());

        Eigen::RowVectorXd cRow = A.row(0);
        Eigen::RowVectorXd cxRow = A.row(1);
        Eigen::RowVectorXd cyRow = A.row(2);
        Eigen::RowVectorXd czRow = A.row(3);

        for (label i = 0; i < A.cols(); ++i)
        {
            c[cellI].append(cRow(i));
            cx[cellI].append(cxRow(i));
            cy[cellI].append(cyRow(i));
            cz[cellI].append(czRow(i));
        }

        c[cellI].shrink();
        cx[cellI].shrink();
        cy[cellI].shrink();
        cz[cellI].shrink();
    }

    forAll(interpCoeffs, cellI)
    {
       const DynamicList<label>& curStencil = stencils[cellI];
       const label Nn = curStencil.size() + cellBoundaryFaces[cellI].size();

       interpCoeffs[cellI].setCapacity(Nn);
       interpGradCoeffs[cellI].setCapacity(Nn);

       for (label I = 0; I < Nn; I++)
       {
           interpCoeffs[cellI].append(c[cellI][I]);
           interpGradCoeffs[cellI].append
           (
               vector(cx[cellI][I], cy[cellI][I], cz[cellI][I])
           );
       }

       interpCoeffs[cellI].shrink();
       interpGradCoeffs[cellI].shrink();
    }

    Info<< "higherOrderGrad::calcCoeffs(): end" << endl;
}


const List<DynamicList<scalar>>& higherOrderGrad::interpCoeffs() const
{
    if (!interpCoeffsPtr_)
    {
        calcCoeffs();
    }

    return *interpCoeffsPtr_;
}


const List<DynamicList<vector>>& higherOrderGrad::interpGradCoeffs() const
{
    if (!interpGradCoeffsPtr_)
    {
        calcCoeffs();
    }

    return *interpGradCoeffsPtr_;
}


const List<DynamicList<label>>& higherOrderGrad::cellBoundaryFaces() const
{
    if (!cellBoundaryFacesPtr_)
    {
        calcCoeffs();
    }

    return *cellBoundaryFacesPtr_;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

higherOrderGrad::higherOrderGrad
(
    const fvMesh& mesh,
    const label N,
    const label nLayers,
    const label k,
    const label maxStencilSize,
    const bool calcConditionNumber
)
:
    mesh_(mesh),
    N_(N),
    nLayers_(nLayers),
    k_(k),
    maxStencilSize_(maxStencilSize),
    calcConditionNumber_(calcConditionNumber),
    condNumberPtr_(),
    stencilsPtr_(),
    interpCoeffsPtr_(),
    interpGradCoeffsPtr_(),
    cellBoundaryFacesPtr_()
{}


// * * * * * * * * * * * * * * * * Destructors * * * * * * * * * * * * * * * //

higherOrderGrad::~higherOrderGrad()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


tmp<volTensorField> higherOrderGrad::grad(const volVectorField& D)
{
    // Info<< "higherOrderGrad::grad(...)" << endl;

    const fvMesh& mesh = mesh_;

    // Prepare the return field
    tmp<volTensorField> tgradD
    (
        new volTensorField
        (
           IOobject
           (
               "grad(" + D.name() + ")",
               mesh.time().timeName(),
               mesh,
               IOobject::NO_READ,
               IOobject::AUTO_WRITE
           ),
           mesh,
           dimensionedTensor("0", dimless, Zero),
           "zeroGradient"
        )
    );
    volTensorField& gradD = tgradD.ref();

    const polyBoundaryMesh& boundaryMesh = mesh.boundaryMesh();
    const List<DynamicList<label>> stencils = this->stencils();
    const List<DynamicList<label>> cellBoundaryFaces =
        this->cellBoundaryFaces();
    //const List<DynamicList<scalar>>& interpCoeffs = this->interpCoeffs();
    const List<DynamicList<vector>>& interpGradCoeffs =
        this->interpGradCoeffs();

    forAll(stencils, cellI)
    {
       const DynamicList<label>& curStencil = stencils[cellI];
       const label Nn = curStencil.size() + cellBoundaryFaces[cellI].size();

       // Loop over stencil and multiply stencil cell values with
       // corresponding interpolation coefficient
       for(label cI = 0; cI < Nn; cI++)
       {
           if (cI < curStencil.size())
           {
               gradD[cellI] += interpGradCoeffs[cellI][cI]*D[curStencil[cI]];
           }
           else
           {
               const label i = cI - curStencil.size();
               const label globalFaceID = cellBoundaryFaces[cellI][i];

               vector boundaryD = vector::zero;

               forAll(boundaryMesh, patchI)
               {
                   const label start = boundaryMesh[patchI].start();
                   const label nFaces = boundaryMesh[patchI].nFaces();

                   if (globalFaceID >= start && globalFaceID < start + nFaces)
                   {
                       const label k = globalFaceID-start;
                       boundaryD = D.boundaryField()[patchI][k];
                   }
               }

               gradD[cellI] += interpGradCoeffs[cellI][cI]*boundaryD;
           }
       }
    }

    gradD.correctBoundaryConditions();

    // Info<< "higherOrderGrad::grad(...): end" << endl;

    return tgradD;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
