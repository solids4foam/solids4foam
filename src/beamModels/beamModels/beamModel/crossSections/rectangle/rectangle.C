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

InClass
    rectangle

\*---------------------------------------------------------------------------*/

#include "rectangle.H"

#include "volFields.H"
#include "fvc.H"
#include "addToRunTimeSelectionTable.H"
#include "primitivePatchInterpolation.H"
#include "scalarMatrices.H"
#include "primitiveFacePatch.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

  defineTypeNameAndDebug(rectangle, 0);
  addToRunTimeSelectionTable(crossSectionModel, rectangle, dictionary);

// * * * * * * * * * * * * Private Member Functions * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

rectangle::rectangle
(
    const word& name,
    const dictionary& dict
)
:
    crossSectionModel
    (
        name,
        dict
    ),
    crossSectionModelDict_(dict.subDict(name + "CrossSectionModelDict")),
    b_
    (
        readScalar(crossSectionModelDict_.lookup("b"))
    ),
    h_
    (
        readScalar(crossSectionModelDict_.lookup("h"))
    ),
    nx_
    (
        readScalar(crossSectionModelDict_.lookup("nx"))
    ),
    ny_
    (
        readScalar(crossSectionModelDict_.lookup("ny"))
    ),
    scalingArea_
    (
        crossSectionModelDict_.lookupOrDefault<scalar>("scalingArea", 1.0)
    ),
    scalingSecondMomentArea_
    (
        crossSectionModelDict_
        .lookupOrDefault<scalar>("scalingSecondMomentArea", 1.0)
    ),
    points_((nx_+1)*(ny_+1), vector::zero),
    quadFaces_(nx_*ny_, quadrilateral())
    // faces_(nx_*ny_, face(4))
{
    // Set local mesh points

    label npx = nx_+1;
    label npy = ny_+1;

    scalar dx = b_/nx_;
    scalar dy = h_/ny_;

    label pI = 0;
    for (label j=0; j<npy; j++)
    {
        for (label i=0; i<npx; i++)
        {
            points_[pI++] =
                vector(-b_/2 + i*dx, -h_/2 + j*dy, 0);
        }
    }

    // Set local mesh faces

    label fI = 0;
    for (label j=0; j<ny_; j++)
    {
        for (label i=0; i<nx_; i++)
        {
            FixedList<label, 4>& fpts = quadFaces_[fI++];
            // labelList& fpts = faces_[fI];

            fpts[0] = j*npx + i;
            fpts[1] = fpts[0] + 1;
            fpts[3] = (j+1)*npx + i;
            fpts[2] = fpts[3] + 1;
        }
    }

    // primitiveFacePatch patch(quadFaces_, points_);
}

// * * * * * * * * * * * * * * Member Functions * * * * * * * * * * * * * * * //

// tmp<vectorField> rectangle::greenLagrangianStrain
// (
//     const vector& Gamma,
//     const vector& K
// ) const
// {
//     tmp<vectorField> tE
//     (
//         new vectorField(points_.size(), vector::zero)
//     );

//     Field<scalar> avgE(6, 0);
//     avgE[0] = Gamma.x();
//     avgE[1] = Gamma.y();
//     avgE[2] = Gamma.z();
//     avgE[3] = K.x();
//     avgE[4] = K.y();
//     avgE[5] = K.z();

//     scalarRectangularMatrix A(3, 6, 0);

//     forAll(points_, pI)
//     {
//         calcA(A, points_[pI]);

//         for (label i=0; i<3; i++)
//         {
//             for (label j=0; j<6; j++)
//             {
//                 tE.ref()[pI].component(i) += A[i][j]*avgE[j];
//             }
//         }
//     }

//     return tE;
// }


// tmp<scalarField> rectangle::resultantTotalForce(const vectorField& S) const
// {
//     tmp<scalarField> trF
//     (
//         new scalarField(6, 0)
//     );

//     // List<scalarField> fS(4);
//     // List<scalarRectangularMatrix> A(4);

//     // forAll(quadFaces_, faceI)
//     // {
//     //     const quadrilateral& fpts = quadFaces_[fI++];

//     //     label i0 = fpts[0];
//     //     label i1 = fpts[1];
//     //     label i2 = fpts[2];
//     //     label i3 = fpts[3];

//     //     fS = scalarField(3, 0);
//     //     A = scalarRectangularMatrix(3, 6, 0);

//     //     fS

//     //     trF() += quadFaces_[faceI].integrate(points_, AT, fldS)
//     // }

//     return trF;
// }

// vector rectangle::resultantForce(const vectorField& S) const
// {
//     vector rF = vector::zero;

//     forAll(quadFaces_, faceI)
//     {
//         rF += quadFaces_[faceI].integral(points_, S);
//     }

//     return rF;
// }


// void rectangle::resultantTangentMatrices
// (
//     const diagTensorField& dSdE,
//     tensor& DQDGamma,
//     tensor& DQDK,
//     tensor& DMDGamma,
//     tensor& DMDK
// ) const
// {
//     // scalarRectangularMatrix A(3, 6, 0);
//     // scalarRectangularMatrix CT(3, 3, 0);

//     List<scalarRectangularMatrix> A(4);
//     // A = scalarRectangularMatrix(3, 6, 0);

//     List<scalarRectangularMatrix> CT(4);
//     // CT = scalarRectangularMatrix(3, 3, 0);

//     // List<scalarRectangularMatrix> ACA(4);
//     // ACA = scalarRectangularMatrix(6, 6, 0);

//     // scalarRectangularMatrix ACA(6, 6, 0);

//     scalarRectangularMatrix DSDE(6, 6, 0);

//     for (label fI=0; fI<quadFaces_.size(); fI++)
//     {
//         const quadrilateral& fpts = quadFaces_[fI];

//         label i0 = fpts[0];
//         label i1 = fpts[1];
//         label i2 = fpts[2];
//         label i3 = fpts[3];

//         A = scalarRectangularMatrix(3, 6, 0);
//         CT = scalarRectangularMatrix(3, 3, 0);

//         // scalarRectangularMatrix ACAm(6, 6, 0);

//         calcA(A[0], points_[i0]);
//         tensorToMatrix(CT[0], dSdE[i0]);
//         // multiply(ACA[0], A.T(), CT, A);
//         // multiply(ACA, A.T(), CT, A);
//         // ACAm = ACAm + 0.25*ACA;

//         calcA(A[1], points_[i1]);
//         tensorToMatrix(CT[1], dSdE[i1]);
//         // multiply(ACA[1], A.T(), CT, A);
//         // multiply(ACA, A.T(), CT, A);
//         // ACAm = ACAm + 0.25*ACA;

//         calcA(A[2], points_[i2]);
//         tensorToMatrix(CT[2], dSdE[i2]);
//         // multiply(ACA[2], A.T(), CT, A);
//         // multiply(ACA, A.T(), CT, A);
//         // ACAm = ACAm + 0.25*ACA;

//         calcA(A[3], points_[i3]);
//         tensorToMatrix(CT[3], dSdE[i3]);
//         // multiply(ACA[3], A.T(), CT, A);
//         // multiply(ACA, A.T(), CT, A);
//         // ACAm = ACAm + 0.25*ACA;

//         DSDE = DSDE + fpts.integral(points_, A, CT);

//         // DSDE = DSDE + fpts.integral(points_, ACA);
//         // DSDE = DSDE + dA*ACAm;
//     }

//     // for (label j=0; j<ny_; j++)
//     // {
//     //     for (label i=0; i<nx_; i++)
//     //     {
//     //         const quadrilateral& fpts = quadFaces_[fI++];

//     //         label i0 = j*npx + i;
//     //         label i1 = i0 + 1;
//     //         label i3 = (j+1)*npx + i;
//     //         label i2 = i3 + 1;

//     //         scalarRectangularMatrix ACAm(6, 6, 0);

//     //         calcA(A, points_[i0]);
//     //         tensorToMatrix(CT, dSdE[i0]);
//     //         multiply(ACA, A.T(), CT, A);
//     //         ACAm = ACAm + 0.25*ACA;

//     //         calcA(A, points_[i1]);
//     //         tensorToMatrix(CT, dSdE[i1]);
//     //         multiply(ACA, A.T(), CT, A);
//     //         ACAm = ACAm + 0.25*ACA;

//     //         calcA(A, points_[i2]);
//     //         tensorToMatrix(CT, dSdE[i2]);
//     //         multiply(ACA, A.T(), CT, A);
//     //         ACAm = ACAm + 0.25*ACA;

//     //         calcA(A, points_[i3]);
//     //         tensorToMatrix(CT, dSdE[i3]);
//     //         multiply(ACA, A.T(), CT, A);
//     //         ACAm = ACAm + 0.25*ACA;

//     //         DSDE = DSDE + dA*ACAm;
//     //     }
//     // }

//     // Info << DSDE << endl;

//     DQDGamma.xx() = DSDE[0][0];
//     DQDGamma.xy() = DSDE[0][1];
//     DQDGamma.xz() = DSDE[0][2];

//     DQDGamma.yx() = DSDE[1][0];
//     DQDGamma.yy() = DSDE[1][1];
//     DQDGamma.yz() = DSDE[1][2];

//     DQDGamma.zx() = DSDE[2][0];
//     DQDGamma.zy() = DSDE[2][1];
//     DQDGamma.zz() = DSDE[2][2];


//     DQDK.xx() = DSDE[0][3];
//     DQDK.xy() = DSDE[0][4];
//     DQDK.xz() = DSDE[0][5];

//     DQDK.yx() = DSDE[1][3];
//     DQDK.yy() = DSDE[1][4];
//     DQDK.yz() = DSDE[1][5];

//     DQDK.zx() = DSDE[2][3];
//     DQDK.zy() = DSDE[2][4];
//     DQDK.zz() = DSDE[2][5];

//     // Info << "DQDK: " << DQDK << endl;
//     DQDK = tensor::zero;

//     DMDGamma.xx() = DSDE[3][0];
//     DMDGamma.xy() = DSDE[3][1];
//     DMDGamma.xz() = DSDE[3][2];

//     DMDGamma.yx() = DSDE[4][0];
//     DMDGamma.yy() = DSDE[4][1];
//     DMDGamma.yz() = DSDE[4][2];

//     DMDGamma.zx() = DSDE[5][0];
//     DMDGamma.zy() = DSDE[5][1];
//     DMDGamma.zz() = DSDE[5][2];

//     // Info << "DMDGamma: " << DMDGamma << endl;
//     DMDGamma = tensor::zero;

//     DMDK.xx() = DSDE[3][3];
//     DMDK.xy() = DSDE[3][4];
//     DMDK.xz() = DSDE[3][5];

//     DMDK.yx() = DSDE[4][3];
//     DMDK.yy() = DSDE[4][4];
//     DMDK.yz() = DSDE[4][5];

//     DMDK.zx() = DSDE[5][3];
//     DMDK.zy() = DSDE[5][4];
//     DMDK.zz() = DSDE[5][5];
// }


// vector rectangle::resultantMoment(const vectorField& S) const
// {
//     vector rM = vector::zero;

//     scalarField Sx(S.component(0));
//     scalarField Sy(S.component(1));
//     scalarField Sz(S.component(2));

//     scalarField Px(points_.component(0));
//     scalarField Py(points_.component(1));

//     forAll(quadFaces_, faceI)
//     {
//         rM.x() += quadFaces_[faceI].integral(points_, Sz, Py);
//         rM.x() -= quadFaces_[faceI].integral(points_, Sy, Px);

//         rM.y() += quadFaces_[faceI].integral(points_, Sx, Px);

//         rM.z() -= quadFaces_[faceI].integral(points_, Sx, Py);
//     }

//     return rM;
// }


// void rectangle::writeVTK
// (
//     const fileName& fn,
//     const scalarField& sf0,
//     const scalarField& sf1,
//     const vectorField& vf
// ) const
// {
//     OFstream vtkFile(fn);

//     // Write header
//     vtkFile << "# vtk DataFile Version 3.0" << endl;
//     vtkFile << "2D scalar data" << endl;
//     vtkFile << "ASCII" << endl;
//     vtkFile << "DATASET POLYDATA" << endl;
//     // vtkFile << "DATASET UNSTRUCTURED_GRID" << endl;


//     // Write points
//     vtkFile << "\nPOINTS " << points_.size() << " float" << endl;
//     for (label i=0; i<points_.size(); i++)
//     {
//         vtkFile << points_[i].x() << " "
//                 << points_[i].y() << " "
//                 << points_[i].z() << endl;
//     }


//     // Write cells
//     label nCells = quadFaces_.size();
//     vtkFile << "\nPOLYGONS " << nCells << " " << 5*nCells << endl;
//     // vtkFile << "\nCELLS " << nCells << " " << 5*nCells << endl;
//     for (label i=0; i<nCells; i++)
//     {
//         vtkFile << quadFaces_[i].size() << " "
//                 << quadFaces_[i][0] << " "
//                 << quadFaces_[i][1] << " "
//                 << quadFaces_[i][2] << " "
//                 << quadFaces_[i][3] << endl;
//     }


//     // // Write cell types
//     // vtkFile << "\nCELL_TYPES " << nCells << endl;
//     // for (label i=0; i<nCells; i++)
//     // {
//     //     vtkFile << 9 << endl;
//     // }


//     //Write data
//     vtkFile << "\nPOINT_DATA " << points_.size() << endl;
//     vtkFile << "SCALARS eqEP float 1" << endl;
//     vtkFile << "LOOKUP_TABLE default" << endl;
//     for (label i=0; i<points_.size(); i++)
//     {
//         vtkFile << sf0[i] << endl;
//     }

//     vtkFile << "\nSCALARS eqS float 1" << endl;
//     vtkFile << "LOOKUP_TABLE default" << endl;
//     for (label i=0; i<points_.size(); i++)
//     {
//         vtkFile << sf1[i] << endl;
//     }

//     vtkFile << "\nVECTORS S float" << endl;
//     // vtkFile << "LOOKUP_TABLE default" << endl;
//     for (label i=0; i<points_.size(); i++)
//     {
//         vtkFile << float(vf[i].x()) << ' '
//                 << float(vf[i].y()) << ' '
//                 << float(vf[i].z()) << ' ';

//         if (i > 0 && (i % 3) == 0)
//         {
//             vtkFile << nl;
//         }
//     }
// }


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
