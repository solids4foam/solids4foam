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
    annularSector

\*---------------------------------------------------------------------------*/

#include "annularSector.H"
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

  defineTypeNameAndDebug(annularSector, 0);
  addToRunTimeSelectionTable(crossSectionModel, annularSector, dictionary);

// * * * * * * * * * * * * Private Member Functions * * * * * * * * * * * * * //

scalar annularSector::area() const
{
    scalar A = 0;

    scalar halfSectorAngle = sectorAngle_/2;

    A = halfSectorAngle*(sqr(outerRadius_) - sqr(innerRadius_));

    Info << "A = " << A << endl;

    return A;
}

void annularSector::calcSecondMomentOfArea(scalar& Ixx, scalar& Iyy)
{
    Ixx = 0;
    Iyy = 0;

    scalar halfSectorAngle = sectorAngle_/2;

    // scalar A = halfSectorAngle*(sqr(outerRadius_) - sqr(innerRadius_));

    // Second moment of area for the sector about the geometric centroid

    // Centroid location (yC0 = 0; symmetry about the local x-axis)
    // scalar xC0 = (2.0/3.0)*(::sin(halfSectorAngle)/halfSectorAngle)
    //    *(
    //         (pow(outerRadius_,3) - pow(innerRadius_,3)
    //     )/
    //     (
    //         sqr(outerRadius_) - sqr(innerRadius_))
    //     );

    Ixx = 0.25*(pow(outerRadius_,4) - pow(innerRadius_,4))
          *(halfSectorAngle - 0.5*::sin(2*halfSectorAngle));

    Iyy = 0.25
       *(
            pow(outerRadius_,4) - pow(innerRadius_,4)
        )
       *(
            halfSectorAngle + 0.5*::sin(2*halfSectorAngle)
        );

    //- A*sqr(xC0);

    Info << "Local Ixx: " << Ixx << endl;
    Info << "Local Iyy: " << Iyy << endl;
}



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

annularSector::annularSector
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
    innerRadius_
    (
        readScalar(crossSectionModelDict_.lookup("innerRadius"))
    ),
    outerRadius_
    (
        readScalar(crossSectionModelDict_.lookup("outerRadius"))
    ),
    sectorAngle_
    (
        mag
        (
            readScalar
            (
                crossSectionModelDict_.lookup("sectorAngleInDeg")
            )*M_PI/180
        )
    ),
    nRadialSegmentsSector_
    (
        readInt(crossSectionModelDict_.lookup("nRadialSegmentsInSector"))
    ),
    A_(),
    Ixx_(),
    Iyy_(),
    nCircumferentialPoints_((nRadialSegmentsSector_+1)*2),
    points_((nRadialSegmentsSector_+1)*2, vector::zero),
    quadFaces_(nRadialSegmentsSector_, quadrilateral())
    // faces_(nx_*ny_, face(4))
{

    // Calculate cross-section area
    A_ = area();

    // Calculate second moment of area
    calcSecondMomentOfArea(Ixx_, Iyy_);

    // Set local mesh points
    label nr = nRadialSegmentsSector_;
    scalar dPhi = sectorAngle_/nr;
    scalar halfSectorAngle = sectorAngle_/2;
    label npr = nr + 1;

    label pI = 0;
    for(label i = 0; i < npr ; i++)
    {
        scalar x_inner = innerRadius_*::cos(-halfSectorAngle + i*dPhi);
        scalar y_inner = innerRadius_*::sin(-halfSectorAngle + i*dPhi);

        points_[pI++] = vector(0, x_inner, y_inner);
    }

    pI = 0;
    for(label i = 0; i < npr; i++)
    {
        scalar x_outer = outerRadius_*::cos(-halfSectorAngle + i*dPhi);
        scalar y_outer = outerRadius_*::sin(-halfSectorAngle + i*dPhi);

        points_[(pI++ + npr)] = vector(0, x_outer, y_outer);
    }

    Info << "Points: " << points_ << endl;

    // Set local mesh faces

    label fI = 0;
    for (label j = 0; j < nr; j++)
    {
            FixedList<label, 4>& fpts = quadFaces_[fI++];
            // labelList& fpts = faces_[fI];

            fpts[0] = j;
            fpts[1] = j + npr;
            fpts[2] = fpts[1] + 1;
            fpts[3] = fpts[0] + 1;

            Info << "Faces : " << fpts << endl;
    }

    // primitiveFacePatch patch(quadFaces_, points_);
}

// * * * * * * * * * * * * * * Member Functions * * * * * * * * * * * * * * * //

tmp<vectorField> annularSector::circumferentialPoints() const
{
    tmp<vectorField> tCircumPoints
    (
        new vectorField
        (
            nCircumferentialPoints_,
            vector::zero
        )
    );

    label nr = nRadialSegmentsSector_;
    scalar dPhi = sectorAngle_/nr;
    scalar halfSectorAngle = sectorAngle_/2;

    label npr = nr + 1;

    label gPointI = 0;
    for(label i = 0; i < npr; i++)
    {
        scalar x_outer = outerRadius_*::cos(halfSectorAngle - i*dPhi);
        scalar y_outer = outerRadius_*::sin(halfSectorAngle - i*dPhi);

        tCircumPoints.ref()[gPointI++]=
            vector(0, x_outer, y_outer);
        // tCircumPoints()[(gPointI++ + npr)]= vector(x_outer, y_outer, 0);
    }

    // label gPointI = 0;
    for(label i = 0; i < npr; i++)
    {
        scalar x_inner = innerRadius_*::cos(-halfSectorAngle + i*dPhi);
        scalar y_inner = innerRadius_*::sin(-halfSectorAngle + i*dPhi);

        tCircumPoints.ref()[gPointI++] =
            vector(0, x_inner, y_inner);
    }


    // Move cross-section to its centroid
    scalar xc =
        (2*::sin(halfSectorAngle)/3/halfSectorAngle)
       *(::pow(outerRadius_,3) - ::pow(innerRadius_,3))
       /(::pow(outerRadius_,2) - ::pow(innerRadius_,2));


    // Info << "xc: " << xc << endl;

    tCircumPoints.ref() -= vector(0, xc, 0);

    // Info << "tCircumPoints: " << tCircumPoints() << endl;

    return tCircumPoints;

}

// tmp<vectorField> annularSector::greenLagrangianStrain
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


// vector annularSector::resultantForce(const vectorField& S) const
// {
//     vector rF = vector::zero;

//     forAll(quadFaces_, faceI)
//     {
//         rF += quadFaces_[faceI].integral(points_, S);
//     }

//     return rF;
// }


// tmp<scalarField> annularSector::resultantTotalForce(const vectorField& S) const
// {
//     tmp<scalarField> trF(new scalarField(6, 0));

//     return trF;
// }


// void annularSector::resultantTangentMatrices
// (
//     const diagTensorField& dSdE,
//     tensor& DQDGamma,
//     tensor& DQDK,
//     tensor& DMDGamma,
//     tensor& DMDK
// ) const
// {
//     // label nr = nRadialSegments_;
//     // label nphi = nSegmentsPerSide_*nSides_;

//     // label npr = nr+1;
//     // label npphi = nphi+1;

//     // scalar dphi = 2*M_PI/nphi;

//     scalarRectangularMatrix A(3, 6, 0);
//     scalarRectangularMatrix CT(3, 3, 0);
//     scalarRectangularMatrix ACA(6, 6, 0);

//     scalarRectangularMatrix DSDE(6, 6, 0);

//     // label fI = 0;
//     forAll(quadFaces_, fI)
//     // for (label j=0; j<nphi; j++)
//     // {
//     //     for (label i=0; i<nr; i++)
//         {
//             const quadrilateral& fpts = quadFaces_[fI];

//             label i0 = fpts[0];
//             label i1 = fpts[1];
//             label i2 = fpts[2];
//             label i3 = fpts[3];

//             scalar dA = fpts.area(points_);

//             // label i0 = j*npr + i;
//             // label i1 = i0 + 1;
//             // label i3 = (j+1)*npr + i;
//             // label i2 = i2 + 1;

//             // label i0 = j*npr + i;
//             // label i1 = i0 + 1;
//             // label i2 = (j+1)*npr + i;
//             // label i3 = i2 + 1;

//             scalarRectangularMatrix ACAm(6, 6, 0);

//             calcA(A, points_[i0]);
//             tensorToMatrix(CT, dSdE[i0]);
//             multiply(ACA, A.T(), CT, A);
//             ACAm = ACAm + 0.25*ACA;

//             calcA(A, points_[i1]);
//             tensorToMatrix(CT, dSdE[i1]);
//             multiply(ACA, A.T(), CT, A);
//             ACAm = ACAm + 0.25*ACA;

//             calcA(A, points_[i2]);
//             tensorToMatrix(CT, dSdE[i2]);
//             multiply(ACA, A.T(), CT, A);
//             ACAm = ACAm + 0.25*ACA;

//             calcA(A, points_[i3]);
//             tensorToMatrix(CT, dSdE[i3]);
//             multiply(ACA, A.T(), CT, A);
//             ACAm = ACAm + 0.25*ACA;

//             DSDE = DSDE + dA*ACAm;
//         // }
//     }

//     // // Info << DSDE << endl;

//     DQDGamma.xx() = DSDE[0][0];
//     // DQDGamma.xy() = DSDE[0][1];
//     // DQDGamma.xz() = DSDE[0][2];
//     // DQDGamma.yx() = DSDE[1][0];
//     DQDGamma.yy() = DSDE[1][1];
//     // DQDGamma.yz() = DSDE[1][2];
//     // DQDGamma.zx() = DSDE[2][0];
//     // DQDGamma.zy() = DSDE[2][1];
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

//     DMDGamma.xx() = DSDE[3][0];
//     DMDGamma.xy() = DSDE[3][1];
//     DMDGamma.xz() = DSDE[3][2];
//     DMDGamma.yx() = DSDE[4][0];
//     DMDGamma.yy() = DSDE[4][1];
//     DMDGamma.yz() = DSDE[4][2];
//     DMDGamma.zx() = DSDE[5][0];
//     DMDGamma.zy() = DSDE[5][1];
//     DMDGamma.zz() = DSDE[5][2];

//     DMDK.xx() = DSDE[3][3];
//     // DMDK.xy() = DSDE[3][4];
//     // DMDK.xz() = DSDE[3][5];
//     // DMDK.yx() = DSDE[4][3];
//     DMDK.yy() = DSDE[4][4];
//     // DMDK.yz() = DSDE[4][5];
//     // DMDK.zx() = DSDE[5][3];
//     // DMDK.zy() = DSDE[5][4];
//     DMDK.zz() = DSDE[5][5];
// }


// vector annularSector::resultantMoment(const vectorField& S) const
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

// void annularSector::writeVTK
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
