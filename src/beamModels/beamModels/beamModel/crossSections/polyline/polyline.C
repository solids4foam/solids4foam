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
    polyline

Description
    Incomplete and not compiled: see the Description in polyline.H.

\*---------------------------------------------------------------------------*/

#include "polyline.H"
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

  defineTypeNameAndDebug(polyline, 0);
  addToRunTimeSelectionTable(crossSectionModel, polyline, dictionary);

// * * * * * * * * * * * * Private Member Functions * * * * * * * * * * * * * //

scalar polyline::area() const
{
    scalar A=0;

    // Calculate cross-section area

    return A;
}


void polyline::calcSecondMomentOfArea(scalar& Ixx, scalar& Iyy)
{
    Ixx = 0;
    Iyy = 0;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

polyline::polyline
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
    circumferentialPoints_
    (
        pointField
        (
            crossSectionModelDict_.lookup
            (
                "polyLinePoints"
            )
        )
    ),
    minRadius_(),
    A_(),
    Ixx_(),
    Iyy_(),
    nCircumferentialPoints_(0),
    points_(0, vector::zero),
    quadFaces_(0, quadrilateral())
{
    Info << circumferentialPoints_ << endl;

    // Set number of circumferential points
    nCircumferentialPoints_ = circumferentialPoints_.size();

    if (nCircumferentialPoints_ == 0)
    {
        FatalErrorIn
	(
            "polyline::polyline(const word& name,const dictionary& dict)"
	)
            << "Zero size polyline" << abort(FatalError);
    }

    // Calculate minimal radius
    minRadius_ = 0;
    forAll(circumferentialPoints_, faceI)
    {
        minRadius_ += mag(circumferentialPoints_[faceI]);
    }
    minRadius_ /= nCircumferentialPoints_;

    // Calculate cross-section area
    A_ = area();

    // Calculate second moment of area
    calcSecondMomentOfArea(Ixx_, Iyy_);

    // Mesh points and faces
    // todo
}

// * * * * * * * * * * * * * * Member Functions * * * * * * * * * * * * * * * //

scalar polyline::radius(scalar theta) const
{
    scalar r=0;

    return r;
}



tmp<vectorField> polyline::circumferentialPoints() const
{
    tmp<vectorField> tCircumPoints
    (
        new vectorField(circumferentialPoints_)
    );

    return tCircumPoints;
}



tmp<vectorField> polyline::greenLagrangianStrain
(
    const vector& Gamma,
    const vector& K
) const
{
    tmp<vectorField> tE
    (
        new vectorField(points_.size(), vector::zero)
    );

    Field<scalar> avgE(6, 0);
    avgE[0] = Gamma.x();
    avgE[1] = Gamma.y();
    avgE[2] = Gamma.z();
    avgE[3] = K.x();
    avgE[4] = K.y();
    avgE[5] = K.z();

    scalarRectangularMatrix A(3, 6, 0);

    forAll(points_, pI)
    {
        calcA(A, points_[pI]);

        for (label i=0; i<3; i++)
        {
            for (label j=0; j<6; j++)
            {
                tE.ref()[pI].component(i) += A[i][j]*avgE[j];
            }
        }
    }

    return tE;
}


tmp<scalarField> polyline::resultantTotalForce(const vectorField& S) const
{
    tmp<scalarField> trF(new scalarField(6, 0));

    return trF;
}


vector polyline::resultantForce(const vectorField& S) const
{
    vector rF = vector::zero;

    forAll(quadFaces_, faceI)
    {
        rF += quadFaces_[faceI].integral(points_, S);
    }

    return rF;
}


void polyline::resultantTangentMatrices
(
    const diagTensorField& dSdE,
    tensor& DQDGamma,
    tensor& DQDK,
    tensor& DMDGamma,
    tensor& DMDK
) const
{
    // label nr = nRadialSegments_;
    // label nphi = nSegmentsPerSide_*nSides_;

    // label npr = nr+1;
    // label npphi = nphi+1;

    // scalar dphi = 2*M_PI/nphi;

    scalarRectangularMatrix A(3, 6, 0);
    scalarRectangularMatrix CT(3, 3, 0);
    scalarRectangularMatrix ACA(6, 6, 0);

    scalarRectangularMatrix DSDE(6, 6, 0);

    // label fI = 0;
    for (label fI=0; fI<quadFaces_.size(); fI++)
    {
        // for (label i=0; i<nr; i++)
        // {
            const quadrilateral& fpts = quadFaces_[fI];

            label i0 = fpts[0];
            label i1 = fpts[1];
            label i2 = fpts[2];
            label i3 = fpts[3];

            scalar dA = fpts.area(points_);

            // label i0 = j*npr + i;
            // label i1 = i0 + 1;
            // label i3 = (j+1)*npr + i;
            // label i2 = i2 + 1;

            // label i0 = j*npr + i;
            // label i1 = i0 + 1;
            // label i2 = (j+1)*npr + i;
            // label i3 = i2 + 1;

            scalarRectangularMatrix ACAm(6, 6, 0);

            calcA(A, points_[i0]);
            tensorToMatrix(CT, dSdE[i0]);
            multiply(ACA, A.T(), CT, A);
            ACAm = ACAm + 0.25*ACA;

            calcA(A, points_[i1]);
            tensorToMatrix(CT, dSdE[i1]);
            multiply(ACA, A.T(), CT, A);
            ACAm = ACAm + 0.25*ACA;

            calcA(A, points_[i2]);
            tensorToMatrix(CT, dSdE[i2]);
            multiply(ACA, A.T(), CT, A);
            ACAm = ACAm + 0.25*ACA;

            calcA(A, points_[i3]);
            tensorToMatrix(CT, dSdE[i3]);
            multiply(ACA, A.T(), CT, A);
            ACAm = ACAm + 0.25*ACA;

            DSDE = DSDE + dA*ACAm;
        // }
    }

    // // Info << DSDE << endl;

    DQDGamma.xx() = DSDE[0][0];
    // DQDGamma.xy() = DSDE[0][1];
    // DQDGamma.xz() = DSDE[0][2];
    // DQDGamma.yx() = DSDE[1][0];
    DQDGamma.yy() = DSDE[1][1];
    // DQDGamma.yz() = DSDE[1][2];
    // DQDGamma.zx() = DSDE[2][0];
    // DQDGamma.zy() = DSDE[2][1];
    DQDGamma.zz() = DSDE[2][2];

    DQDK.xx() = DSDE[0][3];
    DQDK.xy() = DSDE[0][4];
    DQDK.xz() = DSDE[0][5];
    DQDK.yx() = DSDE[1][3];
    DQDK.yy() = DSDE[1][4];
    DQDK.yz() = DSDE[1][5];
    DQDK.zx() = DSDE[2][3];
    DQDK.zy() = DSDE[2][4];
    DQDK.zz() = DSDE[2][5];

    DMDGamma.xx() = DSDE[3][0];
    DMDGamma.xy() = DSDE[3][1];
    DMDGamma.xz() = DSDE[3][2];
    DMDGamma.yx() = DSDE[4][0];
    DMDGamma.yy() = DSDE[4][1];
    DMDGamma.yz() = DSDE[4][2];
    DMDGamma.zx() = DSDE[5][0];
    DMDGamma.zy() = DSDE[5][1];
    DMDGamma.zz() = DSDE[5][2];

    DMDK.xx() = DSDE[3][3];
    // DMDK.xy() = DSDE[3][4];
    // DMDK.xz() = DSDE[3][5];
    // DMDK.yx() = DSDE[4][3];
    DMDK.yy() = DSDE[4][4];
    // DMDK.yz() = DSDE[4][5];
    // DMDK.zx() = DSDE[5][3];
    // DMDK.zy() = DSDE[5][4];
    DMDK.zz() = DSDE[5][5];
}


vector polyline::resultantMoment(const vectorField& S) const
{
    vector rM = vector::zero;

    scalarField Sx(S.component(0));
    scalarField Sy(S.component(1));
    scalarField Sz(S.component(2));

    scalarField Px(points_.component(0));
    scalarField Py(points_.component(1));

    forAll(quadFaces_, faceI)
    {
        rM.x() += quadFaces_[faceI].integral(points_, Sz, Py);
        rM.x() -= quadFaces_[faceI].integral(points_, Sy, Px);

        rM.y() += quadFaces_[faceI].integral(points_, Sx, Px);

        rM.z() -= quadFaces_[faceI].integral(points_, Sx, Py);
    }

    return rM;
}


void polyline::writeVTK
(
    const fileName& fn,
    const scalarField& sf0,
    const scalarField& sf1,
    const vectorField& vf
) const
{
    OFstream vtkFile(fn);

    // Write header
    vtkFile << "# vtk DataFile Version 3.0" << endl;
    vtkFile << "2D scalar data" << endl;
    vtkFile << "ASCII" << endl;
    vtkFile << "DATASET POLYDATA" << endl;
    // vtkFile << "DATASET UNSTRUCTURED_GRID" << endl;


    // Write points
    vtkFile << "\nPOINTS " << points_.size() << " float" << endl;
    for (label i=0; i<points_.size(); i++)
    {
        vtkFile << points_[i].x() << " "
                << points_[i].y() << " "
                << points_[i].z() << endl;
    }


    // Write cells
    label nCells = quadFaces_.size();
    vtkFile << "\nPOLYGONS " << nCells << " " << 5*nCells << endl;
    // vtkFile << "\nCELLS " << nCells << " " << 5*nCells << endl;
    for (label i=0; i<nCells; i++)
    {
        vtkFile << quadFaces_[i].size() << " "
                << quadFaces_[i][0] << " "
                << quadFaces_[i][1] << " "
                << quadFaces_[i][2] << " "
                << quadFaces_[i][3] << endl;
    }


    // // Write cell types
    // vtkFile << "\nCELL_TYPES " << nCells << endl;
    // for (label i=0; i<nCells; i++)
    // {
    //     vtkFile << 9 << endl;
    // }


    //Write data
    vtkFile << "\nPOINT_DATA " << points_.size() << endl;
    vtkFile << "SCALARS eqEP float 1" << endl;
    vtkFile << "LOOKUP_TABLE default" << endl;
    for (label i=0; i<points_.size(); i++)
    {
        vtkFile << sf0[i] << endl;
    }

    vtkFile << "\nSCALARS eqS float 1" << endl;
    vtkFile << "LOOKUP_TABLE default" << endl;
    for (label i=0; i<points_.size(); i++)
    {
        vtkFile << sf1[i] << endl;
    }

    vtkFile << "\nVECTORS S float" << endl;
    // vtkFile << "LOOKUP_TABLE default" << endl;
    for (label i=0; i<points_.size(); i++)
    {
        vtkFile << float(vf[i].x()) << ' '
                << float(vf[i].y()) << ' '
                << float(vf[i].z()) << ' ';

        if (i > 0 && (i % 3) == 0)
        {
            vtkFile << nl;
        }
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
