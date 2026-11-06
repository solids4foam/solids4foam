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
    polygon

Description
    Not compiled yet: see the Description in polygon.H.

\*---------------------------------------------------------------------------*/

#include "polygon.H"
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

  defineTypeNameAndDebug(polygon, 0);
  addToRunTimeSelectionTable(crossSectionModel, polygon, dictionary);

// * * * * * * * * * * * * Private Member Functions * * * * * * * * * * * * * //

scalar polygon::area() const
{
    scalar A;

    // Calculate area of periodic patern
    scalar Aper = 0;
    if (filletAngle_ < SMALL)
    {
        scalar halfDTheta = M_PI/nSides_;

        scalar sideLength =
            2*minRadius_*::tan(halfDTheta);

        Aper = minRadius_*sideLength/2;
    }
    else
    {
        scalar DTheta = 2*M_PI/nSides_ - filletAngle_;
        scalar halfDTheta = DTheta/2;

        scalar sideLength =
            2*minRadius_*::tan(halfDTheta);

        Aper = minRadius_*sideLength/2;

        // Add fillet region
        scalar filletRadius = radius(filletAngle_/2);
        Aper += sqr(filletRadius)*0.5*filletAngle_;
    }

    A = Aper*nSides_;

    // Info << "A = " << A << endl;

    return A;
}


void polygon::calcSecondMomentOfArea(scalar& Ixx, scalar& Iyy)
{
    Ixx = 0;
    Iyy = 0;

    if (filletAngle_ < SMALL)
    {
        scalar DTheta = 2*M_PI/nSides_;
        scalar halfDTheta = DTheta/2;
        scalar sideLength = 2*minRadius_*::tan(halfDTheta);

        // Area of isoscales triangles
        scalar A = minRadius_*sideLength/2;

        // Second moment of area for isoscales triangle
        // whose symmetry axis is aligned with the global x-axis
        scalar IxC0 = minRadius_*::pow(sideLength,3)/48;
        scalar IyC0 = sideLength*::pow(minRadius_,3)/36;

        for (label i=0; i<nSides_; i++)
        {
            // Centroid of rotated isoscales triangle
            scalar thetaC = halfDTheta + i*DTheta;
            scalar rC = radius(thetaC);
            rC *= 2.0/3.0;

            scalar xC = rC*::cos(thetaC);
            scalar yC = rC*::sin(thetaC);

            // Coordinate rotation angle
            scalar rotAngle = -thetaC;

            // Second moment of area transformation due to coordinate rotaton
            scalar IxC = (IxC0+IyC0)/2 + (IxC0-IyC0)*::cos(2*rotAngle)/2;
            scalar IyC = (IxC0+IyC0)/2 - (IxC0-IyC0)*::cos(2*rotAngle)/2;

            // Parallel axis rule
            Ixx += IxC + A*sqr(yC);
            Iyy += IyC + A*sqr(xC);
        }
    }
    else
    {
        // Isoscales triangle parts
        {
            scalar DTheta = 2*M_PI/nSides_ - filletAngle_;
            scalar halfDTheta = DTheta/2;

            scalar sideLength =
                2*minRadius_*::tan(halfDTheta);

            // Area of isoscales triangles
            scalar A = minRadius_*sideLength/2;

            // Second moment of area for isoscales triangle
            // whose symmetry axis is aligned with the global x-axis
            scalar IxC0 = minRadius_*::pow(sideLength,3)/48;
            scalar IyC0 = sideLength*::pow(minRadius_,3)/36;

            for (label i=0; i<nSides_; i++)
            {
                // Centroid of rotated isoscales triangle
                scalar thetaC = halfDTheta + i*DTheta;
                scalar rC = radius(thetaC);
                rC *= 2.0/3.0;

                scalar xC = rC*::cos(thetaC);
                scalar yC = rC*::sin(thetaC);

                // Coordinate rotation angle
                scalar rotAngle = -thetaC;

                // Second moment of area transformation
                // due to coordinate rotaton
                scalar IxC = (IxC0+IyC0)/2 + (IxC0-IyC0)*::cos(2*rotAngle)/2;
                scalar IyC = (IxC0+IyC0)/2 - (IxC0-IyC0)*::cos(2*rotAngle)/2;

                // Parallel axis rule
                Ixx += IxC + A*sqr(yC);
                Iyy += IyC + A*sqr(xC);
            }
        }

        // Circular sector (fillet) parts
        {
            scalar DTheta = 2*M_PI/nSides_;

            scalar halfFilletAngle = filletAngle_/2;

            // Fillet radius
            scalar filletRadius = radius(filletAngle_/2);

            // Area of circular sector
            scalar A = sqr(filletRadius)*halfFilletAngle;

            // Second moment of area for circular sector
            // whose symmetry axis is aligned with the global x-axis
            // https://structx.com/Shape_Formulas_004.html
            scalar IxC0 =
                ::pow(filletRadius, 4)
               *(halfFilletAngle - ::sin(2*halfFilletAngle)/2)/4;
            scalar IyC0 =
                ::pow(filletRadius, 4)
               *(halfFilletAngle + ::sin(2*halfFilletAngle)/2)/4
              - 4*::pow(filletRadius, 4)*pow(::sin(halfFilletAngle),2)
               /(9*halfFilletAngle);

            for (label i=0; i<nSides_; i++)
            {
                // Centroid of rotated isoscales triangle
                scalar thetaC = i*DTheta;
                scalar rC =
                    filletRadius
                   *(2.0/3.0)*(::sin(halfFilletAngle)/halfFilletAngle);

                scalar xC = rC*::cos(thetaC);
                scalar yC = rC*::sin(thetaC);

                // Coordinate rotation angle
                scalar rotAngle = -thetaC;

                // Second moment of area transformation
                // due to coordinate rotaton
                scalar IxC = (IxC0+IyC0)/2 + (IxC0-IyC0)*::cos(2*rotAngle)/2;
                scalar IyC = (IxC0+IyC0)/2 - (IxC0-IyC0)*::cos(2*rotAngle)/2;

                // Parallel axis rule
                Ixx += IxC + A*sqr(yC);
                Iyy += IyC + A*sqr(xC);
            }
        }
    }

    // test
    if (false)
    {
        scalar DTheta = 2*M_PI/nSides_;
        scalar halfDTheta = DTheta/2;
        scalar sideLength = 2*minRadius_*::tan(halfDTheta);

        Info << Ixx_ << ", " << Iyy_ << endl;
        Info << (5*::sqrt(3)/16)*::pow(sideLength,4) << endl;
        Info << M_PI*pow(minRadius_, 4)/4 << endl;
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

polygon::polygon
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
    nSides_
    (
        readInt(crossSectionModelDict_.lookup("nSides"))
    ),
    minRadius_
    (
        readScalar
        (
            crossSectionModelDict_.lookup("minRadius")
        )
    ),
    filletAngle_
    (
        mag
        (
            readScalar
            (
                crossSectionModelDict_.lookup("filletAngleDeg")
            )*M_PI/180
        )
    ),
    A_(),
    Ixx_(),
    Iyy_(),
    nCircumferentialPoints_(0),
    nSegmentsPerSide_
    (
        readInt(crossSectionModelDict_.lookup("nSegmentsPerSide"))
    ),
    nRadialSegments_
    (
        readInt(crossSectionModelDict_.lookup("nRadialSegments"))
    ),
    points_((nSegmentsPerSide_*nSides_+1)*(nRadialSegments_+1), vector::zero),
    quadFaces_(nSegmentsPerSide_*nRadialSegments_*nSides_, quadrilateral())
{

    // Set number of circumferential points
    if (filletAngle_<SMALL)
    {
        nCircumferentialPoints_ = nSides_;
    }
    else
    {
        nCircumferentialPoints_ = nSides_*5;
    }

    // Calculate cross-section area
    A_ = area();

    // Calculate second moment of area
    calcSecondMomentOfArea(Ixx_, Iyy_);


    // Set local mesh points

    label nr = nRadialSegments_;
    label nphi = nSegmentsPerSide_*nSides_;

    label npr = nr+1;
    label npphi = nphi+1;

    scalar dphi = 2*M_PI/nphi;

    label pI = 0;
    for (label j=0; j<npphi; j++)
    {
        scalar phi = j*dphi;
        scalar R = radius(phi);
        scalar dr = R/nr;

        for (label i=0; i<npr; i++)
        {
            scalar r = i*dr;

            scalar x = r*cos(phi);
            scalar y = r*sin(phi);

            points_[pI++] = vector(x, y, 0);
        }
    }


    // Set local mesh faces

    label fI = 0;
    for (label j=0; j<nphi; j++)
    {
        for (label i=0; i<nr; i++)
        {
            FixedList<label, 4>& fpts = quadFaces_[fI++];

            fpts[0] = j*npr + i;
            fpts[1] = fpts[0] + 1;
            fpts[3] = (j+1)*npr + i;
            fpts[2] = fpts[3] + 1;
        }
    }

    // scalarField sf0(points_.size(), 0);
    // scalarField sf1(points_.size(), 0);
    // vectorField vf(points_.size(), vector::zero);
    // fileName fn("test.vtk");
    // writeVTK(fn, sf0, sf1, vf);

    // primitiveFacePatch patch(quadFaces_, points_);
}

// * * * * * * * * * * * * * * Member Functions * * * * * * * * * * * * * * * //

scalar polygon::radius(scalar theta) const
{
    // Polygon entral angel
    scalar DTheta = 2*M_PI/nSides_; //32;

    // Half of central angle
    scalar DPhi = DTheta/2;

    // Distance of the side from the pole
    scalar p = minRadius_;
    // scalar p = sideLength_/(2*::tan(DPhi));

    // Angle inside periodic patern
    label nPerPatern = floor(theta/DTheta);
    scalar thetaPer = theta - nPerPatern*DTheta;

    // Polar equation of line (first side)
    // (http://www.nabla.hr/Z_MemoHU-015.htm)
    scalar r = p/::cos(thetaPer-DPhi);

    if (mag(filletAngle_) > SMALL)
    {
        scalar halfFilletAngle = filletAngle_/2;
        scalar filletRadius = p/::cos(halfFilletAngle-DPhi);

        if
        (
            (thetaPer <= halfFilletAngle)
         || (thetaPer >= (DTheta-halfFilletAngle))
        )
        {
            r = filletRadius;
        }
    }

    return r;
}



tmp<vectorField> polygon::circumferentialPoints() const
{
    tmp<vectorField> tCircumPoints
    (
        new vectorField
        (
            nCircumferentialPoints_,
            vector::zero
        )
    );

    // Central angel
    scalar DTheta = -2*M_PI/nSides_; // Note minus sign
                                     // Important for correct orientation
                                     // of cell-faces

    if (filletAngle_ < SMALL)
    {
        for (label i=0; i<nSides_; i++)
        {
            scalar theta = i*DTheta;
            scalar r = radius(theta);

            scalar x = r*::cos(theta);
            scalar y = r*::sin(theta);

            tCircumPoints.ref()[i] = vector(x, y, 0);
        }
    }
    else
    {
        label gPointI = 0;
        for (label i=0; i<nSides_; i++)
        {
            // 1st point
            {
                scalar theta = i*DTheta;
                scalar r = radius(theta);
                scalar x = r*::cos(theta);
                scalar y = r*::sin(theta);
                tCircumPoints.ref()[gPointI++] = vector(x, y, 0);
            }

            // 2nd point
            {
                scalar theta = i*DTheta - 0.25*filletAngle_;
                scalar r = radius(theta);
                scalar x = r*::cos(theta);
                scalar y = r*::sin(theta);
                tCircumPoints.ref()[gPointI++] = vector(x, y, 0);
            }

            // 3rd point
            {
                scalar theta = i*DTheta - 0.5*filletAngle_;
                scalar r = radius(theta);
                scalar x = r*::cos(theta);
                scalar y = r*::sin(theta);
                tCircumPoints.ref()[gPointI++] = vector(x, y, 0);
            }

            // 4th point
            {
                scalar theta = (i+1)*DTheta + 0.5*filletAngle_;
                scalar r = radius(theta);
                scalar x = r*::cos(theta);
                scalar y = r*::sin(theta);
                tCircumPoints.ref()[gPointI++] = vector(x, y, 0);
            }

            // 5th point
            {
                scalar theta = (i+1)*DTheta + 0.25*filletAngle_;
                scalar r = radius(theta);
                scalar x = r*::cos(theta);
                scalar y = r*::sin(theta);
                tCircumPoints.ref()[gPointI++] = vector(x, y, 0);
            }
        }
    }

    // Info << tCircumPoints() << endl;

    return tCircumPoints;
}



tmp<vectorField> polygon::greenLagrangianStrain
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


tmp<scalarField> polygon::resultantTotalForce(const vectorField& S) const
{
    tmp<scalarField> trF(new scalarField(6, 0));

    return trF;
}


vector polygon::resultantForce(const vectorField& S) const
{
    vector rF = vector::zero;

    forAll(quadFaces_, faceI)
    {
        rF += quadFaces_[faceI].integral(points_, S);
    }

    return rF;
}


void polygon::resultantTangentMatrices
(
    const diagTensorField& dSdE,
    tensor& DQDGamma,
    tensor& DQDK,
    tensor& DMDGamma,
    tensor& DMDK
) const
{
    label nr = nRadialSegments_;
    label nphi = nSegmentsPerSide_*nSides_;

    // label npr = nr+1;
    // label npphi = nphi+1;

    // scalar dphi = 2*M_PI/nphi;

    scalarRectangularMatrix A(3, 6, 0);
    scalarRectangularMatrix CT(3, 3, 0);
    scalarRectangularMatrix ACA(6, 6, 0);

    scalarRectangularMatrix DSDE(6, 6, 0);

    label fI = 0;
    for (label j=0; j<nphi; j++)
    {
        for (label i=0; i<nr; i++)
        {
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
        }
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


vector polygon::resultantMoment(const vectorField& S) const
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


void polygon::writeVTK
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
