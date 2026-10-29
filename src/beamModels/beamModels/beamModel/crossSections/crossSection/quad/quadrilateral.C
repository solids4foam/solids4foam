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

#include "quadrilateral.H"
#include "IOstreams.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * //

quadrilateral::quadrilateral
(
    const label a,
    const label b,
    const label c,
    const label d
)
//:
    // face(4)
    // FixedList<label, 4>()
{
    operator[](0) = a;
    operator[](1) = b;
    operator[](2) = c;
    operator[](3) = d;
}

quadrilateral::quadrilateral()
:
    // face(4)
    FixedList<label, 4>()
    // FixedList<label, 4>(pts)
{}

quadrilateral::quadrilateral(const labelList& pts)
:
    // face(4)
    FixedList<label, 4>(pts)
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline pointField quadrilateral::points(const pointField& points) const
{
    pointField p(4);

    p[0] = points[operator[](0)];
    p[1] = points[operator[](1)];
    p[2] = points[operator[](2)];
    p[3] = points[operator[](3)];

    return p;
}

scalar quadrilateral::area(const pointField& pts) const
{
    scalar A = 0;

    pointField pt = points(pts);

    vector n0 =  0.5*((pt[1] - pt[0])^(pt[2] - pt[0]));
    vector n1 =  0.5*((pt[2] - pt[0])^(pt[3] - pt[0]));

    A = mag(n0) + mag(n1);

    return A;
}

// inline const vector& quadrilateral::a() const
// {
//     return a_;
// }

// inline const vector& quadrilateral::b() const
// {
//     return b_;
// }

// inline const vector& quadrilateral::c() const
// {
//     return c_;
// }

// inline const vector& quadrilateral::d() const
// {
//     return d_;
// }

scalar quadrilateral::Ni(label i, scalar xi, scalar eta) const
{
    scalar N = 0.0;

    switch(i)
    {
        case 0:
            N = (1.0 - xi)*(1 - eta)/4;
            break;
        case 1:
            N = (1.0 + xi)*(1 - eta)/4;
            break;
        case 2:
            N = (1.0 + xi)*(1 + eta)/4;
            break;
        case 3:
            N = (1.0 - xi)*(1 + eta)/4;
            break;
        default:
            FatalErrorIn
            (
                "scalar quadrilateral::Ni(label i, scalar xi, scalar eta) const"
            )   << "Index of quadrilateral shape function is out of range."
                << abort(FatalError);
    }

    return N;
}


scalar quadrilateral::dNidXi(label i, scalar xi, scalar eta) const
{
    scalar dN = 0.0;

    switch(i)
    {
        case 0:
            dN = -(1 - eta)/4;
            break;
        case 1:
            dN =  (1 - eta)/4;
            break;
        case 2:
            dN =  (1 + eta)/4;
            break;
        case 3:
            dN = -(1 + eta)/4;
            break;
        default:
            FatalErrorIn
            (
                "scalar quadrilateral::dNidXi(label i, scalar xi, scalar eta)"
            )   << "Index of quadrilateral shape function is out of range."
                << abort(FatalError);
    }

    return dN;
}


scalar quadrilateral::dNidEta(label i, scalar xi, scalar eta) const
{
    scalar dN = 0.0;

    switch(i)
    {
        case 0:
            dN = -(1.0 - xi)/4;
            break;
        case 1:
            dN = -(1.0 + xi)/4;
            break;
        case 2:
            dN =  (1.0 + xi)/4;
            break;
        case 3:
            dN =  (1.0 - xi)/4;
            break;
        default:
            FatalErrorIn
            (
                "scalar quadrilateral::dNidEta(label i, scalar xi, scalar eta)"
            )   << "Index of quadrilateral shape function is out of range."
                << abort(FatalError);
    }

    return dN;
}


scalar quadrilateral::detJacob
(
    scalar xi,
    scalar eta,
    const pointField& pts
) const
{
    scalar detJ = 0.0;

    tensor2D J = tensor2D::zero;

    pointField locPoints = points(pts);

    forAll (locPoints, pI)
    {
        J.xx() += locPoints[pI].x()*dNidXi(pI, xi, eta);
        J.xy() += locPoints[pI].y()*dNidXi(pI, xi, eta);
        J.yx() += locPoints[pI].x()*dNidEta(pI, xi, eta);
        J.yy() += locPoints[pI].y()*dNidEta(pI, xi, eta);

        // J.xx() += locPoints[pI].x()*dNidXi(pI, xi, eta);
        // J.xy() += locPoints[pI].x()*dNidEta(pI, xi, eta);
        // J.yx() += locPoints[pI].y()*dNidXi(pI, xi, eta);
        // J.yy() += locPoints[pI].y()*dNidEta(pI, xi, eta);
    }

    detJ = det(J);

    return detJ;
}


scalarRectangularMatrix quadrilateral::integral
(
    const Field<vector>& pts,
    const List<scalarRectangularMatrix>& A,
    const List<scalarRectangularMatrix>& CT
) const
{
    scalarRectangularMatrix integ
    (
        A[0].m(),
        A[0].m(),
        0
    );

    // Gauss integration points
    scalarField gXi(4, 0);
    gXi[0] = -1.0/::sqrt(3);
    gXi[1] =  1.0/::sqrt(3);
    gXi[2] =  1.0/::sqrt(3);
    gXi[3] = -1.0/::sqrt(3);

    scalarField gEta(4, 0);
    gEta[0] = -1.0/::sqrt(3);
    gEta[1] = -1.0/::sqrt(3);
    gEta[2] =  1.0/::sqrt(3);
    gEta[3] =  1.0/::sqrt(3);

    List<scalarRectangularMatrix> gA(4);
    gA = scalarRectangularMatrix(3, 6, 0);

    List<scalarRectangularMatrix> gCT(4);
    gCT = scalarRectangularMatrix(3, 3, 0);

    // pointField lp = points(pts);

    // pointField gp(4, vector::zero);
    // forAll(gCT, gpI)
    // {
    //     forAll(*this, pI)
    //     {
    //         gp[gpI] = gp[gpI]
    //           + Ni(pI, gXi[gpI], gEta[gpI])*lp[pI];
    //     }
    // }

    forAll(gCT, gpI)
    {
        forAll(*this, pI)
        {
            // gA[gpI][0][0] = 1;
            // gA[gpI][1][1] = 1;
            // gA[gpI][2][2] = 1;

            // gA[gpI][0][4] =  gp[gpI].x();
            // gA[gpI][0][5] = -gp[gpI].y();
            // gA[gpI][1][3] = -gp[gpI].x();
            // gA[gpI][2][3] =  gp[gpI].y();

            gA[gpI] = gA[gpI]
              + Ni(pI, gXi[gpI], gEta[gpI])*A[pI];
            gCT[gpI] = gCT[gpI]
              + Ni(pI, gXi[gpI], gEta[gpI])*CT[pI];
        }
    }

    Field<scalar> gMagJacob(4, 0);
    forAll(gMagJacob, gpI)
    {
        gMagJacob[gpI] =
            mag(detJacob(gXi[gpI], gEta[gpI], pts));
    }

    // scalar A = 0;
    forAll(gA, gpI)
    {
        scalarRectangularMatrix gACA(6, 6, 0);
        multiply(gACA, gA[gpI].T(), gCT[gpI], gA[gpI]);

        integ = integ + gMagJacob[gpI]*gACA;

        // A += gMagJacob[gpI];
    }

    // Note: dA = (dXi*dEta)*detJacob, dXi*dEta = 2
    // Info << "xxxxxxx: " << area(points) << ", " << A << endl;

    return integ;
}



// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

