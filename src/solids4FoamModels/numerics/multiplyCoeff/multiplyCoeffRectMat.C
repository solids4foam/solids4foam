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

#include "multiplyCoeffRectMat.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::multiplyCoeffRectMat
(
    tensor& coeff,
    const vector& Sf,
    const RectangularMatrix<scalar>& C,
    const vector& g
)
{

    const label XX = 0;
    const label YY = 1;
    const label ZZ = 2;
    const label XY = 3;
    const label YZ = 4;
    const label XZ = 5;

//    Define indexes for readability.
//    scalar XX = 0;
//    scalar YY = 1;
//    scalar ZZ = 2;
//    scalar XY = 3;
//    scalar YZ = 4;
//    scalar ZX = 5;

    // Index notation
    // coeff_ij = Sf_m C_mikl g_k delta_lj
    // where delta is kronecker delta

    // Term is only non-zero for l = j from 1 to 3:
    // coeff_i1 = Sf_m C_mik1 g_k delta_11
    // coeff_i2 = Sf_m C_mik2 g_k delta_22
    // coeff_i3 = Sf_m C_mik3 g_k delta_33

    // coeff_11 = Sf_m C_m1k1 g_k
    //
    //          = Sf_1 C_11k1 g_k
    //            + Sf_2 C_21k1 g_k
    //            + Sf_3 C_31k1 g_k
    //
    //          = Sf_1 C_1111 g_1
    //            + Sf_1 C_1121 g_2
    //            + Sf_1 C_1131 g_3

    //            + Sf_2 C_2111 g_1
    //            + Sf_2 C_2121 g_2
    //            + Sf_2 C_2131 g_3

    //            + Sf_3 C_3111 g_1
    //            + Sf_3 C_3121 g_2
    //            + Sf_3 C_3131 g_3

    // Note: C_2121 == C_1212

    coeff[tensor::XX] =
        (Sf[vector::X]*C(XX,tensor::XX))*g[vector::X]
      + (Sf[vector::X]*C(XX,tensor::YX))*g[vector::Y]
      + (Sf[vector::X]*C(XX,tensor::ZX))*g[vector::Z]

      + (Sf[vector::Y]*C(XY,tensor::XX))*g[vector::X]
      + (Sf[vector::Y]*C(XY,tensor::YX))*g[vector::Y]
      + (Sf[vector::Y]*C(XY,tensor::ZX))*g[vector::Z]

      + (Sf[vector::Z]*C(XZ,tensor::XX))*g[vector::X]
      + (Sf[vector::Z]*C(XZ,tensor::YX))*g[vector::Y]
      + (Sf[vector::Z]*C(XZ,tensor::ZX))*g[vector::Z];

//    coeff[tensor::XX] =
//        Sf[vector::X]*C(XX,XX)*g[vector::X]
//      + Sf[vector::X]*C(XX,XY)*g[vector::Y]
//      + Sf[vector::X]*C(XX,ZX)*g[vector::Z]
//
//      + Sf[vector::Y]*C(XY,XX)*g[vector::X]
//      + Sf[vector::Y]*C(XY,XY)*g[vector::Y]
//      + Sf[vector::Y]*C(XY,ZX)*g[vector::Z]
//
//      + Sf[vector::Z]*C(ZX,XX)*g[vector::X]
//      + Sf[vector::Z]*C(ZX,XY)*g[vector::Y]
//      + Sf[vector::Z]*C(ZX,ZX)*g[vector::Z];

    // Similarly the other components can be calculated as
    // coeff_21 = Sf_m C_m2k1 g_k
    //
    //          = Sf_1 C_12k1 g_k
    //            + Sf_2 C_22k1 g_k
    //            + Sf_3 C_32k1 g_k

    //          = Sf_1 C_1211 g_1
    //            + Sf_1 C_1221 g_2
    //            + Sf_1 C_1231 g_3

    //            + Sf_2 C_2211 g_1
    //            + Sf_2 C_2221 g_2
    //            + Sf_2 C_2231 g_3

    //            + Sf_3 C_3211 g_1
    //            + Sf_3 C_3221 g_2
    //            + Sf_3 C_3231 g_3

    // Note: C_2211 == C_1122 and C_1221 == C_2112 == C_1212 == C_2121

    coeff[tensor::YX] =
        (Sf[vector::X]*C(XY,tensor::XX))*g[vector::X]
      + (Sf[vector::X]*C(XY,tensor::YX))*g[vector::Y]
      + (Sf[vector::X]*C(XY,tensor::ZX))*g[vector::Z]

      + (Sf[vector::Y]*C(YY,tensor::XX))*g[vector::X]
      + (Sf[vector::Y]*C(YY,tensor::YX))*g[vector::Y]
      + (Sf[vector::Y]*C(YY,tensor::ZX))*g[vector::Z]

      + (Sf[vector::Z]*C(YZ,tensor::XX))*g[vector::X]
      + (Sf[vector::Z]*C(YZ,tensor::YX))*g[vector::Y]
      + (Sf[vector::Z]*C(YZ,tensor::ZX))*g[vector::Z];

//    coeff[tensor::YX] =
//        Sf[vector::X]*C(XY,XX)*g[vector::X]
//      + Sf[vector::X]*C(XY,XY)*g[vector::Y]
//      + Sf[vector::X]*C(XY,ZX)*g[vector::Z]
//
//      + Sf[vector::Y]*C(YY,XX)*g[vector::X]
//      + Sf[vector::Y]*C(YY,XY)*g[vector::Y]
//      + Sf[vector::Y]*C(YY,ZX)*g[vector::Z]
//
//      + Sf[vector::Z]*C(YZ,XX)*g[vector::X]
//      + Sf[vector::Z]*C(YZ,XY)*g[vector::Y]
//      + Sf[vector::Z]*C(YZ,ZX)*g[vector::Z];

    // coeff_31 = Sf_m C_m3k1 g_k
    //
    //          = Sf_1 C_13k1 g_k
    //            + Sf_2 C_23k1 g_k
    //            + Sf_3 C_33k1 g_k

    //          = Sf_1 C_1311 g_1
    //            + Sf_1 C_1321 g_2
    //            + Sf_1 C_1331 g_3

    //            + Sf_2 C_2311 g_1
    //            + Sf_2 C_2321 g_2
    //            + Sf_2 C_2331 g_3

    //            + Sf_3 C_3311 g_1
    //            + Sf_3 C_3321 g_2
    //            + Sf_3 C_3331 g_3
    //
    // Note: C_1331 == C_1313 and C_3311 == C_1133

    coeff[tensor::ZX] =
        (Sf[vector::X]*C(XZ,tensor::XX))*g[vector::X]
      + (Sf[vector::X]*C(XZ,tensor::YX))*g[vector::Y]
      + (Sf[vector::X]*C(XZ,tensor::ZX))*g[vector::Z]

      + (Sf[vector::Y]*C(YZ,tensor::XX))*g[vector::X]
      + (Sf[vector::Y]*C(YZ,tensor::YX))*g[vector::Y]
      + (Sf[vector::Y]*C(YZ,tensor::ZX))*g[vector::Z]

      + (Sf[vector::Z]*C(ZZ,tensor::XX))*g[vector::X]
      + (Sf[vector::Z]*C(ZZ,tensor::YX))*g[vector::Y]
      + (Sf[vector::Z]*C(ZZ,tensor::ZX))*g[vector::Z];

//    coeff[tensor::ZX] =
//        Sf[vector::X]*C(ZX,XX)*g[vector::X]
//      + Sf[vector::X]*C(ZX,XY)*g[vector::Y]
//      + Sf[vector::X]*C(ZX,ZX)*g[vector::Z]
//
//      + Sf[vector::Y]*C(YZ,XX)*g[vector::X]
//      + Sf[vector::Y]*C(YZ,XY)*g[vector::Y]
//      + Sf[vector::Y]*C(YZ,ZX)*g[vector::Z]
//
//      + Sf[vector::Z]*C(ZZ,XX)*g[vector::X]
//      + Sf[vector::Z]*C(ZZ,XY)*g[vector::Y]
//      + Sf[vector::Z]*C(ZZ,ZX)*g[vector::Z];

    // coeff_i2 = Sf_m C_mik2 g_k delta_22
    //
    //          = Sf_1 C_1ik2 g_k
    //            + Sf_2 C_2ik2 g_k
    //            + Sf_3 C_3ik2 g_k
    //
    //          = Sf_1 C_1i12 g_1
    //            + Sf_1 C_1i22 g_2
    //            + Sf_1 C_1i32 g_3
    //
    //            + Sf_2 C_2i12 g_1
    //            + Sf_2 C_2i22 g_2
    //            + Sf_2 C_2i32 g_3
    //
    //            + Sf_3 C_3i12 g_1
    //            + Sf_3 C_3i22 g_2
    //            + Sf_3 C_3i32 g_3

    // coeff_12 = Sf_m C_m1k2 g_k
    //
    //          = Sf_1 C_11k2 g_k
    //            + Sf_2 C_21k2 g_k
    //            + Sf_3 C_31k2 g_k

    //          = Sf_1 C_1112 g_1
    //            + Sf_1 C_1122 g_2
    //            + Sf_1 C_1132 g_3

    //            + Sf_2 C_2112 g_1
    //            + Sf_2 C_2122 g_2
    //            + Sf_2 C_2132 g_3

    //            + Sf_3 C_3112 g_1
    //            + Sf_3 C_3122 g_2
    //            + Sf_3 C_3132 g_3

    // coeff_12 = Sf_2 C_2112 g_1 + Sf_1 C_1122 g_2

    coeff[tensor::XY] =
        (Sf[vector::X]*C(XX,tensor::XY))*g[vector::X]
      + (Sf[vector::X]*C(XX,tensor::YY))*g[vector::Y]
      + (Sf[vector::X]*C(XX,tensor::ZY))*g[vector::Z]

      + (Sf[vector::Y]*C(XY,tensor::XY))*g[vector::X]
      + (Sf[vector::Y]*C(XY,tensor::YY))*g[vector::Y]
      + (Sf[vector::Y]*C(XY,tensor::ZY))*g[vector::Z]

      + (Sf[vector::Z]*C(XZ,tensor::XY))*g[vector::X]
      + (Sf[vector::Z]*C(XZ,tensor::YY))*g[vector::Y]
      + (Sf[vector::Z]*C(XZ,tensor::ZY))*g[vector::Z];

//    coeff[tensor::XY] =
//        Sf[vector::X]*C(XX,XY)*g[vector::X]
//      + Sf[vector::X]*C(XX,YY)*g[vector::Y]
//      + Sf[vector::X]*C(XX,YZ)*g[vector::Z]
//
//      + Sf[vector::Y]*C(XY,XY)*g[vector::X]
//      + Sf[vector::Y]*C(XY,YY)*g[vector::Y]
//      + Sf[vector::Y]*C(XY,YZ)*g[vector::Z]
//
//      + Sf[vector::Z]*C(ZX,XY)*g[vector::X]
//      + Sf[vector::Z]*C(ZX,YY)*g[vector::Y]
//      + Sf[vector::Z]*C(ZX,YZ)*g[vector::Z];

    // coeff_22 = Sf_m C_m2k2 g_k
    //
    //          = Sf_1 C_12k2 g_k
    //            + Sf_2 C_22k2 g_k
    //            + Sf_3 C_32k2 g_k

    //          = Sf_1 C_1212 g_1
    //            + Sf_1 C_1222 g_2
    //            + Sf_1 C_1232 g_3
    //
    //            + Sf_2 C_2212 g_1
    //            + Sf_2 C_2222 g_2
    //            + Sf_2 C_2232 g_3
    //
    //            + Sf_3 C_3212 g_1
    //            + Sf_3 C_3222 g_2
    //            + Sf_3 C_3232 g_3
    // Note: C_3232 == C_2323

    coeff[tensor::YY] =
        (Sf[vector::X]*C(XY,tensor::XY))*g[vector::X]
      + (Sf[vector::X]*C(XY,tensor::YY))*g[vector::Y]
      + (Sf[vector::X]*C(XY,tensor::ZY))*g[vector::Z]

      + (Sf[vector::Y]*C(YY,tensor::XY))*g[vector::X]
      + (Sf[vector::Y]*C(YY,tensor::YY))*g[vector::Y]
      + (Sf[vector::Y]*C(YY,tensor::ZY))*g[vector::Z]

      + (Sf[vector::Z]*C(YZ,tensor::XY))*g[vector::X]
      + (Sf[vector::Z]*C(YZ,tensor::YY))*g[vector::Y]
      + (Sf[vector::Z]*C(YZ,tensor::ZY))*g[vector::Z];

//    coeff[tensor::YY] =
//        Sf[vector::X]*C(XY,XY)*g[vector::X]
//      + Sf[vector::X]*C(XY,YY)*g[vector::Y]
//      + Sf[vector::X]*C(XY,YZ)*g[vector::Z]
//
//      + Sf[vector::Y]*C(YY,XY)*g[vector::X]
//      + Sf[vector::Y]*C(YY,YY)*g[vector::Y]
//      + Sf[vector::Y]*C(YY,YZ)*g[vector::Z]
//
//      + Sf[vector::Z]*C(YZ,XY)*g[vector::X]
//      + Sf[vector::Z]*C(YZ,YY)*g[vector::Y]
//      + Sf[vector::Z]*C(YZ,YZ)*g[vector::Z];

    // coeff_32 = Sf_m C_m3k2 g_k
    //
    //          = Sf_1 C_13k2 g_k
    //            + Sf_2 C_23k2 g_k
    //            + Sf_3 C_33k2 g_k

    //          = Sf_1 C_1312 g_1
    //            + Sf_1 C_1322 g_2
    //            + Sf_1 C_1332 g_3
    //
    //            + Sf_2 C_2312 g_1
    //            + Sf_2 C_2322 g_2
    //            + Sf_2 C_2332 g_3
    //
    //            + Sf_3 C_3312 g_1
    //            + Sf_3 C_3322 g_2
    //            + Sf_3 C_3332 g_3
    //
    // Note: C_2332 == C_2323 and C_3322 == C_2233

    coeff[tensor::ZY] =
        (Sf[vector::X]*C(XZ,tensor::XY))*g[vector::X]
      + (Sf[vector::X]*C(XZ,tensor::YY))*g[vector::Y]
      + (Sf[vector::X]*C(XZ,tensor::ZY))*g[vector::Z]

      + (Sf[vector::Y]*C(YZ,tensor::XY))*g[vector::X]
      + (Sf[vector::Y]*C(YZ,tensor::YY))*g[vector::Y]
      + (Sf[vector::Y]*C(YZ,tensor::ZY))*g[vector::Z]

      + (Sf[vector::Z]*C(ZZ,tensor::XY))*g[vector::X]
      + (Sf[vector::Z]*C(ZZ,tensor::YY))*g[vector::Y]
      + (Sf[vector::Z]*C(ZZ,tensor::ZY))*g[vector::Z];

//    coeff[tensor::ZY] =
//        Sf[vector::X]*C(ZX,XY)*g[vector::X]
//      + Sf[vector::X]*C(ZX,YY)*g[vector::Y]
//      + Sf[vector::X]*C(ZX,YZ)*g[vector::Z]
//
//      + Sf[vector::Y]*C(YZ,XY)*g[vector::X]
//      + Sf[vector::Y]*C(YZ,YY)*g[vector::Y]
//      + Sf[vector::Y]*C(YZ,YZ)*g[vector::Z]
//
//      + Sf[vector::Z]*C(ZZ,XY)*g[vector::X]
//      + Sf[vector::Z]*C(ZZ,YY)*g[vector::Y]
//      + Sf[vector::Z]*C(ZZ,YZ)*g[vector::Z];

    // coeff_i3 = Sf_m C_mik3 g_k delta_33
    //
    //          = Sf_1 C_1ik3 g_k
    //            + Sf_2 C_2ik3 g_k
    //            + Sf_3 C_3ik3 g_k

    // coeff_13 = Sf_m C_m1k3 g_k
    //
    //          = Sf_1 C_11k3 g_k
    //            + Sf_2 C_21k3 g_k
    //            + Sf_3 C_31k3 g_k

    // coeff_13 = Sf_1 C_1113 g_1
    //            + Sf_1 C_1123 g_2
    //            + Sf_1 C_1133 g_3
    //
    //            + Sf_2 C_2113 g_1
    //            + Sf_2 C_2123 g_2
    //            + Sf_2 C_2133 g_3
    //
    //            + Sf_3 C_3113 g_1
    //            + Sf_3 C_3123 g_2
    //            + Sf_3 C_3133 g_3

    coeff[tensor::XZ] =
        (Sf[vector::X]*C(XX,tensor::XZ))*g[vector::X]
      + (Sf[vector::X]*C(XX,tensor::YZ))*g[vector::Y]
      + (Sf[vector::X]*C(XX,tensor::ZZ))*g[vector::Z]

      + (Sf[vector::Y]*C(XY,tensor::XZ))*g[vector::X]
      + (Sf[vector::Y]*C(XY,tensor::YZ))*g[vector::Y]
      + (Sf[vector::Y]*C(XY,tensor::ZZ))*g[vector::Z]

      + (Sf[vector::Z]*C(XZ,tensor::XZ))*g[vector::X]
      + (Sf[vector::Z]*C(XZ,tensor::YZ))*g[vector::Y]
      + (Sf[vector::Z]*C(XZ,tensor::ZZ))*g[vector::Z];

//    coeff[tensor::XZ] =
//        Sf[vector::X]*C(XX,ZX)*g[vector::X]
//      + Sf[vector::X]*C(XX,YZ)*g[vector::Y]
//      + Sf[vector::X]*C(XX,ZZ)*g[vector::Z]
//
//      + Sf[vector::Y]*C(XY,ZX)*g[vector::X]
//      + Sf[vector::Y]*C(XY,YZ)*g[vector::Y]
//      + Sf[vector::Y]*C(XY,ZZ)*g[vector::Z]
//
//      + Sf[vector::Z]*C(ZX,ZX)*g[vector::X]
//      + Sf[vector::Z]*C(ZX,YZ)*g[vector::Y]
//      + Sf[vector::Z]*C(ZX,ZZ)*g[vector::Z];

    // coeff_23 = Sf_m C_m2k3 g_k
    //
    //          = Sf_1 C_12k3 g_k
    //            + Sf_2 C_22k3 g_k
    //            + Sf_3 C_32k3 g_k

    //          = Sf_1 C_1213 g_1
    //            + Sf_1 C_1223 g_2
    //            + Sf_1 C_1233 g_3
    //
    //            + Sf_2 C_2213 g_1
    //            + Sf_2 C_2223 g_2
    //            + Sf_2 C_2233 g_3
    //
    //            + Sf_3 C_3213 g_1
    //            + Sf_3 C_3223 g_2
    //            + Sf_3 C_3233 g_3

    coeff[tensor::YZ] =
        (Sf[vector::X]*C(XY,tensor::XZ))*g[vector::X]
      + (Sf[vector::X]*C(XY,tensor::YZ))*g[vector::Y]
      + (Sf[vector::X]*C(XY,tensor::ZZ))*g[vector::Z]

      + (Sf[vector::Y]*C(YY,tensor::XZ))*g[vector::X]
      + (Sf[vector::Y]*C(YY,tensor::YZ))*g[vector::Y]
      + (Sf[vector::Y]*C(YY,tensor::ZZ))*g[vector::Z]

      + (Sf[vector::Z]*C(YZ,tensor::XZ))*g[vector::X]
      + (Sf[vector::Z]*C(YZ,tensor::YZ))*g[vector::Y]
      + (Sf[vector::Z]*C(YZ,tensor::ZZ))*g[vector::Z];

//    coeff[tensor::YZ] =
//        Sf[vector::X]*C(XY,ZX)*g[vector::X]
//      + Sf[vector::X]*C(XY,YZ)*g[vector::Y]
//      + Sf[vector::X]*C(XY,ZZ)*g[vector::Z]
//
//      + Sf[vector::Y]*C(YY,ZX)*g[vector::X]
//      + Sf[vector::Y]*C(YY,YZ)*g[vector::Y]
//      + Sf[vector::Y]*C(YY,ZZ)*g[vector::Z]
//
//      + Sf[vector::Z]*C(YZ,ZX)*g[vector::X]
//      + Sf[vector::Z]*C(YZ,YZ)*g[vector::Y]
//      + Sf[vector::Z]*C(YZ,ZZ)*g[vector::Z];

    // coeff_33 = Sf_m C_m3k3 g_k
    //
    //          = Sf_1 C_13k3 g_k
    //            + Sf_2 C_23k3 g_k
    //            + Sf_3 C_33k3 g_k

    //          = Sf_1 C_1313 g_1
    //            + Sf_1 C_1323 g_2
    //            + Sf_1 C_1333 g_3
    //
    //            + Sf_2 C_2313 g_1
    //            + Sf_2 C_2323 g_2
    //            + Sf_2 C_2333 g_3
    //
    //            + Sf_3 C_3313 g_1
    //            + Sf_3 C_3323 g_2
    //            + Sf_3 C_3333 g_3
    // Note: C_1313 == C_3131

    coeff[tensor::ZZ] =
        (Sf[vector::X]*C(XZ,tensor::XZ))*g[vector::X]
      + (Sf[vector::X]*C(XZ,tensor::YZ))*g[vector::Y]
      + (Sf[vector::X]*C(XZ,tensor::ZZ))*g[vector::Z]

      + (Sf[vector::Y]*C(YZ,tensor::XZ))*g[vector::X]
      + (Sf[vector::Y]*C(YZ,tensor::YZ))*g[vector::Y]
      + (Sf[vector::Y]*C(YZ,tensor::ZZ))*g[vector::Z]

      + (Sf[vector::Z]*C(ZZ,tensor::XZ))*g[vector::X]
      + (Sf[vector::Z]*C(ZZ,tensor::YZ))*g[vector::Y]
      + (Sf[vector::Z]*C(ZZ,tensor::ZZ))*g[vector::Z];

//    coeff[tensor::ZZ] =
//        Sf[vector::X]*C(ZX,ZX)*g[vector::X]
//      + Sf[vector::X]*C(ZX,YZ)*g[vector::Y]
//      + Sf[vector::X]*C(ZX,ZZ)*g[vector::Z]
//
//      + Sf[vector::Y]*C(YZ,ZX)*g[vector::X]
//      + Sf[vector::Y]*C(YZ,YZ)*g[vector::Y]
//      + Sf[vector::Y]*C(YZ,ZZ)*g[vector::Z]
//
//      + Sf[vector::Z]*C(ZZ,ZX)*g[vector::X]
//      + Sf[vector::Z]*C(ZZ,YZ)*g[vector::Y]
//      + Sf[vector::Z]*C(ZZ,ZZ)*g[vector::Z];

}

// ************************************************************************* //
