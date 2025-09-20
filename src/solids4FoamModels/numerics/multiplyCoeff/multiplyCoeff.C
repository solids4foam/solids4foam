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

#include "multiplyCoeff.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::multiplyCoeff
(
    tensor& coeff,
    const vector& Sf,
    const mat66& C,
    const vector& g
)
{
#ifdef FULLDEBUG
    // Check the dimensions of C are correct
    if (C.m() != 6 || C.n() != 6)
    {
        FatalError
            << "The material tangent matrix C has the wrong dimensions!" << nl
            << "It should be 6 x 6 but it is " << C.m() << " x " << C.n()
            << abort(FatalError);
    }
#endif

    // Define matrix indices for readability
    const label XX = symmTensor::XX;
    const label YY = symmTensor::YY;
    const label ZZ = symmTensor::ZZ;
    const label XY = symmTensor::XY;
    const label YZ = symmTensor::YZ;
    const label XZ = symmTensor::XZ;

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
        Sf[vector::X]*C(XX,XX)*g[vector::X]
      + Sf[vector::X]*C(XX,XY)*g[vector::Y]
      + Sf[vector::X]*C(XX,XZ)*g[vector::Z]

      + Sf[vector::Y]*C(XY,XX)*g[vector::X]
      + Sf[vector::Y]*C(XY,XY)*g[vector::Y]
      + Sf[vector::Y]*C(XY,XZ)*g[vector::Z]

      + Sf[vector::Z]*C(XZ,XX)*g[vector::X]
      + Sf[vector::Z]*C(XZ,XY)*g[vector::Y]
      + Sf[vector::Z]*C(XZ,XZ)*g[vector::Z];

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
        Sf[vector::X]*C(XY,XX)*g[vector::X]
      + Sf[vector::X]*C(XY,XY)*g[vector::Y]
      + Sf[vector::X]*C(XY,XZ)*g[vector::Z]

      + Sf[vector::Y]*C(YY,XX)*g[vector::X]
      + Sf[vector::Y]*C(YY,XY)*g[vector::Y]
      + Sf[vector::Y]*C(YY,XZ)*g[vector::Z]

      + Sf[vector::Z]*C(YZ,XX)*g[vector::X]
      + Sf[vector::Z]*C(YZ,XY)*g[vector::Y]
      + Sf[vector::Z]*C(YZ,XZ)*g[vector::Z];

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
        Sf[vector::X]*C(XZ,XX)*g[vector::X]
      + Sf[vector::X]*C(XZ,XY)*g[vector::Y]
      + Sf[vector::X]*C(XZ,XZ)*g[vector::Z]

      + Sf[vector::Y]*C(YZ,XX)*g[vector::X]
      + Sf[vector::Y]*C(YZ,XY)*g[vector::Y]
      + Sf[vector::Y]*C(YZ,XZ)*g[vector::Z]

      + Sf[vector::Z]*C(ZZ,XX)*g[vector::X]
      + Sf[vector::Z]*C(ZZ,XY)*g[vector::Y]
      + Sf[vector::Z]*C(ZZ,XZ)*g[vector::Z];

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
        Sf[vector::X]*C(XX,XY)*g[vector::X]
      + Sf[vector::X]*C(XX,YY)*g[vector::Y]
      + Sf[vector::X]*C(XX,YZ)*g[vector::Z]

      + Sf[vector::Y]*C(XY,XY)*g[vector::X]
      + Sf[vector::Y]*C(XY,YY)*g[vector::Y]
      + Sf[vector::Y]*C(XY,YZ)*g[vector::Z]

      + Sf[vector::Z]*C(XZ,XY)*g[vector::X]
      + Sf[vector::Z]*C(XZ,YY)*g[vector::Y]
      + Sf[vector::Z]*C(XZ,YZ)*g[vector::Z];

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
        Sf[vector::X]*C(XY,XY)*g[vector::X]
      + Sf[vector::X]*C(XY,YY)*g[vector::Y]
      + Sf[vector::X]*C(XY,YZ)*g[vector::Z]

      + Sf[vector::Y]*C(YY,XY)*g[vector::X]
      + Sf[vector::Y]*C(YY,YY)*g[vector::Y]
      + Sf[vector::Y]*C(YY,YZ)*g[vector::Z]

      + Sf[vector::Z]*C(YZ,XY)*g[vector::X]
      + Sf[vector::Z]*C(YZ,YY)*g[vector::Y]
      + Sf[vector::Z]*C(YZ,YZ)*g[vector::Z];

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
        Sf[vector::X]*C(XZ,XY)*g[vector::X]
      + Sf[vector::X]*C(XZ,YY)*g[vector::Y]
      + Sf[vector::X]*C(XZ,YZ)*g[vector::Z]

      + Sf[vector::Y]*C(YZ,XY)*g[vector::X]
      + Sf[vector::Y]*C(YZ,YY)*g[vector::Y]
      + Sf[vector::Y]*C(YZ,YZ)*g[vector::Z]

      + Sf[vector::Z]*C(ZZ,XY)*g[vector::X]
      + Sf[vector::Z]*C(ZZ,YY)*g[vector::Y]
      + Sf[vector::Z]*C(ZZ,YZ)*g[vector::Z];

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
        Sf[vector::X]*C(XX,XZ)*g[vector::X]
      + Sf[vector::X]*C(XX,YZ)*g[vector::Y]
      + Sf[vector::X]*C(XX,ZZ)*g[vector::Z]

      + Sf[vector::Y]*C(XY,XZ)*g[vector::X]
      + Sf[vector::Y]*C(XY,YZ)*g[vector::Y]
      + Sf[vector::Y]*C(XY,ZZ)*g[vector::Z]

      + Sf[vector::Z]*C(XZ,XZ)*g[vector::X]
      + Sf[vector::Z]*C(XZ,YZ)*g[vector::Y]
      + Sf[vector::Z]*C(XZ,ZZ)*g[vector::Z];

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
        Sf[vector::X]*C(XY,XZ)*g[vector::X]
      + Sf[vector::X]*C(XY,YZ)*g[vector::Y]
      + Sf[vector::X]*C(XY,ZZ)*g[vector::Z]

      + Sf[vector::Y]*C(YY,XZ)*g[vector::X]
      + Sf[vector::Y]*C(YY,YZ)*g[vector::Y]
      + Sf[vector::Y]*C(YY,ZZ)*g[vector::Z]

      + Sf[vector::Z]*C(YZ,XZ)*g[vector::X]
      + Sf[vector::Z]*C(YZ,YZ)*g[vector::Y]
      + Sf[vector::Z]*C(YZ,ZZ)*g[vector::Z];

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
        Sf[vector::X]*C(XZ,XZ)*g[vector::X]
      + Sf[vector::X]*C(XZ,YZ)*g[vector::Y]
      + Sf[vector::X]*C(XZ,ZZ)*g[vector::Z]

      + Sf[vector::Y]*C(YZ,XZ)*g[vector::X]
      + Sf[vector::Y]*C(YZ,YZ)*g[vector::Y]
      + Sf[vector::Y]*C(YZ,ZZ)*g[vector::Z]

      + Sf[vector::Z]*C(ZZ,XZ)*g[vector::X]
      + Sf[vector::Z]*C(ZZ,YZ)*g[vector::Y]
      + Sf[vector::Z]*C(ZZ,ZZ)*g[vector::Z];
}


void Foam::multiplyCoeff
(
    tensor& coeff,
    const vector& Sf,
    const mat66& C,
    const mat39& G,
    const symmTensor& sigma,
    const vector& g
)
{
#ifdef FULLDEBUG
    // Check the dimensions of C are correct
    if (C.m() != 6 || C.n() != 6)
    {
        FatalError
            << "The material tangent matrix C has the wrong dimensions!" << nl
            << "It should be 6 x 6 but it is " << C.m() << " x " << C.n()
            << abort(FatalError);
    }

    // Check the dimensions of G are correct
    if (G.m() != 3 || G.n() != 9)
    {
        FatalError
            << "The geometric stiffness G has the wrong dimensions!" << nl
            << "It should be 3 x 9 but it is " << G.m() << " x " << G.n()
            << abort(FatalError);
    }
#endif

    // Define matrix indices for readability (for C only)
    const label XX = symmTensor::XX;
    const label YY = symmTensor::YY;
    const label ZZ = symmTensor::ZZ;
    const label XY = symmTensor::XY, YX = XY;
    const label YZ = symmTensor::YZ, ZY = YZ;
    const label XZ = symmTensor::XZ, ZX = XZ;

    // Index notation
    // coeff_ij = (Sf_m C_mikl + G_mkl sigma_mi) g_k delta_lj
    // where delta is kronecker delta

    // Term is only non-zero for l = j from 1 to 3:
    // coeff_i1 = (Sf_m C_mik1 + G_mk1 sigma_mi) g_k delta_11
    // coeff_i2 = (Sf_m C_mik2 + G_mk2 sigma_mi) g_k delta_22
    // coeff_i3 = (Sf_m C_mik3 + G_mk3 sigma_mi) g_k delta_33

    // coeff_11 = (Sf_m C_m1k1 + G_mk1 sigma_m1)*g_k
    //
    //          = (Sf_1 C_11k1 + G_1k1 sigma_11) g_k
    //          + (Sf_2 C_21k1 + G_2k1 sigma_21) g_k
    //          + (Sf_3 C_31k1 + G_3k1 sigma_31) g_k
    //
    //          = (Sf_1 C_1111 + G_111 sigma_11) g_1
    //          + (Sf_1 C_1121 + G_121 sigma_11) g_2
    //          + (Sf_1 C_1131 + G_131 sigma_11) g_3

    //          + (Sf_2 C_2111 + G_211 sigma_21) g_1
    //          + (Sf_2 C_2121 + G_221 sigma_21) g_2
    //          + (Sf_2 C_2131 + G_231 sigma_21) g_3

    //          + (Sf_3 C_3111 + G_311 sigma_31) g_1
    //          + (Sf_3 C_3121 + G_321 sigma_31) g_2
    //          + (Sf_3 C_3131 + G_331 sigma_31) g_3

    coeff[tensor::XX] =
        (Sf[vector::X]*C(XX, XX)
      + G(vector::X,tensor::XX)*sigma[symmTensor::XX])*g[vector::X]
      + (Sf[vector::X]*C(XX, YX)
      + G(vector::X,tensor::YX)*sigma[symmTensor::XX])*g[vector::Y]
      + (Sf[vector::X]*C(XX, ZX)
      + G(vector::X,tensor::ZX)*sigma[symmTensor::XX])*g[vector::Z]

      + (Sf[vector::Y]*C(XY, XX)
      + G(vector::Y,tensor::XX)*sigma[symmTensor::XY])*g[vector::X]
      + (Sf[vector::Y]*C(XY, YX)
      + G(vector::Y,tensor::YX)*sigma[symmTensor::XY])*g[vector::Y]
      + (Sf[vector::Y]*C(XY, ZX)
      + G(vector::Y,tensor::ZX)*sigma[symmTensor::XY])*g[vector::Z]

      + (Sf[vector::Z]*C(XZ, XX)
      + G(vector::Z,tensor::XX)*sigma[symmTensor::XZ])*g[vector::X]
      + (Sf[vector::Z]*C(XZ, YX)
      + G(vector::Z,tensor::YX)*sigma[symmTensor::XZ])*g[vector::Y]
      + (Sf[vector::Z]*C(XZ, ZX)
      + G(vector::Z,tensor::ZX)*sigma[symmTensor::XZ])*g[vector::Z];

    // Similarly the other components can be calculated as
    // coeff_21 = (Sf_m C_m2k1 + G_mk1 sigma_m2) g_k
    //
    //          = (Sf_1 C_12k1 + G_1k1 sigma_12) g_k
    //          + (Sf_2 C_22k1 + G_2k1 sigma_22) g_k
    //          + (Sf_3 C_32k1 + G_3k1 sigma_32) g_k

    //          = (Sf_1 C_1211 + G_111 sigma_12) g_1
    //          + (Sf_1 C_1221 + G_121 sigma_12) g_2
    //          + (Sf_1 C_1231 + G_131 sigma_12) g_3

    //          + (Sf_2 C_2211 + G_211 sigma_22) g_1
    //          + (Sf_2 C_2221 + G_221 sigma_22) g_2
    //          + (Sf_2 C_2231 + G_231 sigma_22) g_3

    //          + (Sf_3 C_3211 + G_311 sigma_32) g_1
    //          + (Sf_3 C_3221 + G_321 sigma_32) g_2
    //          + (Sf_3 C_3231 + G_331 sigma_32) g_3

    coeff[tensor::YX] =
        (Sf[vector::X]*C(XY, XX)
      + G(vector::X,tensor::XX)*sigma[symmTensor::XY])*g[vector::X]
      + (Sf[vector::X]*C(XY, YX)
      + G(vector::X,tensor::YX)*sigma[symmTensor::XY])*g[vector::Y]
      + (Sf[vector::X]*C(XY, ZX)
      + G(vector::X,tensor::ZX)*sigma[symmTensor::XY])*g[vector::Z]

      + (Sf[vector::Y]*C(YY, XX)
      + G(vector::Y,tensor::XX)*sigma[symmTensor::YY])*g[vector::X]
      + (Sf[vector::Y]*C(YY, YX)
      + G(vector::Y,tensor::YX)*sigma[symmTensor::YY])*g[vector::Y]
      + (Sf[vector::Y]*C(YY, ZX)
      + G(vector::Y,tensor::ZX)*sigma[symmTensor::YY])*g[vector::Z]

      + (Sf[vector::Z]*C(YZ, XX)
      + G(vector::Z,tensor::XX)*sigma[symmTensor::YZ])*g[vector::X]
      + (Sf[vector::Z]*C(YZ, YX)
      + G(vector::Z,tensor::YX)*sigma[symmTensor::YZ])*g[vector::Y]
      + (Sf[vector::Z]*C(YZ, ZX)
      + G(vector::Z,tensor::ZX)*sigma[symmTensor::YZ])*g[vector::Z];

    // coeff_31 = (Sf_m C_m3k1 + G_mk1 sigma_m3) g_k
    //
    //          = (Sf_1 C_13k1 + G_1k1 sigma_13) g_k
    //          + (Sf_2 C_23k1 + G_2k1 sigma_23) g_k
    //          + (Sf_3 C_33k1 + G_3k1 sigma_33) g_k

    //          = (Sf_1 C_1311 + G_111 sigma_13) g_1
    //          + (Sf_1 C_1321 + G_121 sigma_13) g_2
    //          + (Sf_1 C_1331 + G_131 sigma_13) g_3

    //          + (Sf_2 C_2311 + G_211 sigma_23) g_1
    //          + (Sf_2 C_2321 + G_221 sigma_23) g_2
    //          + (Sf_2 C_2331 + G_231 sigma_23) g_3

    //          + (Sf_3 C_3311 + G_311 sigma_33) g_1
    //          + (Sf_3 C_3321 + G_321 sigma_33) g_2
    //          + (Sf_3 C_3331 + G_331 sigma_33) g_3

    coeff[tensor::ZX] =
        (Sf[vector::X]*C(XZ, XX)
      + G(vector::X,tensor::XX)*sigma[symmTensor::XZ])*g[vector::X]
      + (Sf[vector::X]*C(XZ, YX)
      + G(vector::X,tensor::YX)*sigma[symmTensor::XZ])*g[vector::Y]
      + (Sf[vector::X]*C(XZ, ZX)
      + G(vector::X,tensor::ZX)*sigma[symmTensor::XZ])*g[vector::Z]

      + (Sf[vector::Y]*C(YZ, XX)
      + G(vector::Y,tensor::XX)*sigma[symmTensor::YZ])*g[vector::X]
      + (Sf[vector::Y]*C(YZ, YX)
      + G(vector::Y,tensor::YX)*sigma[symmTensor::YZ])*g[vector::Y]
      + (Sf[vector::Y]*C(YZ, ZX)
      + G(vector::Y,tensor::ZX)*sigma[symmTensor::YZ])*g[vector::Z]

      + (Sf[vector::Z]*C(ZZ, XX)
      + G(vector::Z,tensor::XX)*sigma[symmTensor::ZZ])*g[vector::X]
      + (Sf[vector::Z]*C(ZZ, YX)
      + G(vector::Z,tensor::YX)*sigma[symmTensor::ZZ])*g[vector::Y]
      + (Sf[vector::Z]*C(ZZ, ZX)
      + G(vector::Z,tensor::ZX)*sigma[symmTensor::ZZ])*g[vector::Z];

    // coeff_i2 = (Sf_m C_mik2 + G_mk2 sigma_mi) g_k delta_22
    //
    //          = (Sf_1 C_1ik2 + G_1k2 sigma_1i) g_k
    //          + (Sf_2 C_2ik2 + G_2k2 sigma_2i)  g_k
    //          + (Sf_3 C_3ik2 + G_3k2 sigma_3i) g_k
    //
    //          = (Sf_1 C_1i12 + G_112 sigma_1i) g_1
    //          + (Sf_1 C_1i22 + G_122 sigma_1i) g_2
    //          + (Sf_1 C_1i32 + G_132 sigma_1i) g_3
    //
    //          + (Sf_2 C_2i12 + G_212 sigma_2i) g_1
    //          + (Sf_2 C_2i22 + G_222 sigma_2i) g_2
    //          + (Sf_2 C_2i32 + G_232 sigma_2i) g_3
    //
    //          + (Sf_3 C_3i12 + G_312 sigma_3i) g_1
    //          + (Sf_3 C_3i22 + G_322 sigma_3i) g_2
    //          + (Sf_3 C_3i32 + G_332 sigma_3i) g_3

    // coeff_12 = (Sf_m C_m1k2 + G_mk2 sigma_m1) g_k
    //
    //          = (Sf_1 C_11k2 + G_1k2 sigma_11) g_k
    //          + (Sf_2 C_21k2 + G_2k2 sigma_21)  g_k
    //          + (Sf_3 C_31k2 + G_3k2 sigma_31) g_k
    //
    //          = (Sf_1 C_1112 + G_112 sigma_11) g_1
    //          + (Sf_1 C_1122 + G_122 sigma_11) g_2
    //          + (Sf_1 C_1132 + G_132 sigma_11) g_3
    //
    //          + (Sf_2 C_2112 + G_212 sigma_21) g_1
    //          + (Sf_2 C_2122 + G_222 sigma_21) g_2
    //          + (Sf_2 C_2132 + G_232 sigma_21) g_3
    //
    //          + (Sf_3 C_3112 + G_312 sigma_31) g_1
    //          + (Sf_3 C_3122 + G_322 sigma_31) g_2
    //          + (Sf_3 C_3132 + G_332 sigma_31) g_3

    coeff[tensor::XY] =
        (Sf[vector::X]*C(XX, XY)
      + G(vector::X,tensor::XY)*sigma[symmTensor::XX])*g[vector::X]
      + (Sf[vector::X]*C(XX, YY)
      + G(vector::X,tensor::YY)*sigma[symmTensor::XX])*g[vector::Y]
      + (Sf[vector::X]*C(XX, ZY)
      + G(vector::X,tensor::ZY)*sigma[symmTensor::XX])*g[vector::Z]

      + (Sf[vector::Y]*C(XY, XY)
      + G(vector::Y,tensor::XY)*sigma[symmTensor::XY])*g[vector::X]
      + (Sf[vector::Y]*C(XY, YY)
      + G(vector::Y,tensor::YY)*sigma[symmTensor::XY])*g[vector::Y]
      + (Sf[vector::Y]*C(XY, ZY)
      + G(vector::Y,tensor::ZY)*sigma[symmTensor::XY])*g[vector::Z]

      + (Sf[vector::Z]*C(XZ, XY)
      + G(vector::Z,tensor::XY)*sigma[symmTensor::XZ])*g[vector::X]
      + (Sf[vector::Z]*C(XZ, YY)
      + G(vector::Z,tensor::YY)*sigma[symmTensor::XZ])*g[vector::Y]
      + (Sf[vector::Z]*C(XZ, ZY)
      + G(vector::Z,tensor::ZY)*sigma[symmTensor::XZ])*g[vector::Z];

    // coeff_22 = (Sf_m C_m2k2 + G_mk2 sigma_m2) g_k
    //
    //          = (Sf_1 C_12k2 + G_1k2 sigma_12) g_k
    //          + (Sf_2 C_22k2 + G_2k2 sigma_22)  g_k
    //          + (Sf_3 C_32k2 + G_3k2 sigma_32) g_k
    //
    //          = (Sf_1 C_1212 + G_112 sigma_12) g_1
    //          + (Sf_1 C_1222 + G_122 sigma_12) g_2
    //          + (Sf_1 C_1232 + G_132 sigma_12) g_3
    //
    //          + (Sf_2 C_2212 + G_212 sigma_22) g_1
    //          + (Sf_2 C_2222 + G_222 sigma_22) g_2
    //          + (Sf_2 C_2232 + G_232 sigma_22) g_3
    //
    //          + (Sf_3 C_3212 + G_312 sigma_32) g_1
    //          + (Sf_3 C_3222 + G_322 sigma_32) g_2
    //          + (Sf_3 C_3232 + G_332 sigma_32) g_3

    coeff[tensor::YY] =
        (Sf[vector::X]*C(XY, XY)
      + G(vector::X,tensor::XY)*sigma[symmTensor::XY])*g[vector::X]
      + (Sf[vector::X]*C(XY, YY)
      + G(vector::X,tensor::YY)*sigma[symmTensor::XY])*g[vector::Y]
      + (Sf[vector::X]*C(XY, ZY)
      + G(vector::X,tensor::ZY)*sigma[symmTensor::XY])*g[vector::Z]

      + (Sf[vector::Y]*C(YY, XY)
      + G(vector::Y,tensor::XY)*sigma[symmTensor::YY])*g[vector::X]
      + (Sf[vector::Y]*C(YY, YY)
      + G(vector::Y,tensor::YY)*sigma[symmTensor::YY])*g[vector::Y]
      + (Sf[vector::Y]*C(YY, ZY)
      + G(vector::Y,tensor::ZY)*sigma[symmTensor::YY])*g[vector::Z]

      + (Sf[vector::Z]*C(YZ, XY)
      + G(vector::Z,tensor::XY)*sigma[symmTensor::YZ])*g[vector::X]
      + (Sf[vector::Z]*C(YZ, YY)
      + G(vector::Z,tensor::YY)*sigma[symmTensor::YZ])*g[vector::Y]
      + (Sf[vector::Z]*C(YZ, ZY)
      + G(vector::Z,tensor::ZY)*sigma[symmTensor::YZ])*g[vector::Z];

    // coeff_32 = (Sf_m C_m3k2 + G_mk2 sigma_m3) g_k
    //
    //          = (Sf_1 C_13k2 + G_1k2 sigma_13) g_k
    //          + (Sf_2 C_23k2 + G_2k2 sigma_23)  g_k
    //          + (Sf_3 C_33k2 + G_3k2 sigma_33) g_k
    //
    //          = (Sf_1 C_1312 + G_112 sigma_13) g_1
    //          + (Sf_1 C_1322 + G_122 sigma_13) g_2
    //          + (Sf_1 C_1332 + G_132 sigma_13) g_3
    //
    //          + (Sf_2 C_2312 + G_212 sigma_23) g_1
    //          + (Sf_2 C_2322 + G_222 sigma_23) g_2
    //          + (Sf_2 C_2332 + G_232 sigma_23) g_3
    //
    //          + (Sf_3 C_3312 + G_312 sigma_33) g_1
    //          + (Sf_3 C_3322 + G_322 sigma_33) g_2
    //          + (Sf_3 C_3332 + G_332 sigma_33) g_3

    coeff[tensor::ZY] =
        (Sf[vector::X]*C(XZ, XY)
      + G(vector::X,tensor::XY)*sigma[symmTensor::XZ])*g[vector::X]
      + (Sf[vector::X]*C(XZ, YY)
      + G(vector::X,tensor::YY)*sigma[symmTensor::XZ])*g[vector::Y]
      + (Sf[vector::X]*C(XZ, ZY)
      + G(vector::X,tensor::ZY)*sigma[symmTensor::XZ])*g[vector::Z]

      + (Sf[vector::Y]*C(YZ, XY)
      + G(vector::Y,tensor::XY)*sigma[symmTensor::YZ])*g[vector::X]
      + (Sf[vector::Y]*C(YZ, YY)
      + G(vector::Y,tensor::YY)*sigma[symmTensor::YZ])*g[vector::Y]
      + (Sf[vector::Y]*C(YZ, ZY)
      + G(vector::Y,tensor::ZY)*sigma[symmTensor::YZ])*g[vector::Z]

      + (Sf[vector::Z]*C(ZZ, XY)
      + G(vector::Z,tensor::XY)*sigma[symmTensor::ZZ])*g[vector::X]
      + (Sf[vector::Z]*C(ZZ, YY)
      + G(vector::Z,tensor::YY)*sigma[symmTensor::ZZ])*g[vector::Y]
      + (Sf[vector::Z]*C(ZZ, ZY)
      + G(vector::Z,tensor::ZY)*sigma[symmTensor::ZZ])*g[vector::Z];

    // coeff_i3 = (Sf_m C_mik3 + G_mk3 sigma_mi) g_k delta_33
    //
    //          = (Sf_1 C_1ik3 + G_1k3 sigma_1i) g_k
    //          + (Sf_2 C_2ik3 + G_2k3 sigma_2i) g_k
    //          + (Sf_3 C_3ik3 + G_3k3 sigma_3i) g_k
    //
    //          = (Sf_1 C_1i13 + G_113 sigma_1i) g_1
    //          + (Sf_1 C_1i23 + G_123 sigma_1i) g_2
    //          + (Sf_1 C_1i33 + G_133 sigma_1i) g_3
    //
    //          + (Sf_2 C_2i13 + G_213 sigma_2i) g_1
    //          + (Sf_2 C_2i23 + G_223 sigma_2i) g_2
    //          + (Sf_2 C_2i33 + G_233 sigma_2i) g_3
    //
    //          + (Sf_3 C_3i13 + G_313 sigma_3i) g_1
    //          + (Sf_3 C_3i23 + G_323 sigma_3i) g_2
    //          + (Sf_3 C_3i33 + G_333 sigma_3i) g_3

    // coeff_13 = (Sf_m C_m1k3 + G_mk3 sigma_m1) g_k

    //          = (Sf_1 C_11k3 + G_1k3 sigma_11) g_k
    //          + (Sf_2 C_21k3 + G_2k3 sigma_21) g_k
    //          + (Sf_3 C_31k3 + G_3k3 sigma_31) g_k
    //
    //          = (Sf_1 C_1113 + G_113 sigma_11) g_1
    //          + (Sf_1 C_1123 + G_123 sigma_11) g_2
    //          + (Sf_1 C_1133 + G_133 sigma_11) g_3
    //
    //          + (Sf_2 C_2113 + G_213 sigma_21) g_1
    //          + (Sf_2 C_2123 + G_223 sigma_21) g_2
    //          + (Sf_2 C_2133 + G_233 sigma_21) g_3
    //
    //          + (Sf_3 C_3113 + G_313 sigma_31) g_1
    //          + (Sf_3 C_3123 + G_323 sigma_31) g_2
    //          + (Sf_3 C_3133 + G_333 sigma_31) g_3

    coeff[tensor::XZ] =
        (Sf[vector::X]*C(XX, XZ)
      + G(vector::X,tensor::XZ)*sigma[symmTensor::XX])*g[vector::X]
      + (Sf[vector::X]*C(XX, YZ)
      + G(vector::X,tensor::YZ)*sigma[symmTensor::XX])*g[vector::Y]
      + (Sf[vector::X]*C(XX, ZZ)
      + G(vector::X,tensor::ZZ)*sigma[symmTensor::XX])*g[vector::Z]

      + (Sf[vector::Y]*C(XY, XZ)
      + G(vector::Y,tensor::XZ)*sigma[symmTensor::XY])*g[vector::X]
      + (Sf[vector::Y]*C(XY, YZ)
      + G(vector::Y,tensor::YZ)*sigma[symmTensor::XY])*g[vector::Y]
      + (Sf[vector::Y]*C(XY, ZZ)
      + G(vector::Y,tensor::ZZ)*sigma[symmTensor::XY])*g[vector::Z]

      + (Sf[vector::Z]*C(XZ, XZ)
      + G(vector::Z,tensor::XZ)*sigma[symmTensor::XZ])*g[vector::X]
      + (Sf[vector::Z]*C(XZ, YZ)
      + G(vector::Z,tensor::YZ)*sigma[symmTensor::XZ])*g[vector::Y]
      + (Sf[vector::Z]*C(XZ, ZZ)
      + G(vector::Z,tensor::ZZ)*sigma[symmTensor::XZ])*g[vector::Z];

    // coeff_23 = (Sf_m C_m2k3 + G_mk3 sigma_m2) g_k
    //
    //          = (Sf_1 C_12k3 + G_1k3 sigma_12) g_k
    //          + (Sf_2 C_22k3 + G_2k3 sigma_22) g_k
    //          + (Sf_3 C_32k3 + G_3k3 sigma_32) g_k
    //
    //          = (Sf_1 C_1213 + G_113 sigma_12) g_1
    //          + (Sf_1 C_1223 + G_123 sigma_12) g_2
    //          + (Sf_1 C_1233 + G_133 sigma_12) g_3
    //
    //          + (Sf_2 C_2213 + G_213 sigma_22) g_1
    //          + (Sf_2 C_2223 + G_223 sigma_22) g_2
    //          + (Sf_2 C_2233 + G_233 sigma_22) g_3
    //
    //          + (Sf_3 C_3213 + G_313 sigma_32) g_1
    //          + (Sf_3 C_3223 + G_323 sigma_32) g_2
    //          + (Sf_3 C_3233 + G_333 sigma_32) g_3

    coeff[tensor::YZ] =
        (Sf[vector::X]*C(XY, XZ)
      + G(vector::X,tensor::XZ)*sigma[symmTensor::XY])*g[vector::X]
      + (Sf[vector::X]*C(XY, YZ)
      + G(vector::X,tensor::YZ)*sigma[symmTensor::XY])*g[vector::Y]
      + (Sf[vector::X]*C(XY, ZZ)
      + G(vector::X,tensor::ZZ)*sigma[symmTensor::XY])*g[vector::Z]

      + (Sf[vector::Y]*C(YY, XZ)
      + G(vector::Y,tensor::XZ)*sigma[symmTensor::YY])*g[vector::X]
      + (Sf[vector::Y]*C(YY, YZ)
      + G(vector::Y,tensor::YZ)*sigma[symmTensor::YY])*g[vector::Y]
      + (Sf[vector::Y]*C(YY, ZZ)
      + G(vector::Y,tensor::ZZ)*sigma[symmTensor::YY])*g[vector::Z]

      + (Sf[vector::Z]*C(YZ, XZ)
      + G(vector::Z,tensor::XZ)*sigma[symmTensor::YZ])*g[vector::X]
      + (Sf[vector::Z]*C(YZ, YZ)
      + G(vector::Z,tensor::YZ)*sigma[symmTensor::YZ])*g[vector::Y]
      + (Sf[vector::Z]*C(YZ, ZZ)
      + G(vector::Z,tensor::ZZ)*sigma[symmTensor::YZ])*g[vector::Z];

    // coeff_33 = (Sf_m C_m3k3 + G_mk3 sigma_m3) g_k
    //
    //          = (Sf_1 C_13k3 + G_1k3 sigma_13) g_k
    //          + (Sf_2 C_23k3 + G_2k3 sigma_23) g_k
    //          + (Sf_3 C_33k3 + G_3k3 sigma_33) g_k
    //
    //          = (Sf_1 C_1313 + G_113 sigma_13) g_1
    //          + (Sf_1 C_1323 + G_123 sigma_13) g_2
    //          + (Sf_1 C_1333 + G_133 sigma_13) g_3
    //
    //          + (Sf_2 C_2313 + G_213 sigma_23) g_1
    //          + (Sf_2 C_2323 + G_223 sigma_23) g_2
    //          + (Sf_2 C_2333 + G_233 sigma_23) g_3
    //
    //          + (Sf_3 C_3313 + G_313 sigma_33) g_1
    //          + (Sf_3 C_3323 + G_323 sigma_33) g_2
    //          + (Sf_3 C_3333 + G_333 sigma_33) g_3

    coeff[tensor::ZZ] =
        (Sf[vector::X]*C(XZ, XZ)
      + G(vector::X,tensor::XZ)*sigma[symmTensor::XZ])*g[vector::X]
      + (Sf[vector::X]*C(XZ, YZ)
      + G(vector::X,tensor::YZ)*sigma[symmTensor::XZ])*g[vector::Y]
      + (Sf[vector::X]*C(XZ, ZZ)
      + G(vector::X,tensor::ZZ)*sigma[symmTensor::XZ])*g[vector::Z]

      + (Sf[vector::Y]*C(YZ, XZ)
      + G(vector::Y,tensor::XZ)*sigma[symmTensor::YZ])*g[vector::X]
      + (Sf[vector::Y]*C(YZ, YZ)
      + G(vector::Y,tensor::YZ)*sigma[symmTensor::YZ])*g[vector::Y]
      + (Sf[vector::Y]*C(YZ, ZZ)
      + G(vector::Y,tensor::ZZ)*sigma[symmTensor::YZ])*g[vector::Z]

      + (Sf[vector::Z]*C(ZZ, XZ)
      + G(vector::Z,tensor::XZ)*sigma[symmTensor::ZZ])*g[vector::X]
      + (Sf[vector::Z]*C(ZZ, YZ)
      + G(vector::Z,tensor::YZ)*sigma[symmTensor::ZZ])*g[vector::Y]
      + (Sf[vector::Z]*C(ZZ, ZZ)
      + G(vector::Z,tensor::ZZ)*sigma[symmTensor::ZZ])*g[vector::Z];
}


void Foam::multiplyCoeff
(
    vector& coeff,
    const mat39& G,
    const vector& gradP,
    const vector& g
)
{
#ifdef FULLDEBUG
    // Check the dimensions of G are correct
    if (G.m() != 3 || G.n() != 9)
    {
        FatalError
            << "The geometric stiffness G has the wrong dimensions!" << nl
            << "It should be 3 x 9 but it is " << G.m() << " x " << G.n()
            << abort(FatalError);
    }
#endif

    // Index notation
    // coeff_l = gradP_m G_mkl g_k

    // coeff_1 = gradP_m G_mk1 g_k
    //
    //          = gradP_1 G_1k1 g_k
    //          + gradP_2 G_2k1 g_k
    //          + gradP_3 G_3k1 g_k
    //
    //          = gradP_1 G_111 g_1
    //          + gradP_1 G_121 g_2
    //          + gradP_1 G_131 g_3
    //
    //          + gradP_2 G_211 g_1
    //          + gradP_2 G_221 g_2
    //          + gradP_2 G_231 g_3
    //
    //          + gradP_3 G_311 g_1
    //          + gradP_3 G_321 g_2
    //          + gradP_3 G_331 g_3

    coeff[vector::X] =
        gradP[vector::X]*G(vector::X, tensor::XX)*g[vector::X]
      + gradP[vector::X]*G(vector::X, tensor::YX)*g[vector::Y]
      + gradP[vector::X]*G(vector::X, tensor::ZX)*g[vector::Z]

      + gradP[vector::Y]*G(vector::Y, tensor::XX)*g[vector::X]
      + gradP[vector::Y]*G(vector::Y, tensor::YX)*g[vector::Y]
      + gradP[vector::Y]*G(vector::Y, tensor::ZX)*g[vector::Z]

      + gradP[vector::Z]*G(vector::Z, tensor::XX)*g[vector::X]
      + gradP[vector::Z]*G(vector::Z, tensor::YX)*g[vector::Y]
      + gradP[vector::Z]*G(vector::Z, tensor::ZX)*g[vector::Z];

    // coeff_2 = gradP_m G_mk2 g_k
    //
    //          = gradP_1 G_1k2 g_k
    //          + gradP_2 G_2k2 g_k
    //          + gradP_3 G_3k2 g_k
    //
    //          = gradP_1 G_112 g_1
    //          + gradP_1 G_122 g_2
    //          + gradP_1 G_132 g_3
    //
    //          + gradP_2 G_212 g_1
    //          + gradP_2 G_222 g_2
    //          + gradP_2 G_232 g_3
    //
    //          + gradP_3 G_312 g_1
    //          + gradP_3 G_322 g_2
    //          + gradP_3 G_332 g_3

    coeff[vector::Y] =
        gradP[vector::X]*G(vector::X, tensor::XY)*g[vector::X]
      + gradP[vector::X]*G(vector::X, tensor::YY)*g[vector::Y]
      + gradP[vector::X]*G(vector::X, tensor::ZY)*g[vector::Z]

      + gradP[vector::Y]*G(vector::Y, tensor::XY)*g[vector::X]
      + gradP[vector::Y]*G(vector::Y, tensor::YY)*g[vector::Y]
      + gradP[vector::Y]*G(vector::Y, tensor::ZY)*g[vector::Z]

      + gradP[vector::Z]*G(vector::Z, tensor::XY)*g[vector::X]
      + gradP[vector::Z]*G(vector::Z, tensor::YY)*g[vector::Y]
      + gradP[vector::Z]*G(vector::Z, tensor::ZY)*g[vector::Z];

    // coeff_1 = gradP_m G_mk3 g_k
    //
    //          = gradP_1 G_1k3 g_k
    //          + gradP_2 G_2k3 g_k
    //          + gradP_3 G_3k3 g_k
    //
    //          = gradP_1 G_113 g_1
    //          + gradP_1 G_123 g_2
    //          + gradP_1 G_133 g_3
    //
    //          + gradP_2 G_213 g_1
    //          + gradP_2 G_223 g_2
    //          + gradP_2 G_233 g_3
    //
    //          + gradP_3 G_313 g_1
    //          + gradP_3 G_323 g_2
    //          + gradP_3 G_333 g_3

    coeff[vector::Z] =
        gradP[vector::X]*G(vector::X, tensor::XZ)*g[vector::X]
      + gradP[vector::X]*G(vector::X, tensor::YZ)*g[vector::Y]
      + gradP[vector::X]*G(vector::X, tensor::ZZ)*g[vector::Z]

      + gradP[vector::Y]*G(vector::Y, tensor::XZ)*g[vector::X]
      + gradP[vector::Y]*G(vector::Y, tensor::YZ)*g[vector::Y]
      + gradP[vector::Y]*G(vector::Y, tensor::ZZ)*g[vector::Z]

      + gradP[vector::Z]*G(vector::Z, tensor::XZ)*g[vector::X]
      + gradP[vector::Z]*G(vector::Z, tensor::YZ)*g[vector::Y]
      + gradP[vector::Z]*G(vector::Z, tensor::ZZ)*g[vector::Z];
}


// ************************************************************************* //
