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

#include "multiplyCoeffExtended.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::multiplyCoeffExtended
(
    tensor& coeff,
    const vector& Sf,
    const RectangularMatrix<scalar>& C,
    const RectangularMatrix<scalar>& G,
    const symmTensor& sigma,
    const vector& g
)
{
    //Define indexes for readability only for C
    const label XX = 0;
    const label YY = 1;
    const label ZZ = 2;
    const label XY = 3;
    const label YZ = 4;
    const label XZ = 5;

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
    //            + (Sf_2 C_21k1 + G_2k1 sigma_21) g_k
    //            + (Sf_3 C_31k1 + G_3k1 sigma_31) g_k
    //
    //          = (Sf_1 C_1111 + G_111 sigma_11) g_1
    //            + (Sf_1 C_1121 + G_121 sigma_11) g_2
    //            + (Sf_1 C_1131 + G_131 sigma_11) g_3

    //            + (Sf_2 C_2111 + G_211 sigma_21) g_1
    //            + (Sf_2 C_2121 + G_221 sigma_21) g_2
    //            + (Sf_2 C_2131 + G_231 sigma_21) g_3

    //            + (Sf_3 C_3111 + G_311 sigma_31) g_1
    //            + (Sf_3 C_3121 + G_321 sigma_31) g_2
    //            + (Sf_3 C_3131 + G_331 sigma_31) g_3

    coeff[tensor::XX] =
        (Sf[vector::X]*C(XX,tensor::XX) + G(vector::X,tensor::XX)*sigma[symmTensor::XX])*g[vector::X]
      + (Sf[vector::X]*C(XX,tensor::YX) + G(vector::X,tensor::YX)*sigma[symmTensor::XX])*g[vector::Y]
      + (Sf[vector::X]*C(XX,tensor::ZX) + G(vector::X,tensor::ZX)*sigma[symmTensor::XX])*g[vector::Z]

      + (Sf[vector::Y]*C(XY,tensor::XX) + G(vector::Y,tensor::XX)*sigma[symmTensor::XY])*g[vector::X]
      + (Sf[vector::Y]*C(XY,tensor::YX) + G(vector::Y,tensor::YX)*sigma[symmTensor::XY])*g[vector::Y]
      + (Sf[vector::Y]*C(XY,tensor::ZX) + G(vector::Y,tensor::ZX)*sigma[symmTensor::XY])*g[vector::Z]

      + (Sf[vector::Z]*C(XZ,tensor::XX) + G(vector::Z,tensor::XX)*sigma[symmTensor::XZ])*g[vector::X]
      + (Sf[vector::Z]*C(XZ,tensor::YX) + G(vector::Z,tensor::YX)*sigma[symmTensor::XZ])*g[vector::Y]
      + (Sf[vector::Z]*C(XZ,tensor::ZX) + G(vector::Z,tensor::ZX)*sigma[symmTensor::XZ])*g[vector::Z];

    // Similarly the other components can be calculated as
    // coeff_21 = (Sf_m C_m2k1 + G_mk1 sigma_m2) g_k
    //
    //          = (Sf_1 C_12k1 + G_1k1 sigma_12) g_k
    //            + (Sf_2 C_22k1 + G_2k1 sigma_22) g_k
    //            + (Sf_3 C_32k1 + G_3k1 sigma_32) g_k

    //          = (Sf_1 C_1211 + G_111 sigma_12) g_1
    //            + (Sf_1 C_1221 + G_121 sigma_12) g_2
    //            + (Sf_1 C_1231 + G_131 sigma_12) g_3

    //            + (Sf_2 C_2211 + G_211 sigma_22) g_1
    //            + (Sf_2 C_2221 + G_221 sigma_22) g_2
    //            + (Sf_2 C_2231 + G_231 sigma_22) g_3

    //            + (Sf_3 C_3211 + G_311 sigma_32) g_1
    //            + (Sf_3 C_3221 + G_321 sigma_32) g_2
    //            + (Sf_3 C_3231 + G_331 sigma_32) g_3

    coeff[tensor::YX] =
        (Sf[vector::X]*C(XY,tensor::XX) + G(vector::X,tensor::XX)*sigma[symmTensor::XY])*g[vector::X]
      + (Sf[vector::X]*C(XY,tensor::YX) + G(vector::X,tensor::YX)*sigma[symmTensor::XY])*g[vector::Y]
      + (Sf[vector::X]*C(XY,tensor::ZX) + G(vector::X,tensor::ZX)*sigma[symmTensor::XY])*g[vector::Z]

      + (Sf[vector::Y]*C(YY,tensor::XX) + G(vector::Y,tensor::XX)*sigma[symmTensor::YY])*g[vector::X]
      + (Sf[vector::Y]*C(YY,tensor::YX) + G(vector::Y,tensor::YX)*sigma[symmTensor::YY])*g[vector::Y]
      + (Sf[vector::Y]*C(YY,tensor::ZX) + G(vector::Y,tensor::ZX)*sigma[symmTensor::YY])*g[vector::Z]

      + (Sf[vector::Z]*C(YZ,tensor::XX) + G(vector::Z,tensor::XX)*sigma[symmTensor::YZ])*g[vector::X]
      + (Sf[vector::Z]*C(YZ,tensor::YX) + G(vector::Z,tensor::YX)*sigma[symmTensor::YZ])*g[vector::Y]
      + (Sf[vector::Z]*C(YZ,tensor::ZX) + G(vector::Z,tensor::ZX)*sigma[symmTensor::YZ])*g[vector::Z];

    // coeff_31 = (Sf_m C_m3k1 + G_mk1 sigma_m3) g_k
    //
    //          = (Sf_1 C_13k1 + G_1k1 sigma_13) g_k
    //            + (Sf_2 C_23k1 + G_2k1 sigma_23) g_k
    //            + (Sf_3 C_33k1 + G_3k1 sigma_33) g_k

    //          = (Sf_1 C_1311 + G_111 sigma_13) g_1
    //            + (Sf_1 C_1321 + G_121 sigma_13) g_2
    //            + (Sf_1 C_1331 + G_131 sigma_13) g_3

    //            + (Sf_2 C_2311 + G_211 sigma_23) g_1
    //            + (Sf_2 C_2321 + G_221 sigma_23) g_2
    //            + (Sf_2 C_2331 + G_231 sigma_23) g_3

    //            + (Sf_3 C_3311 + G_311 sigma_33) g_1
    //            + (Sf_3 C_3321 + G_321 sigma_33) g_2
    //            + (Sf_3 C_3331 + G_331 sigma_33) g_3

    coeff[tensor::ZX] =
        (Sf[vector::X]*C(XZ,tensor::XX) + G(vector::X,tensor::XX)*sigma[symmTensor::XZ])*g[vector::X]
      + (Sf[vector::X]*C(XZ,tensor::YX) + G(vector::X,tensor::YX)*sigma[symmTensor::XZ])*g[vector::Y]
      + (Sf[vector::X]*C(XZ,tensor::ZX) + G(vector::X,tensor::ZX)*sigma[symmTensor::XZ])*g[vector::Z]

      + (Sf[vector::Y]*C(YZ,tensor::XX) + G(vector::Y,tensor::XX)*sigma[symmTensor::YZ])*g[vector::X]
      + (Sf[vector::Y]*C(YZ,tensor::YX) + G(vector::Y,tensor::YX)*sigma[symmTensor::YZ])*g[vector::Y]
      + (Sf[vector::Y]*C(YZ,tensor::ZX) + G(vector::Y,tensor::ZX)*sigma[symmTensor::YZ])*g[vector::Z]

      + (Sf[vector::Z]*C(ZZ,tensor::XX) + G(vector::Z,tensor::XX)*sigma[symmTensor::ZZ])*g[vector::X]
      + (Sf[vector::Z]*C(ZZ,tensor::YX) + G(vector::Z,tensor::YX)*sigma[symmTensor::ZZ])*g[vector::Y]
      + (Sf[vector::Z]*C(ZZ,tensor::ZX) + G(vector::Z,tensor::ZX)*sigma[symmTensor::ZZ])*g[vector::Z];

    // coeff_i2 = (Sf_m C_mik2 + G_mk2 sigma_mi) g_k delta_22
    //
    //          = (Sf_1 C_1ik2 + G_1k2 sigma_1i) g_k
    //            + (Sf_2 C_2ik2 + G_2k2 sigma_2i)  g_k
    //            + (Sf_3 C_3ik2 + G_3k2 sigma_3i) g_k
    //
    //          = (Sf_1 C_1i12 + G_112 sigma_1i) g_1
    //            + (Sf_1 C_1i22 + G_122 sigma_1i) g_2
    //            + (Sf_1 C_1i32 + G_132 sigma_1i) g_3
    //
    //            + (Sf_2 C_2i12 + G_212 sigma_2i) g_1
    //            + (Sf_2 C_2i22 + G_222 sigma_2i) g_2
    //            + (Sf_2 C_2i32 + G_232 sigma_2i) g_3
    //
    //            + (Sf_3 C_3i12 + G_312 sigma_3i) g_1
    //            + (Sf_3 C_3i22 + G_322 sigma_3i) g_2
    //            + (Sf_3 C_3i32 + G_332 sigma_3i) g_3

    // coeff_12 = (Sf_m C_m1k2 + G_mk2 sigma_m1) g_k
    //
    //          = (Sf_1 C_11k2 + G_1k2 sigma_11) g_k
    //            + (Sf_2 C_21k2 + G_2k2 sigma_21)  g_k
    //            + (Sf_3 C_31k2 + G_3k2 sigma_31) g_k
    //
    //          = (Sf_1 C_1112 + G_112 sigma_11) g_1
    //            + (Sf_1 C_1122 + G_122 sigma_11) g_2
    //            + (Sf_1 C_1132 + G_132 sigma_11) g_3
    //
    //            + (Sf_2 C_2112 + G_212 sigma_21) g_1
    //            + (Sf_2 C_2122 + G_222 sigma_21) g_2
    //            + (Sf_2 C_2132 + G_232 sigma_21) g_3
    //
    //            + (Sf_3 C_3112 + G_312 sigma_31) g_1
    //            + (Sf_3 C_3122 + G_322 sigma_31) g_2
    //            + (Sf_3 C_3132 + G_332 sigma_31) g_3

    coeff[tensor::XY] =
        (Sf[vector::X]*C(XX,tensor::XY) + G(vector::X,tensor::XY)*sigma[symmTensor::XX])*g[vector::X]
      + (Sf[vector::X]*C(XX,tensor::YY) + G(vector::X,tensor::YY)*sigma[symmTensor::XX])*g[vector::Y]
      + (Sf[vector::X]*C(XX,tensor::ZY) + G(vector::X,tensor::ZY)*sigma[symmTensor::XX])*g[vector::Z]

      + (Sf[vector::Y]*C(XY,tensor::XY) + G(vector::Y,tensor::XY)*sigma[symmTensor::XY])*g[vector::X]
      + (Sf[vector::Y]*C(XY,tensor::YY) + G(vector::Y,tensor::YY)*sigma[symmTensor::XY])*g[vector::Y]
      + (Sf[vector::Y]*C(XY,tensor::ZY) + G(vector::Y,tensor::ZY)*sigma[symmTensor::XY])*g[vector::Z]

      + (Sf[vector::Z]*C(XZ,tensor::XY) + G(vector::Z,tensor::XY)*sigma[symmTensor::XZ])*g[vector::X]
      + (Sf[vector::Z]*C(XZ,tensor::YY) + G(vector::Z,tensor::YY)*sigma[symmTensor::XZ])*g[vector::Y]
      + (Sf[vector::Z]*C(XZ,tensor::ZY) + G(vector::Z,tensor::ZY)*sigma[symmTensor::XZ])*g[vector::Z];

    // coeff_22 = (Sf_m C_m2k2 + G_mk2 sigma_m2) g_k
    //
    //          = (Sf_1 C_12k2 + G_1k2 sigma_12) g_k
    //            + (Sf_2 C_22k2 + G_2k2 sigma_22)  g_k
    //            + (Sf_3 C_32k2 + G_3k2 sigma_32) g_k
    //
    //          = (Sf_1 C_1212 + G_112 sigma_12) g_1
    //            + (Sf_1 C_1222 + G_122 sigma_12) g_2
    //            + (Sf_1 C_1232 + G_132 sigma_12) g_3
    //
    //            + (Sf_2 C_2212 + G_212 sigma_22) g_1
    //            + (Sf_2 C_2222 + G_222 sigma_22) g_2
    //            + (Sf_2 C_2232 + G_232 sigma_22) g_3
    //
    //            + (Sf_3 C_3212 + G_312 sigma_32) g_1
    //            + (Sf_3 C_3222 + G_322 sigma_32) g_2
    //            + (Sf_3 C_3232 + G_332 sigma_32) g_3

    coeff[tensor::YY] =
        (Sf[vector::X]*C(XY,tensor::XY) + G(vector::X,tensor::XY)*sigma[symmTensor::XY])*g[vector::X]
      + (Sf[vector::X]*C(XY,tensor::YY) + G(vector::X,tensor::YY)*sigma[symmTensor::XY])*g[vector::Y]
      + (Sf[vector::X]*C(XY,tensor::ZY) + G(vector::X,tensor::ZY)*sigma[symmTensor::XY])*g[vector::Z]

      + (Sf[vector::Y]*C(YY,tensor::XY) + G(vector::Y,tensor::XY)*sigma[symmTensor::YY])*g[vector::X]
      + (Sf[vector::Y]*C(YY,tensor::YY) + G(vector::Y,tensor::YY)*sigma[symmTensor::YY])*g[vector::Y]
      + (Sf[vector::Y]*C(YY,tensor::ZY) + G(vector::Y,tensor::ZY)*sigma[symmTensor::YY])*g[vector::Z]

      + (Sf[vector::Z]*C(YZ,tensor::XY) + G(vector::Z,tensor::XY)*sigma[symmTensor::YZ])*g[vector::X]
      + (Sf[vector::Z]*C(YZ,tensor::YY) + G(vector::Z,tensor::YY)*sigma[symmTensor::YZ])*g[vector::Y]
      + (Sf[vector::Z]*C(YZ,tensor::ZY) + G(vector::Z,tensor::ZY)*sigma[symmTensor::YZ])*g[vector::Z];

    // coeff_32 = (Sf_m C_m3k2 + G_mk2 sigma_m3) g_k
    //
    //          = (Sf_1 C_13k2 + G_1k2 sigma_13) g_k
    //            + (Sf_2 C_23k2 + G_2k2 sigma_23)  g_k
    //            + (Sf_3 C_33k2 + G_3k2 sigma_33) g_k
    //
    //          = (Sf_1 C_1312 + G_112 sigma_13) g_1
    //            + (Sf_1 C_1322 + G_122 sigma_13) g_2
    //            + (Sf_1 C_1332 + G_132 sigma_13) g_3
    //
    //            + (Sf_2 C_2312 + G_212 sigma_23) g_1
    //            + (Sf_2 C_2322 + G_222 sigma_23) g_2
    //            + (Sf_2 C_2332 + G_232 sigma_23) g_3
    //
    //            + (Sf_3 C_3312 + G_312 sigma_33) g_1
    //            + (Sf_3 C_3322 + G_322 sigma_33) g_2
    //            + (Sf_3 C_3332 + G_332 sigma_33) g_3

    coeff[tensor::ZY] =
        (Sf[vector::X]*C(XZ,tensor::XY) + G(vector::X,tensor::XY)*sigma[symmTensor::XZ])*g[vector::X]
      + (Sf[vector::X]*C(XZ,tensor::YY) + G(vector::X,tensor::YY)*sigma[symmTensor::XZ])*g[vector::Y]
      + (Sf[vector::X]*C(XZ,tensor::ZY) + G(vector::X,tensor::ZY)*sigma[symmTensor::XZ])*g[vector::Z]

      + (Sf[vector::Y]*C(YZ,tensor::XY) + G(vector::Y,tensor::XY)*sigma[symmTensor::YZ])*g[vector::X]
      + (Sf[vector::Y]*C(YZ,tensor::YY) + G(vector::Y,tensor::YY)*sigma[symmTensor::YZ])*g[vector::Y]
      + (Sf[vector::Y]*C(YZ,tensor::ZY) + G(vector::Y,tensor::ZY)*sigma[symmTensor::YZ])*g[vector::Z]

      + (Sf[vector::Z]*C(ZZ,tensor::XY) + G(vector::Z,tensor::XY)*sigma[symmTensor::ZZ])*g[vector::X]
      + (Sf[vector::Z]*C(ZZ,tensor::YY) + G(vector::Z,tensor::YY)*sigma[symmTensor::ZZ])*g[vector::Y]
      + (Sf[vector::Z]*C(ZZ,tensor::ZY) + G(vector::Z,tensor::ZY)*sigma[symmTensor::ZZ])*g[vector::Z];

    // coeff_i3 = (Sf_m C_mik3 + G_mk3 sigma_mi) g_k delta_33
    //
    //          = (Sf_1 C_1ik3 + G_1k3 sigma_1i) g_k
    //            + (Sf_2 C_2ik3 + G_2k3 sigma_2i) g_k
    //            + (Sf_3 C_3ik3 + G_3k3 sigma_3i) g_k
	//
	//			= (Sf_1 C_1i13 + G_113 sigma_1i) g_1
	//			  + (Sf_1 C_1i23 + G_123 sigma_1i) g_2
	//			  + (Sf_1 C_1i33 + G_133 sigma_1i) g_3
	//
	//			  + (Sf_2 C_2i13 + G_213 sigma_2i) g_1
	//			  + (Sf_2 C_2i23 + G_223 sigma_2i) g_2
	//			  + (Sf_2 C_2i33 + G_233 sigma_2i) g_3
	//
	//			  + (Sf_3 C_3i13 + G_313 sigma_3i) g_1
	//			  + (Sf_3 C_3i23 + G_323 sigma_3i) g_2
	//			  + (Sf_3 C_3i33 + G_333 sigma_3i) g_3

    // coeff_13 = (Sf_m C_m1k3 + G_mk3 sigma_m1) g_k

    //          = (Sf_1 C_11k3 + G_1k3 sigma_11) g_k
    //            + (Sf_2 C_21k3 + G_2k3 sigma_21) g_k
    //            + (Sf_3 C_31k3 + G_3k3 sigma_31) g_k
	//
	//			= (Sf_1 C_1113 + G_113 sigma_11) g_1
	//			  + (Sf_1 C_1123 + G_123 sigma_11) g_2
	//			  + (Sf_1 C_1133 + G_133 sigma_11) g_3
	//
	//			  + (Sf_2 C_2113 + G_213 sigma_21) g_1
	//			  + (Sf_2 C_2123 + G_223 sigma_21) g_2
	//			  + (Sf_2 C_2133 + G_233 sigma_21) g_3
	//
	//			  + (Sf_3 C_3113 + G_313 sigma_31) g_1
	//			  + (Sf_3 C_3123 + G_323 sigma_31) g_2
	//			  + (Sf_3 C_3133 + G_333 sigma_31) g_3

    coeff[tensor::XZ] =
        (Sf[vector::X]*C(XX,tensor::XZ) + G(vector::X,tensor::XZ)*sigma[symmTensor::XX])*g[vector::X]
      + (Sf[vector::X]*C(XX,tensor::YZ) + G(vector::X,tensor::YZ)*sigma[symmTensor::XX])*g[vector::Y]
      + (Sf[vector::X]*C(XX,tensor::ZZ) + G(vector::X,tensor::ZZ)*sigma[symmTensor::XX])*g[vector::Z]

      + (Sf[vector::Y]*C(XY,tensor::XZ) + G(vector::Y,tensor::XZ)*sigma[symmTensor::XY])*g[vector::X]
      + (Sf[vector::Y]*C(XY,tensor::YZ) + G(vector::Y,tensor::YZ)*sigma[symmTensor::XY])*g[vector::Y]
      + (Sf[vector::Y]*C(XY,tensor::ZZ) + G(vector::Y,tensor::ZZ)*sigma[symmTensor::XY])*g[vector::Z]

      + (Sf[vector::Z]*C(XZ,tensor::XZ) + G(vector::Z,tensor::XZ)*sigma[symmTensor::XZ])*g[vector::X]
      + (Sf[vector::Z]*C(XZ,tensor::YZ) + G(vector::Z,tensor::YZ)*sigma[symmTensor::XZ])*g[vector::Y]
      + (Sf[vector::Z]*C(XZ,tensor::ZZ) + G(vector::Z,tensor::ZZ)*sigma[symmTensor::XZ])*g[vector::Z];

    // coeff_23 = (Sf_m C_m2k3 + G_mk3 sigma_m2) g_k
    //
    //          = (Sf_1 C_12k3 + G_1k3 sigma_12) g_k
    //            + (Sf_2 C_22k3 + G_2k3 sigma_22) g_k
    //            + (Sf_3 C_32k3 + G_3k3 sigma_32) g_k
	//
	//			= (Sf_1 C_1213 + G_113 sigma_12) g_1
	//			  + (Sf_1 C_1223 + G_123 sigma_12) g_2
	//			  + (Sf_1 C_1233 + G_133 sigma_12) g_3
	//
	//			  + (Sf_2 C_2213 + G_213 sigma_22) g_1
	//			  + (Sf_2 C_2223 + G_223 sigma_22) g_2
	//			  + (Sf_2 C_2233 + G_233 sigma_22) g_3
	//
	//			  + (Sf_3 C_3213 + G_313 sigma_32) g_1
	//			  + (Sf_3 C_3223 + G_323 sigma_32) g_2
	//			  + (Sf_3 C_3233 + G_333 sigma_32) g_3

    coeff[tensor::YZ] =
        (Sf[vector::X]*C(XY,tensor::XZ) + G(vector::X,tensor::XZ)*sigma[symmTensor::XY])*g[vector::X]
      + (Sf[vector::X]*C(XY,tensor::YZ) + G(vector::X,tensor::YZ)*sigma[symmTensor::XY])*g[vector::Y]
      + (Sf[vector::X]*C(XY,tensor::ZZ) + G(vector::X,tensor::ZZ)*sigma[symmTensor::XY])*g[vector::Z]

      + (Sf[vector::Y]*C(YY,tensor::XZ) + G(vector::Y,tensor::XZ)*sigma[symmTensor::YY])*g[vector::X]
      + (Sf[vector::Y]*C(YY,tensor::YZ) + G(vector::Y,tensor::YZ)*sigma[symmTensor::YY])*g[vector::Y]
      + (Sf[vector::Y]*C(YY,tensor::ZZ) + G(vector::Y,tensor::ZZ)*sigma[symmTensor::YY])*g[vector::Z]

      + (Sf[vector::Z]*C(YZ,tensor::XZ) + G(vector::Z,tensor::XZ)*sigma[symmTensor::YZ])*g[vector::X]
      + (Sf[vector::Z]*C(YZ,tensor::YZ) + G(vector::Z,tensor::YZ)*sigma[symmTensor::YZ])*g[vector::Y]
      + (Sf[vector::Z]*C(YZ,tensor::ZZ) + G(vector::Z,tensor::ZZ)*sigma[symmTensor::YZ])*g[vector::Z];

    // coeff_33 = (Sf_m C_m3k3 + G_mk3 sigma_m3) g_k
    //
    //          = (Sf_1 C_13k3 + G_1k3 sigma_13) g_k
    //            + (Sf_2 C_23k3 + G_2k3 sigma_23) g_k
    //            + (Sf_3 C_33k3 + G_3k3 sigma_33) g_k
	//
	//			= (Sf_1 C_1313 + G_113 sigma_13) g_1
	//			  + (Sf_1 C_1323 + G_123 sigma_13) g_2
	//			  + (Sf_1 C_1333 + G_133 sigma_13) g_3
	//
	//			  + (Sf_2 C_2313 + G_213 sigma_23) g_1
	//			  + (Sf_2 C_2323 + G_223 sigma_23) g_2
	//			  + (Sf_2 C_2333 + G_233 sigma_23) g_3
	//
	//			  + (Sf_3 C_3313 + G_313 sigma_33) g_1
	//			  + (Sf_3 C_3323 + G_323 sigma_33) g_2
	//			  + (Sf_3 C_3333 + G_333 sigma_33) g_3

    coeff[tensor::ZZ] =
        (Sf[vector::X]*C(XZ,tensor::XZ) + G(vector::X,tensor::XZ)*sigma[symmTensor::XZ])*g[vector::X]
      + (Sf[vector::X]*C(XZ,tensor::YZ) + G(vector::X,tensor::YZ)*sigma[symmTensor::XZ])*g[vector::Y]
      + (Sf[vector::X]*C(XZ,tensor::ZZ) + G(vector::X,tensor::ZZ)*sigma[symmTensor::XZ])*g[vector::Z]

      + (Sf[vector::Y]*C(YZ,tensor::XZ) + G(vector::Y,tensor::XZ)*sigma[symmTensor::YZ])*g[vector::X]
      + (Sf[vector::Y]*C(YZ,tensor::YZ) + G(vector::Y,tensor::YZ)*sigma[symmTensor::YZ])*g[vector::Y]
      + (Sf[vector::Y]*C(YZ,tensor::ZZ) + G(vector::Y,tensor::ZZ)*sigma[symmTensor::YZ])*g[vector::Z]

      + (Sf[vector::Z]*C(ZZ,tensor::XZ) + G(vector::Z,tensor::XZ)*sigma[symmTensor::ZZ])*g[vector::X]
      + (Sf[vector::Z]*C(ZZ,tensor::YZ) + G(vector::Z,tensor::YZ)*sigma[symmTensor::ZZ])*g[vector::Y]
      + (Sf[vector::Z]*C(ZZ,tensor::ZZ) + G(vector::Z,tensor::ZZ)*sigma[symmTensor::ZZ])*g[vector::Z];

}


// ************************************************************************* //
