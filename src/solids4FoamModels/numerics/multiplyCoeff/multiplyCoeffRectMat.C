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
        Sf[vector::X]*C(0,0)*g[vector::X]
      + Sf[vector::Y]*C(0,3)*g[vector::Y]
      + Sf[vector::Z]*C(0,5)*g[vector::Z]
      
      + Sf[vector::X]*C(3,0)*g[vector::X]
      + Sf[vector::Y]*C(3,3)*g[vector::Y]
      + Sf[vector::Z]*C(3,5)*g[vector::Z]
      
      + Sf[vector::X]*C(5,0)*g[vector::X]
      + Sf[vector::Y]*C(5,3)*g[vector::Y]
      + Sf[vector::Z]*C(5,5)*g[vector::Z];

    // Similarly the other components can be calculated as
    // coeff_21 = Sf_1 C_1211 g_1
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
        Sf[vector::X]*C(3,0)*g[vector::X]
      + Sf[vector::Y]*C(3,3)*g[vector::Y]
      + Sf[vector::Z]*C(3,5)*g[vector::Z]
      
      + Sf[vector::X]*C(1,0)*g[vector::X]
      + Sf[vector::Y]*C(1,3)*g[vector::Y]
      + Sf[vector::Z]*C(1,5)*g[vector::Z]
      
      + Sf[vector::X]*C(4,0)*g[vector::X]
      + Sf[vector::Y]*C(4,3)*g[vector::Y]
      + Sf[vector::Z]*C(4,5)*g[vector::Z];    

    // coeff_31 = Sf_1 C_1311 g_1
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
        Sf[vector::X]*C(5,0)*g[vector::X]
      + Sf[vector::Y]*C(5,3)*g[vector::Y]
      + Sf[vector::Z]*C(5,5)*g[vector::Z]
      
      + Sf[vector::X]*C(4,0)*g[vector::X]
      + Sf[vector::Y]*C(4,3)*g[vector::Y]
      + Sf[vector::Z]*C(4,5)*g[vector::Z]
      
      + Sf[vector::X]*C(2,0)*g[vector::X]
      + Sf[vector::Y]*C(2,3)*g[vector::Y]
      + Sf[vector::Z]*C(2,5)*g[vector::Z];

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

    // coeff_12 = Sf_1 C_1112 g_1
    //            + Sf_1 C_1122 g_2
    //            + Sf_1 C_1132 g_3
    //
    //            + Sf_2 C_2112 g_1
    //            + Sf_2 C_2122 g_2
    //            + Sf_2 C_2132 g_3
    //
    //            + Sf_3 C_3112 g_1
    //            + Sf_3 C_3122 g_2
    //            + Sf_3 C_3132 g_3

    // coeff_12 = Sf_2 C_2112 g_1 + Sf_1 C_1122 g_2
    
    coeff[tensor::XY] = 
        Sf[vector::X]*C(0,3)*g[vector::X]
      + Sf[vector::Y]*C(0,1)*g[vector::Y]
      + Sf[vector::Z]*C(0,4)*g[vector::Z]
      
      + Sf[vector::X]*C(3,3)*g[vector::X]
      + Sf[vector::Y]*C(3,1)*g[vector::Y]
      + Sf[vector::Z]*C(3,4)*g[vector::Z]
      
      + Sf[vector::X]*C(5,3)*g[vector::X]
      + Sf[vector::Y]*C(5,1)*g[vector::Y]
      + Sf[vector::Z]*C(5,4)*g[vector::Z];
      
    // coeff_22 = Sf_1 C_1212 g_1
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
        Sf[vector::X]*C(3,3)*g[vector::X]
      + Sf[vector::Y]*C(3,1)*g[vector::Y]
      + Sf[vector::Z]*C(3,4)*g[vector::Z]
      
      + Sf[vector::X]*C(1,3)*g[vector::X]
      + Sf[vector::Y]*C(1,1)*g[vector::Y]
      + Sf[vector::Z]*C(1,4)*g[vector::Z]
      
      + Sf[vector::X]*C(4,3)*g[vector::X]
      + Sf[vector::Y]*C(4,1)*g[vector::Y]
      + Sf[vector::Z]*C(4,4)*g[vector::Z];

    // coeff_32 = Sf_1 C_1312 g_1
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
    
    coeff[tensor::YZ] = 
        Sf[vector::X]*C(5,3)*g[vector::X]
      + Sf[vector::Y]*C(5,1)*g[vector::Y]
      + Sf[vector::Z]*C(5,4)*g[vector::Z]
      
      + Sf[vector::X]*C(3,3)*g[vector::X]
      + Sf[vector::Y]*C(3,1)*g[vector::Y]
      + Sf[vector::Z]*C(3,4)*g[vector::Z]
      
      + Sf[vector::X]*C(2,3)*g[vector::X]
      + Sf[vector::Y]*C(2,1)*g[vector::Y]
      + Sf[vector::Z]*C(2,4)*g[vector::Z];

    // coeff_i3 = Sf_m C_mik3 g_k delta_33
    //
    //          = Sf_1 C_1ik3 g_k
    //          + Sf_2 C_2ik3 g_k
    //          + Sf_3 C_3ik3 g_k
    //
    // coeff_13 = Sf_1 C_1113 g_1
    //          + Sf_1 C_1123 g_2
    //          + Sf_1 C_1133 g_3
    //
    //          + Sf_2 C_2113 g_1
    //          + Sf_2 C_2123 g_2
    //          + Sf_2 C_2133 g_3
    //
    //          + Sf_3 C_3113 g_1
    //          + Sf_3 C_3123 g_2
    //          + Sf_3 C_3133 g_3
       
    coeff[tensor::XZ] = 
        Sf[vector::X]*C(0,5)*g[vector::X]
      + Sf[vector::Y]*C(0,4)*g[vector::Y]
      + Sf[vector::Z]*C(0,2)*g[vector::Z]
      
      + Sf[vector::X]*C(3,5)*g[vector::X]
      + Sf[vector::Y]*C(3,4)*g[vector::Y]
      + Sf[vector::Z]*C(3,2)*g[vector::Z]
      
      + Sf[vector::X]*C(5,5)*g[vector::X]
      + Sf[vector::Y]*C(5,4)*g[vector::Y]
      + Sf[vector::Z]*C(5,2)*g[vector::Z];
      
    // coeff_23 = Sf_1 C_1213 g_1
    //          + Sf_1 C_1223 g_2
    //          + Sf_1 C_1233 g_3
    //
    //          + Sf_2 C_2213 g_1
    //          + Sf_2 C_2223 g_2
    //          + Sf_2 C_2233 g_3
    //
    //          + Sf_3 C_3213 g_1
    //          + Sf_3 C_3223 g_2
    //          + Sf_3 C_3233 g_3

    coeff[tensor::XZ] = 
        Sf[vector::X]*C(3,5)*g[vector::X]
      + Sf[vector::Y]*C(3,4)*g[vector::Y]
      + Sf[vector::Z]*C(3,2)*g[vector::Z]
      
      + Sf[vector::X]*C(1,5)*g[vector::X]
      + Sf[vector::Y]*C(1,4)*g[vector::Y]
      + Sf[vector::Z]*C(1,2)*g[vector::Z]
      
      + Sf[vector::X]*C(4,5)*g[vector::X]
      + Sf[vector::Y]*C(4,4)*g[vector::Y]
      + Sf[vector::Z]*C(4,2)*g[vector::Z];

    // coeff_33 = Sf_1 C_1313 g_1
    //          + Sf_1 C_1323 g_2
    //          + Sf_1 C_1333 g_3
    //
    //          + Sf_2 C_2313 g_1
    //          + Sf_2 C_2323 g_2
    //          + Sf_2 C_2333 g_3
    //
    //          + Sf_3 C_3313 g_1
    //          + Sf_3 C_3323 g_2
    //          + Sf_3 C_3333 g_3
    // Note: C_1313 == C_3131
    
    coeff[tensor::ZZ] = 
        Sf[vector::X]*C(5,5)*g[vector::X]
      + Sf[vector::Y]*C(5,4)*g[vector::Y]
      + Sf[vector::Z]*C(5,2)*g[vector::Z]
      
      + Sf[vector::X]*C(4,5)*g[vector::X]
      + Sf[vector::Y]*C(4,4)*g[vector::Y]
      + Sf[vector::Z]*C(4,2)*g[vector::Z]
      
      + Sf[vector::X]*C(2,5)*g[vector::X]
      + Sf[vector::Y]*C(2,4)*g[vector::Y]
      + Sf[vector::Z]*C(2,2)*g[vector::Z];
      
}

// ************************************************************************* //
