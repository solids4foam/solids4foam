/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     |
    \\  /    A nd           | For copyright notice see file Copyright
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Application
    lduToCsrMatrix.H

Description
    Function to convert the LDU matix to a CSR (compressed row) matrix format.

    This code is adapted from SpeedIT plugin CSR_convert.H.

Author
    SpeedIT adapted by Philip Cardiff UCD/UT

\*---------------------------------------------------------------------------*/


#ifndef CSR_CONVERT_H
#define CSR_CONVERT_H

#include <vector>
#include "lduMatrix.H"


//
// Conversion from OpenFOAM lduMatrix to classic Compressed Sparse Row matrix.
//
template<class REAL_TYPE>
void ldu2csr
(
    //const Foam::fvScalarMatrix& matrix,
    const scalarField& diag,
    const scalarField& upper,
    const scalarField& lower,
    const unallocLabelList& upperAddr,
    const unallocLabelList& lowerAddr,
    std::vector< int >& c_idx,
    std::vector< int >& r_idx,
    std::vector< REAL_TYPE >& vals
)
{
    // Add boundary condition contributions to the diag
    //scalarField diag = matrix.diag();
    // THIS SHOULD BE DONE BEFORE CALLING THIS FUNCTION
    //matrix.addBoundaryDiag(diag, 0);


    int n_rows = diag.size();

    //
    // Calculate each row size. Sizes are shifted, because in the next part
    // array r_idx is modified twice.
    //
    for (int i=0; i < upperAddr.size(); i++)
    {
        int ri1 = lowerAddr[i];
        int ri2 = upperAddr[i];
        r_idx[ri1+2] ++;
        r_idx[ri2+2] ++;
    };

    //
    // Compute row offsets. Offsets are shifted by one positions,
    //  because they are used as START positions while filling values.
    //
    r_idx[0] = 0;
    r_idx[1] = 0;
    for(int i=1; i<n_rows; i++)
    {
        r_idx[i+1] += r_idx[i];
    };

    //
    // Fill in CSR matrix.
    // Order below is important to keep column indices sorted.
    //

    // lower triangle
    for (int i=0; i < lowerAddr.size(); i++)
    {
        int row    = upperAddr[i] +1;
        int column = lowerAddr[i];

        int idx = r_idx[row];
        vals[idx] = lower[i];
        c_idx[idx] = column;
        r_idx[row]++;
    };
    // diagonal
    for (int i=0; i<diag.size(); i++)
    {
        int idx = r_idx[i+1];
        vals[idx] = diag[i];
        c_idx[idx] = i; // i is row and column index
        r_idx[i+1]++;
    };
    // upper triangle
    for (int i=0; i < upperAddr.size(); i++)
    {
        int row    = lowerAddr[i] +1;
        int column = upperAddr[i];

        int idx = r_idx[row];
        vals[idx] = upper[i];
        c_idx[idx] = column;
        r_idx[row]++;
    };
};


#endif
