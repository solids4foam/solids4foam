/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     |
    \\  /    A nd           | For copyright notice see file Copyright
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

Description
    Transform vectorField into diagTensorField

Author
    Zeljko Tukovic, FSB Zagreb

\*---------------------------------------------------------------------------*/

#include "diagVectorField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

    tmp<diagTensorField> diag(const vectorField& vf)
{
    // Prepare the result field
    tmp<diagTensorField> tresult
    (
        new diagTensorField(vf.size(), diagTensor::zero)
    );
    diagTensorField& dtf = tresult();

    forAll(dtf, i)
    {
        dtf[i].xx() = vf[i].x();
        dtf[i].yy() = vf[i].y();
        dtf[i].zz() = vf[i].z();
    }

    return tresult;
}

} // End namespace Foam

// ************************************************************************* //
