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

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Type Foam::quadrilateral::integral
(
    const Field<vector>& points,
    const Field<Type>& fld
) const
{
    Type integ = pTraits<Type>::zero;

    // Gauss integration for quad face with 4 integration points
    // See: https://link.springer.com/content/pdf
    // /bbm%3A978-3-540-32609-0%2F1.pdf
    // Page: 475

    // Gauss integration points
    scalarField gXi(4, 0);
    gXi[0] = -1.0/::sqrt(3);
    gXi[1] =  1.0/::sqrt(3);
    gXi[2] = -1.0/::sqrt(3);
    gXi[3] =  1.0/::sqrt(3);

    scalarField gEta(4, 0);
    gEta[0] = -1.0/::sqrt(3);
    gEta[1] = -1.0/::sqrt(3);
    gEta[2] =  1.0/::sqrt(3);
    gEta[3] =  1.0/::sqrt(3);

    Field<Type> gFld(4, pTraits<Type>::zero);
    forAll(gFld, gpI)
    {
        forAll(*this, pI)
        {
            gFld[gpI] +=
                fld[operator[](pI)]*Ni(pI, gXi[gpI], gEta[gpI]);
        }
    }

    Field<scalar> gMagJacob(4, 0);
    forAll(gMagJacob, gpI)
    {
        gMagJacob[gpI] =
            mag(detJacob(gXi[gpI], gEta[gpI], points));
    }

    forAll(gFld, gpI)
    {
        integ += gFld[gpI]*gMagJacob[gpI];
    }

    return integ;
}

template<class Type>
Type Foam::quadrilateral::integral
(
    const Field<vector>& points,
    const Field<Type>& fld0,
    const Field<scalar>& fld1
) const
{
    Type integ = pTraits<Type>::zero;

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

    // Info << gXi << endl;
    // Info << gEta << endl;

    Field<Type> gFld0(4, pTraits<Type>::zero);
    forAll(gFld0, gpI)
    {
        forAll(*this, pI)
        {
            gFld0[gpI] +=
                fld0[operator[](pI)]*Ni(pI, gXi[gpI], gEta[gpI]);
        }
    }

    Field<scalar> gFld1(4, 0);
    forAll(gFld1, gpI)
    {
        forAll(*this, pI)
        {
            gFld1[gpI] +=
                fld1[operator[](pI)]*Ni(pI, gXi[gpI], gEta[gpI]);
        }
    }

    Field<scalar> gMagJacob(4, 0);
    forAll(gMagJacob, gpI)
    {
        gMagJacob[gpI] =
            mag(detJacob(gXi[gpI], gEta[gpI], points));
    }

    forAll(gFld0, gpI)
    {
        integ += gFld0[gpI]*gFld1[gpI]*gMagJacob[gpI];

        // Info << gFld[gpI] << ", " << gMagJacob[gpI] << endl;
    }

    return integ;
}

// ************************************************************************* //
