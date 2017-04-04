/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     3.2
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
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

\*---------------------------------------------------------------------------*/

#include "newLeastSquaresVolPointInterpolation.H"

namespace Foam
{

template <class Type>
void volPointInterpolateField
(
    const IOobject& header,
    const fvMesh& mesh,
    const pointMesh& pMesh,
    const newLeastSquaresVolPointInterpolation& volToPointInterp,
    bool& processed
)
{
    typedef GeometricField<Type, fvPatchField, volMesh> fieldType;
    typedef GeometricField<Type, pointPatchField, pointMesh> pFieldType;

    if (header.headerClassName() == fieldType::typeName)
    {
        Info<< "    Reading " << header.name() << endl;
        fieldType vf(header, mesh);

        // Create output point field
        pFieldType pf
        (
            IOobject
            (
                "volToPoint(" + header.name() + ')',
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ
            ),
            pMesh,
            dimensioned<Type>("zero", vf.dimensions(), pTraits<Type>::zero)
        );

        volToPointInterp.interpolate(vf, pf);

        Info<< "    Writing " << pf.name() << endl;
        pf.write();

        processed = true;
    }
}

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
