/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     |
    \\  /    A nd           | For copyright notice see file Copyright
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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

#include "pointGaussLsDivSigmaScheme.H"

#include "blockTangentialCoeffs.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{


// * * * * * * * * * * Private member functions * * * * * * * * * * * * * * * //


void pointGaussLsDivSigmaScheme::insertCoeffsTang
(
    const solidPolyMesh& solidMesh,
    const surfaceScalarField& muf,
    const surfaceScalarField& lambdaf,
    const GeometricField<vector, fvPatchField, volMesh>& vf,
    Field<vector>& blockB,
    BlockLduMatrix<vector>& blockM
)
{
    if (debug)
    {
        Info<< nl << type()
            << nl << "    start time: " << solidMesh().time().elapsedClockTime()
            << " s" << endl;
    }

    blockFvmInsertCoeffsTang
    (
        solidMesh,
        muf,
        lambdaf,
        vf,
        blockB,
        blockM,
        volToPointInterp(),
        3    // laplacian + transpose + trace
    );

    if (debug)
    {
        Info<< nl << type()
            << "    end time: " << solidMesh().time().elapsedClockTime() << " s"
            << nl << endl;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
