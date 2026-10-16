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

Application
    Test-cellMoments

Description
    Prints the volume-normalised central first-, second- and third-order
    moments for every cell of an arbitrary two- or three-dimensional case
    mesh.

    The component ordering is

        first order : (x y z)
        second order: (xx xy xz yy yz zz)
        third order : (xxx xxy xxz xyy xyz xzz yyy yyz yzz zzz)

    For an axis-aligned hexahedron that is symmetric about its cell centre,
    the expected first-order moments are zero. With side lengths h_x, h_y and
    h_z, the expected second-order moments are

        (h_x^2/12 0 0 h_y^2/12 0 h_z^2/12)

    and all third-order moments are zero. For a two-dimensional mesh with z as
    the empty direction, all components containing z are zero. A rotated
    hexahedron can have non-zero second-order cross moments when expressed in
    the global coordinate system.

Author
    Ivan Batistic, UCD.
    Philip Cardiff, UCD.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "compatibilityFunctions.H"
#include "fvMeshQuadrature.H"
#include "globalIndex.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    if (mesh.nGeometricD() != 2 && mesh.nGeometricD() != 3)
    {
        FatalErrorInFunction
            << "Only two- and three-dimensional meshes are supported"
            << abort(FatalError);
    }

    if (mesh.nGeometricD() == 2 && mesh.solutionD()[vector::Z] != -1)
    {
        FatalErrorInFunction
            << "For a two-dimensional mesh, the empty direction must be "
            << "vector::Z" << abort(FatalError);
    }

    // Third-order quadrature provides all requested moment fields
    fvMeshQuadrature quadrature(mesh, 3, 2, true);

    const List<vector>& firstMoments = quadrature.firstOrderCellMoments();
    const List<symmTensor>& secondMoments =
        quadrature.secondOrderCellMoments();
    const List<symmTensor3rdOrder>& thirdMoments =
        quadrature.thirdOrderCellMoments();
    auto& cellQuadWeights =
        compactListListCRef(quadrature.cellQuadWeights());

    const vectorField& cellCentres = mesh.C();
    const scalarField& cellVolumes = mesh.V();

    scalarField cellVolumeRelativeErrors(mesh.nCells(), 0.0);
    scalarField cellQuadWeightSums(mesh.nCells(), 0.0);

    scalar sumFirstOrderError = 0.0;
    scalar maxFirstOrderError = 0.0;
    scalar sumCellVolumeRelativeError = 0.0;
    scalar maxCellVolumeRelativeError = 0.0;
    label nCells = firstMoments.size();

    forAll(firstMoments, cellI)
    {
        const scalar firstOrderError = mag(firstMoments[cellI]);

        sumFirstOrderError += firstOrderError;
        maxFirstOrderError = max(maxFirstOrderError, firstOrderError);

        forAll(cellQuadWeights[cellI], qI)
        {
            cellQuadWeightSums[cellI] += cellQuadWeights[cellI][qI];
        }

        cellVolumeRelativeErrors[cellI] =
            mag(cellQuadWeightSums[cellI] - cellVolumes[cellI])
           /cellVolumes[cellI];

        sumCellVolumeRelativeError += cellVolumeRelativeErrors[cellI];
        maxCellVolumeRelativeError =
            max
            (
                maxCellVolumeRelativeError,
                cellVolumeRelativeErrors[cellI]
            );
    }

    reduce(sumFirstOrderError, sumOp<scalar>());
    reduce(maxFirstOrderError, maxOp<scalar>());
    reduce(sumCellVolumeRelativeError, sumOp<scalar>());
    reduce(maxCellVolumeRelativeError, maxOp<scalar>());
    reduce(nCells, sumOp<label>());

    const scalar averageFirstOrderError =
        nCells > 0 ? sumFirstOrderError/nCells : 0.0;
    const scalar averageCellVolumeRelativeError =
        nCells > 0 ? sumCellVolumeRelativeError/nCells : 0.0;

    Info<< nl
        << "Cell moment component ordering:" << nl
        << "    firstOrder  (x y z)" << nl
        << "    secondOrder (xx xy xz yy yz zz)" << nl
        << "    thirdOrder  (xxx xxy xxz xyy xyz xzz yyy yyz yzz zzz)"
        << nl << endl;

    forAll(cellCentres, cellI)
    {
        Info<< " cell: " << cellI << nl
            << " centre " << cellCentres[cellI] << nl
            << " cellVolume " << cellVolumes[cellI] << nl
            << " cellQuadWeightSum " << cellQuadWeightSums[cellI] << nl
            << " cellVolumeRelativeError "
            << cellVolumeRelativeErrors[cellI] * 100 << "%" << nl
            << " firstOrder " << firstMoments[cellI] << nl
            << " secondOrder " << secondMoments[cellI] << nl
            << " thirdOrder " << thirdMoments[cellI] << nl << endl;
    }

    Info<< "First-order central moment errors:" << nl
        << "    average magnitude: " << averageFirstOrderError << nl
        << "    maximum magnitude: " << maxFirstOrderError << nl << endl;

    Info<< "Cell quadrature weight volume errors:" << nl
        << "    average relative error: "
        << averageCellVolumeRelativeError * 100 << "%" << nl
        << "    maximum relative error: "
        << maxCellVolumeRelativeError * 100 << "%" <<nl << endl;

    return 0;
}


// ************************************************************************* //
