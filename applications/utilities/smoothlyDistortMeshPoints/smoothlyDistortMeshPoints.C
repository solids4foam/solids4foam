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
    perturbMeshPoints

Description
    Move mesh points by the smooth function:

        pointD(x,y) =
            (
                Ax*sin(Bx*pi*x)*sin(Cx*pi*y),
                Ay*sin(By*pi*x)*sin(Cy*pi*y),
                Az*sin(Bz*pi*x)*sin(Cz*pi*y)
            )

    where A[x-z], B[x-z] and C[x-z] are user-provided parameters.

    Patches which should not move can be defined via the fixedPatches entry.

    The inputs are defined in $FOAM_CASE/system/smoothlyDistortMeshPointsDict

    This utility is useful for creating distorted grids for testing
    discretisations.

Author
    Philip Cardiff, UCD.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "argList.H"
#include "twoDPointCorrector.H"
#include "unitConversion.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
#   include "setRootCase.H"
#   include "createTime.H"
#   include "createMesh.H"

    argList::noParallel();

    // Read dictionary
    Info<< "Reading smoothlyDistortMeshPointsDict dictionary" << nl << endl;
    IOdictionary perturbDict
    (
        IOobject
        (
            "smoothlyDistortMeshPointsDict",
            runTime.system(),
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    // Read inputs
    const vector A(perturbDict.lookup("A"));
    const vector B(perturbDict.lookup("B"));
    const vector C(perturbDict.lookup("C"));
    const wordList fixedPatchesList(perturbDict.lookup("fixedPatches"));

    // Convert fixedPatches list to a set
    HashSet<word> fixedPatches;
    forAll(fixedPatchesList, pI)
    {
        fixedPatches.insert(fixedPatchesList[pI]);
    }

    // Store original points
    const pointField oldPoints = mesh.points();

    // Calculate a mask to identify fixed points
    boolList fixedPoint(oldPoints.size(), false);
    forAll(mesh.boundary(), patchI)
    {
        const word& patchName = mesh.boundary()[patchI].name();
        if (fixedPatches.found(patchName))
        {
            const labelList& meshPoints =
                mesh.boundaryMesh()[patchI].meshPoints();

            forAll(meshPoints, mpI)
            {
                const label pointID = meshPoints[mpI];
                fixedPoint[pointID] = true;
            }
        }
    }

    // Calculate new points
    pointField newPoints(oldPoints);

    const scalar Ax = A.x();
    const scalar Ay = A.y();
    const scalar Az = A.z();
    const scalar Bx = B.x();
    const scalar By = B.y();
    const scalar Bz = B.z();
    const scalar Cx = C.x();
    const scalar Cy = C.y();
    const scalar Cz = C.z();
    const scalar pi = constant::mathematical::pi;

    forAll(newPoints, pointI)
    {
        if (!fixedPoint[pointI])
        {
            const scalar x = oldPoints[pointI][vector::X];
            const scalar y = oldPoints[pointI][vector::Y];

            newPoints[pointI] +=
                vector
                (
                    Ax*Foam::sin(Bx*pi*x)*Foam::sin(Cx*pi*y),
                    Ay*Foam::sin(By*pi*x)*Foam::sin(Cy*pi*y),
                    Az*Foam::sin(Bz*pi*x)*Foam::sin(Cz*pi*y)
                );
        }
    }

    // Remove the normal component on boundary patches
    forAll(mesh.boundary(), patchI)
    {
        const pointField& pointNormals =
            mesh.boundaryMesh()[patchI].pointNormals();
        const labelList& meshPoints =
            mesh.boundaryMesh()[patchI].meshPoints();

        forAll(pointNormals, pI)
        {
            const vector& n = pointNormals[pI];
            const label pointID = meshPoints[pI];
            const vector disp = newPoints[pointID] - oldPoints[pointID];

            newPoints[pointID] = oldPoints[pointID] + ((I - sqr(n)) & disp);
        }
    }

    // Correct points for 2-D
    twoDPointCorrector twoD(mesh);
    twoD.correctPoints(newPoints);

    // Move the mesh
    Info<< "Applying the perturbation to the points" << endl;
    mesh.movePoints(newPoints);

    // Write the mesh
    Info<< "Writing the mesh" << endl;
    mesh.setInstance(mesh.polyMesh::instance());
    mesh.write();

    Info<< nl << "End" << nl << endl;

    return(0);
}


// ************************************************************************* //
