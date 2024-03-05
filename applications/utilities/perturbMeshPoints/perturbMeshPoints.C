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
    Add a random perturbation to each mesh points. Boundary points are not
    changed except for empty patches.

    This utility is useful for creating distorteg grids for testing
    discretisations.

    The inputs are defined in $FOAM_CASE/systemm/perturbMeshPointsDict, and
    consist of a seed (for the random number generator) and a scaling factor
    to scale the perturbations. The scaling factor is a vector to allow
    different scalings in different directions; for example, for 2-D, the Z
    component should be set to 0.0.

Author
    Philip Cardiff, UCD.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "argList.H"
#include "Random.H"
#include "twoDPointCorrector.H"

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
    Info<< "Reading perturbMeshPointsDict dictionary" << nl << endl;
    IOdictionary perturbDict
    (
        IOobject
        (
            "perturbMeshPointsDict",
            runTime.system(),
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    // Read inputs
    const scalar seed(readScalar(perturbDict.lookup("seed")));
    const vector scaleFactor(perturbDict.lookup("scaleFactor"));
#ifdef OPENFOAM_COM
    const Switch Gaussian(perturbDict.lookup("Gaussian"));
#endif

    // Create random number generator
    Random rnd(seed);

    // Calculate new points
    pointField newPoints(mesh.points());

    // Calculate a mask to identify boundary points, excluding points on empty
    // and wedge patches
    boolList boundaryPoint(newPoints.size(), false);
    forAll(mesh.boundary(), patchI)
    {
        if
        (
            mesh.boundary()[patchI].type() != "empty"
         && mesh.boundary()[patchI].type() != "wedge"
        )
        {
            const labelList& meshPoints =
                mesh.boundaryMesh()[patchI].meshPoints();

            forAll(meshPoints, mpI)
            {
                const label pointID = meshPoints[mpI];
                boundaryPoint[pointID] = true;
            }
        }
    }

    forAll(newPoints, pointI)
    {
        if (!boundaryPoint[pointI])
        {
#ifdef OPENFOAM_COM
            if (Gaussian)
            {
                // Gaussian distribution
                newPoints[pointI] +=
                    vector
                    (
                        scaleFactor.x()*rnd.GaussNormal<scalar>(),
                        scaleFactor.y()*rnd.GaussNormal<scalar>(),
                        scaleFactor.z()*rnd.GaussNormal<scalar>()
                    );
            }
            else
#endif
            {
                // Uniform distribution
                newPoints[pointI] +=
                    vector
                    (
#ifdef FOAMEXTEND
                        scaleFactor.x()*(2.0*rnd.scalar01() - 1.0),
                        scaleFactor.y()*(2.0*rnd.scalar01() - 1.0),
                        scaleFactor.z()*(2.0*rnd.scalar01() - 1.0)
#else
                        scaleFactor.x()*(2.0*rnd.sample01<scalar>() - 1.0),
                        scaleFactor.y()*(2.0*rnd.sample01<scalar>() - 1.0),
                        scaleFactor.z()*(2.0*rnd.sample01<scalar>() - 1.0)
#endif
                    );
            }
        }
    }

    // Correct points for 2-D
    twoDPointCorrector twoD(mesh);
    twoD.correctPoints(newPoints);

    // Move the mesh
    mesh.movePoints(newPoints);

    // Write the mesh
    Info<< "Writing the mesh" << endl;
    mesh.setInstance(mesh.polyMesh::instance());
    mesh.write();

    Info<< nl << "End" << nl << endl;

    return(0);
}


// ************************************************************************* //
