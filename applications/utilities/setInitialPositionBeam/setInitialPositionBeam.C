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

Application
    setInitialPositionBeam

Description
    This utility sets the initial configuration of the beams:

    There are two configurations of beam:
    1. Reference Config:
    - Beam mesh when created always lies such that the length direction of beam is aligned to x-axis.
    - The tangents at the face centres of the beam cells are set to (1 0 0) direction.

    2. Initial Configuration of beam:
    - To start with any other orientation of the beam, one needs to set the position vector of the beam curve.
    - The tangents at the face centres of the beam need to be set from (1 0 0) to the cross-section orientation
        of the beam spatial curve.

    * valid arguments are:
    * -cellZone : this command grabs the specififed subset of mesh cells
    * -translate vector
        Translate the points by the given vector before rotations
    * -rotateAngle (vector angle)
        Rotate angle degrees about vector axis.
    * -rotate : -rotate (vector vector)
        Rotate the points from the first vector to the second


\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "twoDPointCorrector.H"
#include "boundBox.H"
#include "pointFields.H"

//   newly added
#include "argList.H"
#include "ReadFields.H"
#include "IStringStream.H"
#include "transform.H"
#include "axisAngleRotation.H"


using namespace Foam;
using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "addRegionOption.H"
    argList::validOptions.insert("cellZone", "name");
    argList::validOptions.insert("translate", "vector");
    argList::validOptions.insert("rotateAngle", "(vector angleInDegree)");
    argList::validOptions.insert("rotate", "(vector vector)");

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    if (args.options().empty())
    {
        FatalErrorIn(args.executable())
            << "No options supplied, please use "
            << "-cellZone and -translate options "
			<< "to set the initial position of the beam\n"
			<< "-rotateAngle and rotate option is optional to use"
            << exit(FatalError);
    }

    surfaceVectorField  refTangent
    (
        IOobject
        (
            "refTangent",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedVector("x", dimless, vector(1, 0, 0))
    );

    // Beam mean line displacement
    surfaceVectorField refWf
    (
        IOobject
        (
            "refWf",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedVector("zero", dimLength, vector::zero)
    );

    // Beam mean line displacement
    volVectorField refW
    (
        IOobject
        (
            "refW",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedVector("zero", dimLength, vector::zero)
    );

    // Beam cross-section rotation
    surfaceTensorField refLambdaf
    (
        IOobject
        (
            "refLambdaf",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
         ),
        mesh,
        dimensionedTensor("I", dimless, tensor::I)
    );

    volTensorField refLambda
    (
        IOobject
        (
            "refLambda",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
         ),
        mesh,
        dimensionedTensor("I", dimless, tensor::I)
    );


    // Read cell zone
    word cellZoneName;

    if (args.found("cellZone"))
    {
        cellZoneName = args.get<word>("cellZone");

        Info << "Beam cell zone: " << cellZoneName << endl;
    }
    else
    {
        FatalErrorIn(args.executable())
            << "Option cellZone is not specified"
            << exit(FatalError);
    }

    // Translation and rotation options
    if (args.found("translate"))
    {
        vector transVector(args.lookup("translate")());

	    tensor T = tensor::I;

        if (args.found("rotate"))
        {
            Pair<vector> n1n2(args.lookup("rotate")());
            n1n2[0] /= mag(n1n2[0]);
            n1n2[1] /= mag(n1n2[1]);
            T = rotationTensor(n1n2[0], n1n2[1]);

            Info<< "Rotating points by " << T << endl;
        }
        if (args.found("rotateAngle"))
        {
            const Tuple2<vector, scalar> rotAxisAngle
            (
                args.lookup("rotateAngle")()
            );

            const vector& axis = rotAxisAngle.first();
            const scalar angle = rotAxisAngle.second();

            Info<< "Rotating points " << nl
                << "    about " << axis << nl
                << "    angle " << angle << nl;

            T = Foam::coordinateRotations::axisAngle::rotation(axis, angle, true);
            Info << "T=" << T << endl;
        }

        const label zoneID = mesh.cellZones().findZoneID(cellZoneName);

	    Info << "zoneID" << zoneID << endl;

	    if (zoneID == -1)
        {
            FatalErrorIn
            (
                "setInitialPositionBeam application utility"
            )
            << "Problem in beam cellZone"
            << "\nzoneID of beam: " << cellZoneName << " is " << zoneID
            << "\nProvide the beam name without the hyphen in front "
            << "e.g. beam_0"
            << abort(FatalError);
        }

        vectorField& refWfI = refWf.primitiveFieldRef();
        vectorField& refWI = refW.primitiveFieldRef();
        vectorField& refTangentI = refTangent.primitiveFieldRef();
        tensorField& refLambdafI = refLambdaf.primitiveFieldRef();
        tensorField& refLambdaI = refLambda.primitiveFieldRef();
        const vectorField& CfI = mesh.Cf().primitiveField();
        const vectorField& CI = mesh.C().primitiveField();

        vectorField newCfI ((T & CfI));
        const labelList& nei = mesh.neighbour();

        //- Set all the internal face fields
        forAll(refWfI, faceI)
        {
            label I = mesh.cellZones().whichZone(nei[faceI]);

            if (I == zoneID)
            {
                refWfI[faceI] =  newCfI[faceI] + transVector - CfI[faceI];
                refLambdafI[faceI] = T;
                // PC: refTang shoudl be wrt to the reference frame, not
                // incremental; otherwise, refTang would be inconsistent with
                // refWf and refW
                //refTangentI[faceI] = (refLambdafI[faceI] & refTangentI[faceI]);
                refTangentI[faceI] = (refLambdafI[faceI] & vector(1, 0, 0));
            }
        }

        vectorField newCI ((T & CI));

        //- Set all the cell-centre reference fields
        forAll(refWI, cellI)
        {
            label I  = mesh.cellZones().whichZone(cellI);

            if (I == zoneID)
            {
                refWI[cellI]  = newCI[cellI] + transVector - CI[cellI];

                refLambdaI[cellI] = T;
            }
        }

        //- Set the reference fields at the boundaries
        forAll(refWf.boundaryField(), patchI)
        {
            vectorField& pRefWf = refWf.boundaryFieldRef()[patchI];
            vectorField& pRefW = refW.boundaryFieldRef()[patchI];
            vectorField& pRefTangent = refTangent.boundaryFieldRef()[patchI];
            tensorField& pRefLambdaf = refLambdaf.boundaryFieldRef()[patchI];
            tensorField& pRefLambda = refLambda.boundaryFieldRef()[patchI];

            const vectorField& pCf = mesh.Cf().boundaryField()[patchI];

            const vectorField newpCf ((T & pCf));

            const labelList faceCells =
            mesh.boundary()[patchI].faceCells();

            forAll(pRefWf, faceI)
            {
                label I = mesh.cellZones().whichZone(faceCells[faceI]);

                if (I == zoneID)
                {
                    pRefWf[faceI] =  newpCf[faceI] + transVector - pCf[faceI];
                    pRefW[faceI] = pRefWf[faceI];

                    pRefLambdaf[faceI] = T;
                    pRefLambda[faceI] = pRefLambdaf[faceI];

                            // PC: see comment above
                    // pRefTangent[faceI] = (pRefLambdaf[faceI] & pRefTangent[faceI]);
                    pRefTangent[faceI] = (pRefLambdaf[faceI] & vector(1, 0, 0));
                }
            }
        }
    }
    else
    {
        FatalErrorIn(args.executable())
            << "Option translate is not specified"
            << exit(FatalError);
    }

    refWf.write();
    refW.write();
    refLambdaf.write();
    refTangent.write();
    refLambda.write();

    #include "updateMeshPoints.H"

    Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
        << "  ClockTime = " << runTime.elapsedClockTime() << " s"
        << nl << endl;

    Info<< "End\n" << endl;

    return(0);
}


// ************************************************************************* //
