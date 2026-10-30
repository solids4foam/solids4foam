{
    const vectorField& refWfI = refWf.internalField();
    const tensorField& refLambdafI = refLambdaf.internalField();
    const faceList& faces = mesh.faces();
    fileName meshDir = polyMesh::meshSubDir;

    pointIOField points
    (
        IOobject
        (
            "points",
            runTime.constant(),
            meshDir,
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    );

    pointMesh pMesh(mesh);

    pointVectorField pointW
    (
        IOobject
        (
            "pointW",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        pMesh,
        dimensionedVector("0", dimLength, vector::zero)
    );

    vectorField& pointWI = pointW.primitiveFieldRef();

    forAll(refWfI, faceI)
    {
        const face& curFace = faces[faceI];

        vector C0 = curFace.centre(points);

        forAll(curFace, pointI)
        {
            label curPoint = curFace[pointI];

            vector oldR = points[curPoint] - C0;

            vector newR = C0 + refWfI[faceI] + (refLambdafI[faceI] & oldR);

            pointWI[curPoint] = newR - points[curPoint];
        }
    }

    forAll(refWf.boundaryField(), patchI)
    {
        const vectorField& pRefWf =
            refWf.boundaryField()[patchI];

        const tensorField& pRefLambdaf =
            refLambdaf.boundaryField()[patchI];

        const label start =
            mesh.boundaryMesh()[patchI].start();

        forAll(pRefWf, faceI)
        {
            const face& curFace = faces[start + faceI];
            vector C0 = curFace.centre(points);

            forAll(curFace, pointI)
            {
                label curPoint = curFace[pointI];

                vector oldR = points[curPoint] - C0;

                vector newR = C0 + pRefWf[faceI] + (pRefLambdaf[faceI] & oldR);

                pointWI[curPoint] = newR - points[curPoint];
            }
        }
    }

    pointW.write();
    Info << "Done setting Initial Position for Beam" << endl;
}
