Info<< "Create solid dynamic mesh for time = "
    << runTime.timeName() << nl << endl;

autoPtr<dynamicFvMesh> solidMeshPtr
(
    dynamicFvMesh::New
    (
            IOobject
            (
                "solid",
                runTime.timeName(),
                runTime,
                IOobject::MUST_READ
            )
    )
);

dynamicFvMesh& solidMesh = solidMeshPtr();
