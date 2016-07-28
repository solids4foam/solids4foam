// Check if the dynamicMeshDict is found, if not, then we will add it setting
// the mesh to be a staticFvMesh
{
    IOobject dynamicMeshDictHeader
    (
        "dynamicMeshDict",
        runTime.constant(),
        runTime,
        IOobject::MUST_READ
    );

    if (!dynamicMeshDictHeader.headerOk())
    {
        WarningIn("checkDynamicMeshDict.H")
            << "dynamicMeshDict not found!" << nl
            << "    Adding a dynamicMeshDict and setting dynamicFvMesh to "
            << "staticFvMesh" << endl;

        // Create the dict
        IOdictionary dynamicMeshDict
        (
            IOobject
            (
                "dynamicMeshDict",
                runTime.constant(),
                runTime,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            )
        );

        dynamicMeshDict.add("dynamicFvMesh", "staticFvMesh");

        dynamicMeshDict.regIOobject::write();
    }
}
