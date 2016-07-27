Info << " Reading solid mesh" << endl;

fvMesh solidMesh
(
    IOobject
    (
        "solid",
        runTime.timeName(),
        runTime,
        IOobject::MUST_READ
    )
);
