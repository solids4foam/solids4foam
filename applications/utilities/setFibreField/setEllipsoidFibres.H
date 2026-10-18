// The ellipsoid method of setFibreField: the fibre (f0), sheet (s0) and
// sheet-normal (n0) fields of a truncated ellipsoid, from the transmural
// coordinate t, after Arostica et al. (2025). See setFibreField.C.
//
// Created by Claude Opus 5.5, supervised by Chanikya Valeti, UCD.
{
    ellipsoidParams p;
    p.rsEndo = readScalar(dict.lookup("rShortEndo"));
    p.rlEndo = readScalar(dict.lookup("rLongEndo"));
    p.rsEpi = readScalar(dict.lookup("rShortEpi"));
    p.rlEpi = readScalar(dict.lookup("rLongEpi"));
    p.alphaEndo = degToRad(readScalar(dict.lookup("alphaEndo")));
    p.alphaEpi = degToRad(readScalar(dict.lookup("alphaEpi")));

    Info<< "Ellipsoid: endo (" << p.rsEndo << ", " << p.rlEndo << "), epi ("
        << p.rsEpi << ", " << p.rlEpi << "), alpha from "
        << radToDeg(p.alphaEndo) << " to " << radToDeg(p.alphaEpi)
        << " degrees" << endl;

    const scalarField& tI = Foam::primitiveField(t);

    Info<< "    t in [" << gMin(tI) << ", " << gMax(tI) << "]" << endl;

    volVectorField f0
    (
        IOobject
        (
            "f0",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedVector("zero", dimless, vector::zero),
        "calculated"
    );
    volVectorField s0("s0", f0);
    volVectorField n0("n0", f0);

    const volVectorField& C = mesh.C();

    forAll(C, cellI)
    {
        ellipsoidBasis(p, C[cellI], tI[cellI], f0[cellI], s0[cellI], n0[cellI]);
    }

    auto& f0B = boundaryFieldRef(f0);
    auto& s0B = boundaryFieldRef(s0);
    auto& n0B = boundaryFieldRef(n0);

    forAll(mesh.boundary(), patchI)
    {
        const vectorField& Cf = C.boundaryField()[patchI];
        const scalarField& tP = t.boundaryField()[patchI];

        forAll(Cf, faceI)
        {
            ellipsoidBasis
            (
                p,
                Cf[faceI],
                tP[faceI],
                f0B[patchI][faceI],
                s0B[patchI][faceI],
                n0B[patchI][faceI]
            );
        }
    }

    Info<< "Writing f0, s0 and n0" << endl;
    f0.write();
    s0.write();
    n0.write();
}
