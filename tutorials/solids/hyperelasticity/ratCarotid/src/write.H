    // if (runTime.outputTime())
    {
        // Cell fields

        volVectorField Er
        (
            IOobject
            (
                "Er",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            fvc::grad(psi0)
        );
        Er /= mag(Er);

        volVectorField Ea
        (
            IOobject
            (
                "Ea",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            fvc::grad(psi1)
        );
        Ea /= mag(Ea);

        volVectorField Ec
        (
            IOobject
            (
                "Ec",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            Ea^Er
        );
        Ec /= mag(Ec);

        Ea.write();
        Ec.write();
        Er.write();

        // Surface fields

        surfaceVectorField Erf
        (
            IOobject
            (
                "Erf",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            fvc::interpolate(Er)
        );
        Erf /= mag(Erf);

        surfaceVectorField Eaf
        (
            IOobject
            (
                "Eaf",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            fvc::interpolate(Ea)
        );
        Eaf /= mag(Eaf);

        surfaceVectorField Ecf
        (
            IOobject
            (
                "Ecf",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            fvc::interpolate(Ec)
        );
        Ecf /= mag(Ecf);

        Eaf.write();
        Ecf.write();
        Erf.write();

        psi0.write();
        psi1.write();

        // runTime.write();
    }
