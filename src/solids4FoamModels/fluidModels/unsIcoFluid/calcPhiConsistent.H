volVectorField HUcorr("HUcorr", HUEqn.H());

if (true)
{
    HUcorr.internalField() -=
        (AU.internalField() - AUcorr.internalField())
        *U().internalField();

    HUcorr.correctBoundaryConditions();
}

phi() = rAUf_*(fvc::interpolate(HUcorr, "interpolate(U)") & mesh.Sf());

// Add temporal term
word ddtScheme
(
    mesh.schemesDict().ddtScheme("ddt(" + U().name() +')')
);

if
(
    (
        ddtScheme
     == fv::EulerDdtScheme<vector>::typeName
    )
)
{
    if (mesh.moving())
    {
        volScalarField V0oV
        (
            IOobject
            (
                "V0oV",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimless,
            zeroGradientFvPatchScalarField::typeName
        );

        V0oV.internalField() = mesh.V0()/mesh.V();
        V0oV.correctBoundaryConditions();

        phi() += rAUf_*fvc::interpolate(V0oV)
                   *(Uf_.oldTime() & mesh.Sf())
                   /runTime().deltaT();
    }
    else
    {
        phi() += rAUf_*phi().oldTime()/runTime().deltaT();
    }
}
else if
(
    ddtScheme
 == fv::backwardDdtScheme<vector>::typeName
)
{
    if (mesh.moving())
    {
        volScalarField V0oV
        (
            IOobject
            (
                "V0oV",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimless,
            zeroGradientFvPatchScalarField::typeName
        );

        V0oV.internalField() = mesh.V0()/mesh.V();
        V0oV.correctBoundaryConditions();

        volScalarField V00oV
        (
            IOobject
            (
                "V00oV",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimless,
            zeroGradientFvPatchScalarField::typeName
        );

        V00oV.internalField() = mesh.V00()/mesh.V();
        V00oV.correctBoundaryConditions();

        if(runTime().timeIndex() == 1)
        {
            phi() += rAUf_*fvc::interpolate(V0oV)
               *(Uf_.oldTime() & mesh.Sf())
               /runTime().deltaT();
        }
        else
        {
            dimensionedScalar Cn("Cn", dimless/dimTime, 0.0);
            dimensionedScalar Co("Co", dimless/dimTime, 0.0);
            dimensionedScalar Coo("Coo", dimless/dimTime, 0.0);

            scalar deltaT = runTime().deltaT().value();
            scalar deltaT0 = runTime().deltaT0().value();

            Cn.value() = 1 + deltaT/(deltaT + deltaT0);
            Coo.value() = deltaT*deltaT/(deltaT0*(deltaT + deltaT0));
            Co.value() = Cn.value() + Coo.value();

            Cn.value() /= deltaT;
            Co.value() /= deltaT;
            Coo.value() /= deltaT;

            phi() +=
                rAUf_
               *(
                    Co*fvc::interpolate(V0oV)
                   *(Uf_.oldTime() & mesh.Sf())
                  - Coo*fvc::interpolate(V00oV)
                   *(Uf_.oldTime().oldTime() & mesh.Sf())
                );
        }
    }
    else
    {
        if(runTime().timeIndex() == 1)
        {
            phi() += rAUf_*phi().oldTime()/runTime().deltaT();
        }
        else
        {
            dimensionedScalar Cn("Cn", dimless/dimTime, 0.0);
            dimensionedScalar Co("Co", dimless/dimTime, 0.0);
            dimensionedScalar Coo("Coo", dimless/dimTime, 0.0);

            scalar deltaT = runTime().deltaT().value();
            scalar deltaT0 = runTime().deltaT0().value();

            Cn.value() = 1 + deltaT/(deltaT + deltaT0);
            Coo.value() = deltaT*deltaT/(deltaT0*(deltaT + deltaT0));
            Co.value() = Cn.value() + Coo.value();

            Cn.value() /= deltaT;
            Co.value() /= deltaT;
            Coo.value() /= deltaT;

            phi() +=
                rAUf_
               *(
                   Co*phi().oldTime()
                 - Coo*phi().oldTime().oldTime()
                );
        }
    }
}
else
{
    phi() +=
        rAUf_
       *(
           fvc::interpolate(ddtUEqn.H(), "interpolate(U)") & mesh.Sf()
        );
}

forAll(phi().boundaryField(), patchI)
{
    if
    (
        U().boundaryField()[patchI].fixesValue()
    )
    {
        //Info << mesh.boundary()[patchI].name() << endl;

        phi().boundaryField()[patchI] =
        (
            U().boundaryField()[patchI]
          & mesh.Sf().boundaryField()[patchI]
        );
    }
    else if
    (
        isA<slipFvPatchVectorField>(U().boundaryField()[patchI])
    )
    {
        //Info << mesh.boundary()[patchI].name() << endl;

        phi().boundaryField()[patchI] = 0;
    }
}
