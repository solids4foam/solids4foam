// Calculate ADt
{
    scalarField& ADtI = ADt_.internalField();
    ADtI = 0;

    // Add unsteady term contribution
    if (!steadyState)
    {
        ADtI =
            (1.0-backward.value())*sqr(coeff_.value())*
            rho().internalField()*rDeltaT.value()*
            rDeltaT.value()*mesh().V().field()

          + backward.value()*rho().internalField()*
            rDeltaT.value()*Cn.value()*
            Cn.value()*rDeltaT.value()*mesh().V().field();
    }

    ADtI /= mesh().V().field();
    ADt_.correctBoundaryConditions();
}
