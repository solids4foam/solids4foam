# Phase field fracture material-point checks

`Test-phaseFieldFractureMechanicalConstitutiveLaw` drives
`linearElasticMisesPlasticPhaseFieldFracture` directly, at one point, with the
phase field supplied as the input the solid model would solve for. It checks
that:

1. With the phase field zero, the law is `linearElasticMisesPlastic`, step
   for step.
2. The stress is degraded by `g(d)`, plastic flow included, except for its
   volumetric part in compression; and a crack keeps the residual stiffness.
3. On a path that yields and unloads, the history of the tensile elastic
   energy, the plastic work, and the source and coefficient of the phase
   field equation are those recomputed here; and a residual stiffness `k`
   scales the source by `1 - k`.
4. For the brittle model, with the phase field iterated to the solution of
   its equation for a uniform field, `A/(1 + A)`, the uniaxial strain
   response is the analytical one, peaking at `A = 1/3` at 9/16 of the
   intact stress.
5. Under `logStrain`, with a phase field of 0.3, on a path whose principal
   directions stay fixed, the Kirchhoff stress equals the small-strain law's
   at the logarithmic strain.
6. Inconsistent entries are refused.
7. The law declares the phase field, with its coefficient.

The laws read `constant/linearHardening`, `constant/nonlinearHardening` and
`constant/elastic`, whose yield stress is out of reach, for the brittle
checks.

Run it with:

```bash
./regressionTest.sh
```

`tutorials/Alltest-regression` runs it with the tutorial regression tests.
