# Lemaitre damage material-point checks

`Test-LemaitreDamageMechanicalConstitutiveLaw` drives
`linearElasticMisesPlasticLemaitreDamage` directly, at one point, over strain
paths of many steps, so that its damage and plastic history can be read back.
It checks that:

1. With the damage threshold out of reach, the law is
   `linearElasticMisesPlastic`, step for step.
2. For each `lemaitreType`, on a damaging path, the returned state satisfies
   the discrete yield condition and damage update, recomputed from the
   stress with the energy release rate in the triaxiality form of Whelan et
   al. (2025), equation 52.
3. Under uniaxial stress, held by solving for the lateral strain at each
   step, the damage converges at first order to the analytical solution of
   Doghri for linear hardening and `b = 1`.
4. Crack closure has no effect when every principal stress is tensile and
   reduces the damage when one is compressive.
5. The triaxiality cut-off keeps the damage at zero under uniaxial
   compression of the strain.
6. A point whose damage passes `Dc` fails at 0.99 and stays failed.
7. The damage converges at first order as the step is refined.
8. Under `logStrain`, on a path whose principal directions stay fixed, the
   Kirchhoff stress equals the small-strain law's at the logarithmic strain.
9. Inconsistent entries are refused.

The laws read `constant/linearHardening`, with two rows, so that the return
mapping is in closed form and the equations of checks 2 and 3 hold to
round-off, and `constant/nonlinearHardening`, which runs the Newton loop.

Run it with:

```bash
./regressionTest.sh
```

`tutorials/Alltest-regression` runs it with the tutorial regression tests.
