# Prestress case for `pulmonaryArteries`

This is the solid-only prestress case of the `pulmonaryArteries` fluid-solid
interaction tutorial; see `../README.md` for the full description.

It computes the prestress `Fpre` of the imaged, loaded arterial wall with
repeated passes of `solids4Foam` and `updatePrestress`, and copies the result
to `../0/solid/Fpre` for the fluid-solid interaction run. It is normally run by
`../Allrun`. To rerun it on its own, after `../Allmesh` has copied the solid
mesh here:

```bash
./Allclean
./Allrun parallel
```

The wall material is read from `../constant/solid/mechanicalProperties`, and
the number of passes is set by `nPasses` in `Allrun`. Each pass appends its
convergence measures to `prestressHistory`.
