# Build helper scripts

This directory contains shared helpers used by the solids4foam build system.

## `wmake-options`

`wmake-options` is a makefile fragment included by the solids4foam
`Make/options` files. It detects whether the active environment is
OpenFOAM.com, OpenFOAM.org, or foam-extend and defines the corresponding
preprocessor and library settings. It also selects the OpenFOAM module output
directories when available, falling back to `FOAM_USER_APPBIN` and
`FOAM_USER_LIBBIN`.

The file is consumed automatically during `wmake`; it is not run directly.
New `Make/options` files should include it using the repository-relative
`SOLIDS4FOAM_ROOT`, following the existing build files.

## `check-petsc.sh`

`check-petsc.sh` validates PETSc before the solids4foam models library is
built. If `PETSC_DIR` is unset, it reports that PETSc support is disabled. If
PETSc is enabled, it reads `petscversion.h` and `petscconf.h`, reports the
detected version, and checks that PETSc was configured with MUMPS and HYPRE.
`PETSC_ARCH` is also honoured for non-installed PETSc builds.

The check runs automatically from `src/solids4FoamModels/Allwmake`. It can also
be run manually from the repository root:

```bash
./etc/check-petsc.sh
```

A missing MUMPS or HYPRE configuration causes a non-zero exit status with
reconfiguration guidance. To intentionally build against an incomplete PETSc
installation, set the explicit override:

```bash
export S4F_ALLOW_INCOMPLETE_PETSC=1
```

The script's standard output is a PETSc configuration fingerprint used by the
incremental build logic; informational messages are written to standard error.
