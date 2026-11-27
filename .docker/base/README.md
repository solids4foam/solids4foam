# Base Docker images

The solids4foam CI workflows and release images (`.docker/Dockerfile.*`) are
built on base images that hold OpenFOAM and, optionally, PETSc and Gmsh. This
folder holds the Dockerfiles of the base images that are built on Ubuntu 24.04.

- `Dockerfile.openfoamv2512-ubuntu24.04`: OpenFOAM-v2512 on Ubuntu 24.04,
  published as `philippic/openfoam-v2512:ubuntu-24.04`.
- `Dockerfile.petsc`: adds PETSc with MUMPS and HYPRE, without OpenMP and with
  the system OpenBLAS, published with the tag suffix
  `-petsc-no-openmp-system-blas`.
- `Dockerfile.gmsh`: adds Gmsh, for the regression tests, published with the
  tag suffix `-petsc-no-openmp-system-blas-gmsh`.

Each layer is built on the previous one, given as the `BASE` build argument:

```bash
cd .docker/base
docker build -t philippic/openfoam-v2512:ubuntu-24.04 \
    -f Dockerfile.openfoamv2512-ubuntu24.04 .
docker build --build-arg BASE=philippic/openfoam-v2512:ubuntu-24.04 \
    -t philippic/openfoam-v2512:ubuntu-24.04-petsc-no-openmp-system-blas \
    -f Dockerfile.petsc .
BASE=philippic/openfoam-v2512:ubuntu-24.04-petsc-no-openmp-system-blas
docker build --build-arg BASE="${BASE}" \
    -t "${BASE}-gmsh" \
    -f Dockerfile.gmsh .
```

The `docker-base-images.yml` workflow does the same. It runs on pull requests
that change this folder, to check that the images still build, and pushes them
to Docker Hub when it is run manually with `push` ticked.

The base images of the older OpenFOAM versions (v2412, v2312, OpenFOAM-9 and
foam-extend-4.1) are still on Ubuntu 22.04 or 20.04 and are not yet built from
this folder.
