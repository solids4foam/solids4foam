# immersedSurfaceTraction quadrature checks

`Test-immersedSurfaceTraction` checks the quadratures of the surface traction
of the `immersedBoundaryForce` option (`tractionQuadrature clustered`, `none`
and `greedy`, see `src/immersedBoundary/README.md`) without a flow solution.
The pressure is set to constant and linear fields, whose least squares wall
fits are exact, and the velocity to that of the body, so that the force is the
quadrature of `-p n` alone. For a closed surface enclosing the volume `V`,

```text
integral(n dA) = 0,   integral(x_j n_i dA) = V delta_ij,
```

so the force of `p = 1` is zero, that of `p = x_j - c_j` is `-V e_j`, and a
constant added to the pressure does not change the force.

The surfaces are written by the application, in three facet orders (as
generated, reversed, and permuted with rotated vertices), and evaluated
static, rotated, translated and bent (`quadraticBend`):

- `twoD` and `twoDGraded`: a 128-sided cylinder, a rotated ellipse, a folded
  thin strip and two stacked thin plates 1.5 mm apart, extruded through a
  two-dimensional mesh, uniform (`twoD`) or graded in `x` (`twoDGraded`), on
  which the motions move the bodies into cells of a different width;
- `threeD`: a rotated ellipsoid in a three-dimensional mesh.

The checks (PASS or FAIL lines in `log.Test-immersedSurfaceTraction`) are:

1. `none` and `clustered`: the force of `p = 1` is zero to `1e-12` of the
   area;
2. `clustered`: a constant added to the pressure does not change the force;
3. `none`: the first moments of the area are `V I`, and `V` is the enclosed
   volume for the rigid motions;
4. `clustered`: the first moments and the torque of `p = x - c` are those of
   `none` to `1e-2` (the second-order error of the clusters);
5. `clustered`: the total and face areas are those of `none`;
6. `clustered`: the traction at each point is the physical traction,
   `|t| = p(x)` for `p = 1 + x - c_x`, while the force integrates the
   pressure with the vector area of each cluster (check 1);
7. `clustered`, two dimensions: at most 0.6 times the points of `none`;
8. `none` and `clustered`: the forces, torques and numbers of points do not
   depend on the facet order;
9. `clustered`: the clusters made after the motion, as on a restart, are
   those made in the reference configuration, also on the graded mesh.

The `greedy` quadrature is reported (INFO lines) only: its forces depend on
the facet order.

Run it with:

```bash
./regressionTest.sh
```

`tutorials/Alltest-regression` runs it with the tutorial regression tests. The
script requires the expected number of checks on each mesh, so that a check
that stopped, or a surface or motion dropped from the application, fails the
test. It is skipped with foam-extend and the OpenFOAM Foundation versions,
where the immersed boundary library is not built.
