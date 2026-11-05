# immersedBoundary

Hybrid fictitious domain-immersed boundary (HFDIB) method for bodies with a
prescribed motion, or driven by a solids4foam fluid-solid interaction,
immersed in an incompressible flow, as the finite volume option
`immersedBoundaryForce`. It is currently only available for OpenFOAM.com,
where `libimmersedBoundary` is built with solids4foam.

The option works with the solids4foam `pimpleFluid` fluid model and with the
standard `pimpleFoam` solver. See `tutorials/fluids/immersedBoundary`.

The bodies are closed surfaces with a thickness, represented by the cells of
the fluid mesh whose centre is inside them: **the mesh needs about three cells
or more across the thickness of a body**, everywhere and at all times. Thinner
bodies, or thinner parts of a body, are not blocked reliably and let the flow
through. Zero-thickness (open) surfaces are not supported. For thin moving
bodies, such as the leaflets of a heart valve, refine the mesh around the
region they sweep (see Mesh refinement below).

## Usage

In `system/controlDict`:

```c++
libs (immersedBoundary);
```

In `constant/fvOptions` (for a fluid-only solids4foam case), or
`constant/fluid/fvOptions` (for the fluid region of a fluid-solid interaction
case):

```c++
immersedBoundary
{
    type            immersedBoundaryForce;

    // Forcing method, optional, default cutLink
    method          cutLink;

    bodies
    {
        cylinder
        {
            // Closed surface, relative to constant/triSurface
            surface     "cylinder.stl";

            // Prescribed motion, optional, default static
            motion
            {
                type        sinusoidalTranslation;
                amplitude   0.25;
                period      4;
                direction   (1 0 0);
            }
        }
    }
}
```

The complete list of entries is in the header of
`fvOptions/immersedBoundaryForce/immersedBoundaryForce.H`. The body motions
are `static` (the default), `sinusoidalTranslation`, `uniformTranslation`
(constant velocity) and `solidBodyRotation` (constant angular velocity about
an axis), the deforming `quadraticBend` and `customProfileBend` (bending
of a beam about its base), and the heart valve motions `valveSliceAxis` (the
collapse of the leaflets towards the valve axis, with the vertices on the
annulus held fixed) and `valveMorph` (the interpolation between a closed and
an open valve surface with the same vertices), with the periodic time laws of
`valveTimeLaw`, in `immersedBodyMotion`. A deforming motion gives
the position and velocity of each vertex of the surface, and the velocity of
the body at a point is interpolated from the vertices of the nearest triangle.
The `fsiDriven` motion is set by a fluid-solid interface (see below).

## Fluid-solid interaction

A body whose motion is `fsiDriven` is driven by a solids4foam fluid-solid
interface: in `constant/fsiProperties`, the solid interface patch is paired
with `fluidPatch none;` and an `immersedInterfaces` entry names the body
(see `src/solids4FoamModels/fluidSolidInterfaces/README.md`). The body needs
no `surface` entry: its surface is built by the interface from the solid
interface patch and, optionally, the closure patches of the solid (e.g. a
clamped root), whose union must be a closed surface; in a two-dimensional
fluid mesh the surface is extended through the mesh and capped in the empty
direction. In each coupling iteration the interface sets the vertex
positions and velocities of the body, and the bodies are re-positioned in
the mesh (occupancy, links, apertures) when the momentum equation is next
assembled, keeping the forcing, momentum and apertures of the previous time
step. The traction returned to the solid is the surface traction (`cutLink`
only) averaged over the quadrature points of each interface face, with every
candidate point (see "Surface traction" below); the log
compares its total with the momentum exchange between the fluid and the
whole body, which also includes the inertia of the fluid inside the body and
the closure patches, so the two agree only for a slowly accelerating body.
The option is found by the interface through the header-only
`fsiImmersedBoundary` interface of solids4FoamModels, so that
`libimmersedBoundary` does not depend on solids4FoamModels. See the
`immersedHronTurekFsi2` tutorial.

## Method

Each body is a closed surface that covers some cells of the fluid mesh with an
occupancy (solid volume fraction) `lambda`. By default, `lambda` is found from
the signed distance of the cell centres to the surface, relative to the width
of the cells normal to the surface, so that it varies continuously as the
body moves (`occupancy signedDistance`). Alternatively, as in openHFDIB-DEM,
it is half the fraction of the cell vertices inside the surface, plus a half
if the cell centre is inside the surface (`occupancy vertexFraction`).

Four forcing methods are available:

- `penalty`: the implicit volume penalisation
  `kappa*(Ui - U)`, where `Ui` is the velocity of the body, is added to the
  momentum equation, with `-kappa*U` in the matrix. The penalisation rate is
  `K/deltaT` in the fully covered cells, where `K` is `penaltyCoeff` (default
  1e3). In the partially covered cells, `kappa*deltaT = lambda/(1 - lambda)`
  (`weighting volumeFraction`, the default), so that the penalised velocity is
  `lambda*Ui + (1 - lambda)*U`; or `kappa = K*lambda/deltaT` (`weighting
  occupancy`), which blocks every covered cell. With the `volumeFraction`
  weighting, the rate in the partially covered cells is also limited to
  `lambda/(1 - lambda)*C*(nu/w^2 + max(|Ui|, |U|)/w)`, where `w` is the
  cell width and `C` is `surfaceRateCoeff` (default 3). The maximum preserves
  the body-motion scale and includes fluid advection for static bodies.
  This numerical velocity scale is relative to the fixed mesh, and the cap
  has no explicit time-step dependence. It is recomputed from the current
  fluid velocity at each momentum assembly; timestep and iteration
  independence must therefore be checked for each case. Local fluid speed
  is itself damped by the penalty, so this is not an estimate of incident
  flow speed or a guarantee of a Reynolds-independent boundary location.
  With `surfaceRateCoeff 0`, the partially covered cells are relaxed towards
  `Ui` once per time step.
- `ghostCell`: a sharp interface variant of `penalty`, for static bodies.
  The cells whose centre is inside a body are penalised with the rate
  `K/deltaT`. In the cells next to the fluid (the ghost cells), the target
  velocity is extrapolated linearly along the surface normal, through the
  body velocity at the nearest surface point and the velocity at an image
  point in the fluid, `imageDistance` (default 1.5) cell widths beyond the
  surface. The image point velocity is interpolated in its donor cell
  (`imageInterpolation cellPoint`, the default) or is the donor cell value
  (`imageInterpolation cell`). The donor may be on another processor. The
  fluid velocity then equals the body velocity on the surface, so the
  effective wall is not displaced into the body as it is with `penalty`.
  The results do not depend on `K`, the time step, the number of pressure
  correctors or the number of processors. For moving bodies, the cells that
  the body leaves (fresh cells) disturb the force: `cutLink` or `penalty`
  is recommended.
- `cutLink` (the default): a sharp interface penalisation for static and
  moving bodies.
  The cells whose centre is inside a body are penalised towards the body
  velocity, with the rate `min(K/deltaT, C*(nu/w^2 + |Ub|/w))`, where `C` is
  `pinnedRateCoeff` (default 100). Each fluid cell next to them is penalised
  towards the body velocity at the surface with the rate
  `sum(nu*|Sf|*deltaCoeff*(|d|/phi - 1))/V` over its faces with a penalised
  neighbour, where `|d|` is the distance between the cell centres and `phi`
  the distance from the fluid cell centre to the surface along the line
  between them, from the intersection with the surface; with the explicit
  correction of the flux to the penalised cell (`linkCorrection`, default
  yes), the viscous flux through the face is that of a linear profile
  through the body velocity on the surface (Shortley-Weller), so that the
  wall is on the surface to second order. The rates do not depend on the
  time step and vary continuously as the body moves. With
  `apertureCoupling` (default yes), the option also registers the fluid
  fraction of the area of the faces cut by the bodies, from the geometric
  cutting of the faces (`cutFaceIso`), and the flux of the body velocity
  through their solid part; the `pimpleFluid` fluid model uses them in its
  pressure equation, so that it changes continuously as the body moves, and
  does not impose continuity in the cells entirely inside a body, where the
  flux of the velocity of a deforming body (for example a bending beam,
  whose sections stretch on one side and shorten on the other) need not be
  zero, and would otherwise be a source of fluid that leaks into the flow
  (this has no effect for rigid bodies). The face
  apertures are sampled at the mesh points, which miss a body thinner than
  the point spacing for a part of the time; an experimental switch for such
  bodies, `thinBodyBlocking` (default `off`), is described in
  [Experimental: thin bodies](#experimental-thin-bodies). The
  force is the momentum exchange (`forceEstimator momentumExchange`, the
  default), which includes the inertia of the fluid inside the body; the
  alternative surface traction (`forceEstimator surfaceTraction`) is from
  quadratic least squares fits of the pressure and velocity at quadrature
  points moving with the surface. A third estimate (`forceEstimator forcing`) is
  the forcing applied in these cells plus the inertia of the fluid inside
  the body taken as rigid. The estimates that are not selected are written
  after the inertia: the momentum exchange, unless it is selected, then the
  surface traction, with `writeSurfaceTraction yes`, unless it is selected.
  The momentum exchange assumes a laminar flow of constant viscosity on a
  static mesh.
- `incremental`: the direct forcing of the openHFDIB-DEM `pimpleHFDIBFoam`
  solver. An explicit forcing `f` is added to the momentum equation. After
  each pressure corrector, `f` is increased by `couplingCoeff*(Ui - U)/deltaT`
  in the covered cells, and the updated forcing is used by the next pressure
  corrector; at the start of each time step, `f` is multiplied by the new
  occupancy. The results depend on `couplingCoeff`, the time step and the
  number of pressure correctors. This method is kept for comparison: with
  `occupancy vertexFraction` it reproduces `pimpleHFDIBFoam`.
  `couplingCoeff` controls how much of the velocity error is corrected after
  each pressure corrector. The default `0.8` gives strong enforcement without
  the oscillation that full correction can produce. Reduce it when forcing
  or pressure corrections oscillate or diverge; a smaller value enforces
  body velocity more slowly and can require more pressure correctors or a
  smaller timestep for the same accuracy.

The force `-rho*sum(f*V)` and torque on each body, where `f` is the forcing
(an acceleration) exerted on the fluid, are written every time step to
`postProcessing/<option name>/<start time>/<body>.dat`, together with the
inertia of the fluid inside the body, `rho*d/dt(sum(lambda*U*V))`. For a
moving body, the hydrodynamic force is the sum of the two (Uhlmann, 2005).

The fields `<option name>:lambda`, `<option name>:Ui`, `<option name>:kappa`
and `<option name>:f`
are written at write times; with the `incremental` method,
`<option name>:f` is read on restart.

### Accuracy

Before adding fluid speed to the surface-rate cap, the two tutorials
with 10, 20 and 40 cells across the cylinder gave the following results. The
static cylinder drag coefficient (reference 5.57-5.59), and the root mean
square difference of the oscillating cylinder drag coefficient, including the
inertia of the fluid inside the cylinder, from a moving body-fitted mesh
solution and from Wan and Turek (2007), over 0.25 < t < 7.5 s (where the root
mean square of the drag coefficient is 2.05), are:

| Settings | Static `Cd` | Oscillating: body-fitted | Wan and Turek |
| -------- | ----------- | ------------------------ | ------------- |
| A | 5.24, 5.38, 5.48 | 0.11, 0.05, 0.02 | 0.13, 0.09, 0.08 |
| B | 5.45, 5.61, 5.61 | 2.00, 0.44, 0.29 | 1.42, 0.33, 0.20 |
| C | 6.04, 6.42, 6.60 | 1.73, 0.60, 0.53 | 1.32, 0.52, 0.48 |
| D | 5.42, 5.51, 5.55 | 0.80, 0.30, - | 0.59, 0.27, - |
| E | 5.61, 5.58, 5.58 | 0.19, 0.07, 0.06 | 0.17, 0.07, 0.06 |
| F | 5.36, 5.46, 5.48 | 0.11, 0.04, 0.02 | 0.12, 0.08, 0.08 |

where the settings are:

- A: `penalty`, `volumeFraction`, `signedDistance`, with the previous
  body-speed-only surface-rate cap;
- B: `penalty` with `weighting occupancy` and `occupancy vertexFraction`;
- C: `incremental` with `occupancy vertexFraction`;
- D: `ghostCell`;
- E: `cutLink` (the setting of the static tutorial), with the force from the
  momentum exchange (the default);
- F: `cutLink`, with the force from the surface traction.

For an immersed plane Poiseuille flow between walls that are not aligned with
the cell faces, with 20, 40 and 80 cells across the channel, the flow rate
differs from the exact solution by 20%, 11% and 7% with setting A, by
0.9%, 0.03% and 0.02% with setting D, and by 0.25%, 0.23% and 0.05% with
setting E. For a cylinder translating at constant velocity through the mesh
(the `translatingCylinderInChannel` tutorial), `cutLink` gives the drag of the
static cylinder at the same position to within 0.2% on all the meshes, with
fluctuations of 1%, 0.4% and 0.3% as the cylinder crosses the cells, whereas
the `penalty` drag is low by 5%, 3% and 1.5%. For the immersed Stokes layer
(the `oscillatingWallStokesLayer` tutorial), the `cutLink` velocity converges
at second order. For a cylinder rotating inside an immersed annulus (the
`immersedTaylorCouette` tutorial), the `cutLink` torque differs from the exact
torque by 0.8%, 0.3% and 0.2%, and the `penalty` torque by -25%, -16% and
-10%. For a beam bending about its base (the `bendingBeamInChannel`
tutorial), the root mean square difference of the `cutLink` drag from a
body-fitted solution with a deforming mesh is 7.4%, 4.6% and 3.5% of the
root mean square drag (momentum exchange), or 8.6%, 4.9% and 3.5% (surface
traction), with 3, 6 and 11 cells across the beam. The momentum exchange is
the force that the fluid receives, and is the more accurate for static and
steadily moving bodies; the surface traction is the more accurate for bodies
accelerating normal to their surface (oscillating cylinder: 0.11, 0.04 and
0.02); the `forcing` estimate is as
accurate as the momentum exchange for static and steadily moving bodies, and
gives the wall shear stress of the Stokes layer to 0.1% with 8 cells across
the layer, but over-predicts the oscillating cylinder drag amplitude by 3%.

The cylinder forces converge at about first order in the cell size for all the
methods, except the static `cutLink` drag, which is within 0.6% of the
reference from 10 cells across the cylinder (the ghost cell method was not
run with 40 cells across the moving cylinder). The differences from Wan and
Turek (2007) stop decreasing at about 0.08, the difference between the
body-fitted mesh solution and Wan and Turek (2007), whose coefficients,
digitised from a figure of Erzincanli and Sahin (2013), lead the converged
solutions by about 0.018 s (see the `oscillatingCylinderInChannel`
tutorial).

### Fluid-speed surface-rate check

With the maximum of body and fluid speed in the cap, OpenFOAM v2512 checks
with fixed timesteps gave static `Cd = 5.3027, 5.3944` at 10 and 20 cells
across the cylinder (at 10 s). On the coarse mesh, reducing the timestep
from 0.01 to 0.0025 s gave `Cd = 5.3029`, and three instead of one PIMPLE
outer correctors gave `Cd = 5.3027`.

For the oscillating cylinder, full 8 s runs at a fixed timestep of 0.0025 s
gave RMS drag differences from the body-fitted reference of `0.1039, 0.0454`
at 10 and 20 cells across the cylinder, over `0.25 < t < 7.5 s`, including
the inertia of the internal fluid. On the coarse mesh, a timestep of
0.000625 s gave `0.1363`, while three outer correctors gave `0.1018`.
These checks show measurable timestep sensitivity in the moving case;
absence of explicit timestep dependence in the cap does not imply timestep
independence of the solution. The finest mesh was not reverified with this
cap; the three-level table above records the earlier body-speed-only cap.

At one-tenth viscosity, the RMS fluid speed in the partially covered cells
at 10 s fell from 0.131 to 0.065 m/s compared with the body-speed-only cap.
This is a damping diagnostic, not a boundary-location or drag-accuracy test:
that flow is unsteady and has a different Reynolds number from the supplied
static reference. In a short zero-viscosity diagnostic, all 36 partially
covered cells had positive penalty rates, compared with zero for the old
cap. This does not establish inviscid no-slip accuracy.

### Mesh refinement

The `immersedBodyRefinementCells` utility writes the cell set of the cells
within a given distance of the immersed bodies of the `immersedBoundaryForce`
options, over their motion between two times, for `refineMesh`. The surfaces
are sampled wherever a vertex has moved by half the distance since the
previous sample, so that the bands around the samples overlap. Two levels of
refinement, halving the band at each level:

```bash
immersedBodyRefinementCells -distance 0.002 -startTime 0 -endTime 0.8
refineMesh -overwrite
immersedBodyRefinementCells -distance 0.001 -startTime 0 -endTime 0.8
refineMesh -overwrite
```

with `system/refineMeshDict` refining the cell set `refineCells` in the three
directions (see the `heartValveInDuct` tutorial, `REFINE_LEVELS=1 ./Allrun`).
The option `-bodies "(name ...)"` restricts the refinement to some bodies.

For the valves of solids4foam/cardiacFoam#20 (leaflets about 0.9 mm thick) in
a mesh of 1.3 mm cells, no cell centre is inside the open leaflets for most of
the cycle, and the leaflets are not seen; after two levels of refinement
around the swept region (0.34 mm cells, from 213 000 to 779 000 cells), 4 000
to 19 000 cells are inside the mitral leaflets and 1 700 to 12 000 inside the
aortic leaflets at all times. The `valveSliceAxis` motion shrinks the
thickness of the leaflets with their radius, so the free edge of a closed
valve with `Fmax` close to 1 can remain thinner than the refined cells.

## Surface traction

With the `cutLink` method, `forceEstimator surfaceTraction` (the force of the
body) or `writeSurfaceTraction yes` (columns 11-13 of the force file) evaluates
the traction exerted by the fluid on the surface (`immersedSurfaceTraction`).
At points on the surface, the pressure and the velocity relative to the body
are fitted by weighted quadratic least squares over the fluid cells within
three cell widths, and the kinematic traction is

```text
t = -p_w n + nu (I - n n) . (grad(U - U_b) . n)
```

with the wall pressure `p_w` and the unit normal `n` (with the surface
gradient terms of the body velocity for a deforming body); the force is
`rho` times its quadrature over the surface.

### Quadrature

The candidate points are made on the surface triangles at a spacing of about
the mean cell width: in two dimensions along the intersection of the surface
with the mid-plane of the mesh, each representing the thickness of the mesh,
and in three dimensions at the centroids of a uniform subdivision of each
triangle. Each candidate is defined by its triangle and barycentric
coordinates, so that it follows the motion or deformation of the surface.
`tractionQuadrature` selects how they are used:

- `clustered` (the default): one point per cluster of candidates. It does not
  depend on the order of the triangles, the force of a constant pressure is
  that of every candidate, and the moment error is of second order;
- `none`: every candidate, the reference, with the most points;
- `greedy` (legacy, the default of earlier versions): the candidates thinned
  to 0.7 mean cell widths, in the order of the triangles. It depends on that
  order, and the moment error is of first order.

```c++
// cutLink, with the surface traction only
tractionQuadrature   clustered;    // default
tractionPointSpacing 0.5;          // clustered only, default 0.5
```

Setting either entry without the surface traction, or with another method,
is an error, as is `tractionPointSpacing` without `clustered`. Bodies driven
by a fluid-solid interface always use `none` (see below).

The greedy thinning keeps a candidate if no kept candidate is closer than
the spacing, in the order of the triangles, and gives the area of each
removed candidate to the nearest kept one, which keeps its own position and
normal. The kept points are therefore offset from the area they represent,
in the direction of the facet order: a pressure varying along the drag
direction produces a lift whose sign depends on the order of the facets
(for the moving cylinder of `oscillatingCylinderInChannel`, a lift 140%,
149% and 140% off the body-fitted reference with 10, 20 and 40 cells across
the cylinder, reversed with the facets of the STL), and on thin bodies the
two faces can be merged. It is kept, as
`tractionQuadrature greedy`, to reproduce earlier results: it gives the force
files of earlier versions unchanged.

### Clustered quadrature

The candidates are made in the reference configuration of the surface and
grouped by

- the cell of a uniform grid of spacing `s` in the solved directions, with
  its origin at the origin of the coordinates, where `s` is
  `tractionPointSpacing` times the smallest width of the cells containing the
  candidates in the reference configuration;
- the dominant signed axis of their normal, so that the two sides of a thin
  body are never merged;
- whether they are inside the mesh (in the reference configuration);

and each such group is split into the connected components of the triangles
of its candidates, so that a cluster never joins surface regions that are not
connected locally (the layers of a fold, nearby bodies). The clusters are
numbered by these geometric keys, and the triangle of each cluster's point
is chosen by quantised distances and positions, so that nothing depends on
the order of the triangles or on rounding.

At each time step, each cluster k is represented by one point, its current
area centroid projected on the plane of one of its triangles, with the area
`A_k = sum_i A_i`, the vector area `a_k = sum_i A_i n_i` of its candidates on
their current triangles, and the unit normal `n_k = a_k/|a_k|`. The force is

```text
F = rho sum_k (-p_k a_k + A_k tau_k)
```

with the wall pressure `p_k` and the viscous traction `tau_k` of the fits at
the point: the pressure is integrated with the vector area, and the viscous
traction with the area. Then

- the force of a constant pressure is `-p sum_k a_k = -p sum_i A_i n_i`, that
  of the candidates: zero for a closed surface in the mesh, so that the
  pressure level, which is arbitrary in an incompressible flow, gives no
  force. With `A_k n_k` instead of `a_k`, the pressure force of a curved
  cluster would be too large by `p_k (A_k - |a_k|)` along `n_k` (since
  `|a_k| <= A_k`), and a constant pressure would give a force;
- the total area is that of the candidates;
- the first moment of the area of each cluster about its point has no
  tangential part, so that the quadrature error of a smooth traction is of
  second order in the cluster size (in the manufactured test, the first
  moments differ from those of every candidate by at most 8e-4 of the
  volume).

The traction at the points (`immersedSurfaceTraction::traction`) is the
physical traction `-p_k n_k + tau_k`: the vector area is a quadrature weight
of the force only. `writeSurfaceTraction` writes the integrated force, not
values at the points. `faceValues` gives each candidate's triangle the value
of its cluster with the candidate's area, so that the face values of the
traction integrate the pressure with the area rather than the vector area of
each cluster (a difference of second order in the cluster size); its only
user is the fluid-solid interface, which uses `none`.

### Invariance and its limits

- Order of the triangles: the clusters and their points do not depend on it
  (manufactured test: forces equal to 3e-14 for three facet orders, 2.5e-11
  on the graded mesh, where the fits are solved with a larger rounding error;
  moving cylinder: original and reversed STL equal to 2e-11).
- Decomposition: every processor holds the whole surface and the cell widths
  are reduced over the processors, so that the clusters are the same on every
  processor (moving cylinder: 1 and 2 processors as `none`).
- Rigid motion and deformation: the clusters are made once, in the reference
  configuration, and move with the surface; their points, areas, vector
  areas and normals are found again at each time step.
- Restart: the clusters are a function of the reference surface, the mesh and
  `tractionPointSpacing` only (the spacing and the inside or outside of the
  mesh are those of the reference configuration), so that a restart, or
  clusters made at any other time, give the same clusters on the same mesh.
  The manufactured test checks this on a mesh graded so that the motions move
  the bodies into cells of a different width.
- Local refinement: the whole surface gets the spacing of the finest cells
  containing it in the reference configuration, and the spacing is not
  changed as the body moves, so that a body moving into finer cells keeps
  clusters that are coarser relative to the local cells, and one moving into
  coarser cells more points than needed. On a mesh refined around the swept
  region (`immersedBodyRefinementCells`), the reference configuration should
  be in the refined region.
- Changing meshes: the clusters are made with the mesh at their creation
  (the first force evaluation of a run); a mesh whose cells change, e.g. by
  adaptive refinement, or a restart on a different mesh, can give different
  clusters.

### Cost

The clustered quadrature never has more points than `none`. With the default
spacing (half the smallest cell width), for the moving cylinder of
`oscillatingCylinderInChannel` (fixed time steps, 8 s):

| Cells across D | Points (greedy, clustered, none) | CPU hours (same order) |
| -------------- | -------------------------------- | ---------------------- |
| 10 | 32, 64, 256 | 0.023, 0.029, 0.054 |
| 20 | 64, 136, 256 | 0.238, 0.304, 0.426 |
| 40 | 128, 232, 256 | 3.39, 3.76, 4.40 |

that is 11-28% more than `greedy`. In three dimensions the default spacing
gives more points relative to `greedy` (798 against 193 for the ellipsoid of
the manufactured test, 1280 for `none`); a larger `tractionPointSpacing`
reduces them. Since the spacing follows the smallest cell width at the
surface, on locally refined meshes the whole surface gets the fine spacing.

### Bodies driven by a fluid-solid interface

A body driven by a fluid-solid interface (`fsiDriven`) always uses `none`.
The traction is returned to the solid face by face: `faceValues` averages the
values of the points of each surface triangle, and the interface averages the
triangles of each interface face. This needs points on every triangle, at
their own positions, with their own areas. The greedy quadrature gives the
area of the removed points to kept points of other triangles, possibly of
other interface faces, and can leave faces without points. The clustered
quadrature keeps the area of every triangle, but gives each triangle the
traction of its cluster, along the cluster normal: the force returned to the
solid then differs from that of the cluster by `p_k (A_k - |a_k|) n_k`, unless
the pressure is applied with the normal of each triangle (`-p_k n_i A_i`, which
sums to `-p_k a_k`). This transfer is not implemented, and the
fluid-solid interface with the immersed boundary (`immersedHronTurekFsi2`,
which requires PETSc) has not been run with the clustered quadrature, so its
use there is left as a separate development.

### Verification

`applications/test/immersedSurfaceTraction` (see its `case/README.md`)
checks the quadratures without a flow solution, against the exact integrals
of constant and linear pressures on closed surfaces: 1294 checks on a
uniform and a graded two-dimensional mesh and on a three-dimensional mesh,
for three facet orders, static, rotated, translated and bent. Run it with

```bash
applications/test/immersedSurfaceTraction/case/regressionTest.sh
```

For the moving cylinder (`oscillatingCylinderInChannel`), the root mean
square differences of the surface-traction drag and lift from a moving
body-fitted solution over 0.25 < t < 7.5 s, relative to the root mean squares
of the body-fitted coefficients, are:

| Cells across D | Drag (greedy, clustered, none) | Lift (same order) |
| -------------- | ------------------------------ | ----------------- |
| 10 | 5.36, 5.40, 5.38% | 140, 13.2, 13.2% |
| 20 | 1.83, 1.86, 1.87% | 149, 6.7, 6.7% |
| 40 | 1.10, 1.08, 1.10% | 140, 4.9, 4.8% |

The body-fitted solution has 901 120 cells. The `greedy` values come from
runs of earlier versions; with 80 cells across the cylinder, the greedy
thinning keeps all 256 candidate points and its lift differs by 4.6%, a
property of that mesh and surface, not a convergence of the greedy
quadrature. The `clustered` values come from 8 s runs made while the
quadrature was developed. This version was compared with those runs with
10 cells across the cylinder, and with 20 and 40 cells only over the first
0.5 s; no clustered run was made with 80 cells.

With this version, `clustered` (the default) against `greedy`:

- `bendingBeamInChannel`: the root mean square difference of the surface
  traction drag from the body-fitted solution over 1 < t < 8 s, relative to
  that of the body-fitted drag, is 7.7%, 4.6% and 3.1% for `MESH_LEVEL`
  1-3, against 8.6%, 4.9% and 3.1% (the momentum exchange gives 7.4%, 4.6%
  and 3.3%);
- `staticCylinderInChannel`: the surface-traction lift at 10, 20 and 40 cells
  across the cylinder is 0.0026, 0.0078 and 0.0084 (`none`: 0.0083 at 40
  cells), against -0.060, -0.063 and -0.048 (reference 0.0104-0.0110). The
  drag is 5.35, 5.46 and 5.49 (`none`: 5.48 at 40 cells), against 5.35, 5.46
  and 5.50: the greedy drag at 40 cells is closer to the reference only
  through its quadrature error, 0.3% from that of every candidate point;
- `oscillatingWallStokesLayer`: the wall shear stress from the surface
  traction is the same with both;
- `tractionQuadrature greedy` gives the force files of earlier versions
  unchanged (checked with `MESH_LEVEL=1` for `oscillatingCylinderInChannel`,
  `staticCylinderInChannel`, `translatingCylinderInChannel`,
  `bendingBeamInChannel` and `oscillatingWallStokesLayer`).

## Experimental: thin bodies

**Experimental, work in progress, off by default.** This section describes
`thinBodyBlocking` of `immersedBoundaryForce` (`method cutLink` with
`apertureCoupling yes`). It is not a validated method for closed valves or
for thin curved shells in general, and with the default `off` the code path
and results are unchanged (bit-identical on the four immersed-boundary
tutorials listed below).

### Problem

With `apertureCoupling` the face apertures are computed from the signed
distance at the mesh points (`cutFaceIso`) and the cells inside the bodies
from the signed distance at the cell centres. A closed body of thickness `t`
smaller than the point spacing `h` is therefore seen by the points only for
the fraction `t/h` of the time, as it moves through the mesh, and the fluid
leaks through it when it is not seen.

### Crossing test

For every internal face with owner `P` and neighbour `N`, the segment
`C_P C_N` between the two cell centres is tested against the triangles of the
surface (`triSurfaceSearch::findLine`): the test is exact for any thickness.
If the segment crosses the surface, the face is blocked: its aperture is set
to 0 and the flux through it is that of the body velocity, as for the faces
that are inside a body, so that the pressure equation of `pimpleFluid` treats
it as a solid face. Coupled (processor) faces are tested in the same way with
the segment to the centre of the cell on the other side, so that both sides
of a processor face block the same faces. Only the faces whose cell centres
are inside the bounding box of a body, grown by two cell sizes, are tested.

```text
thinBodyBlocking  off;           // default: unchanged behaviour
thinBodyBlocking  openFaces;     // also block the faces that the points leave
                                 // fully open and a segment crosses
thinBodyBlocking  crossedFaces;  // openFaces, and also the partly cut internal
                                 // faces between two cells whose centres are
                                 // outside the bodies
thinBodyDiffusivity  zero;       // default with thinBodyBlocking; or keep
```

The blocking is by whole faces (binary), not by a partial aperture.

### Pressure equation of a blocked face

In `pimpleFluid` the corrected flux on a face is

```text
phi_f = phiHbyA_f - rAtU_f |S_f| snGrad(p)_f,
phiHbyA_f = alpha_f phiHbyA_f^* + wallFlux_f,
rAtU_f = alpha_f rAtU_f^* + (1 - alpha_f) min(rAtU_P, rAtU_N),
```

where `phiHbyA_f^*` is the flux of the interpolated `HbyA`, `alpha_f` the
aperture and `wallFlux_f` the flux of the body velocity through the solid part
of the face, and `p` solves `sum_f rAtU_f |S_f| snGrad(p)_f = sum_f phiHbyA_f`
(cells inside a body excepted). For a blocked face, `alpha_f = 0` gives
`phiHbyA_f = wallFlux_f`, so the flux is the body flux only if the second term
of `phi_f` vanishes. It does not: `rAtU_f = min(rAtU_P, rAtU_N)` is that of the
fluid when no cell centre is inside the body (the penalty makes `rAtU` small
only in the cells whose centre is inside), so a pressure difference across the
plate drives a flux through a face that was meant to be impermeable. The
leakage of a static plate under a pressure drop was 2-5 % of the open flow
(`thinBodyDiffusivity keep`), and a plate that contains a mesh point but no
cell centre leaked completely.

With `thinBodyDiffusivity zero` (default):

- the faces blocked by a crossing, and the faces that are entirely inside a
  body between two cells that are not entirely solid, get `rAtU_f = 0`: the
  pressure matrix stays a graph Laplacian (symmetric, non-positive
  off-diagonal, zero row sums away from the boundary), and `phi_f` equals
  `wallFlux_f` exactly on them, so the body flux is conserved;
- the faces with a zero diffusivity separate the mesh into regions
  (`regionSplit`, across processors). A region that has no fixed-pressure
  boundary face is singular: the net volume imbalance of its sources is
  spread uniformly over its cells (compatibility, logged as the largest net
  imbalance removed) and one of its cells gets a reference pressure. A region
  with a fixed pressure is untouched. The cells inside a plate form such a
  region, with a zero imbalance;
- the cell-centred pressure gradient, used in the momentum predictor and the
  velocity correction, uses the cell value on each side of a blocked face
  instead of the interpolated face value, as `p` is discontinuous across a
  closed wall; without it a plate thinner than 0.2 cells and no cell centre
  inside it was unstable (static, slip walls, 0.1 cells).

The default `off` path does not call any of this (the registered field is
absent), checked by bit-identity.

### Evidence

All results are with OpenFOAM.com v2412. The test is
`tutorials/fluids/immersedBoundary/thinPlateInDuct`: a closed box plate
spanning a 2 m x 0.2 m duct (cells of h = 10 mm, one cell in the
empty direction, 200 x 20 x 1 cells, laminar, nu = 1e-3, impulsive start,
open ends, fixed pressure at both ends) translated along the duct at
w = 0.1 m/s by `uniformTranslation`. An impermeable plate pushes the fluid at
its own speed, so the volume flux through every cross-section is exactly `w A`
(`A = 0.2 x 0.01 m2`). Q is the outlet flux, and Q/(wA) is the average of the
values of every time step between 2 and 8 s (the start-up transient is
excluded; the instantaneous values vary around it with the position of the
plate in the mesh).

Time-averaged Q/(wA), 2-8 s, plate moving at w (`crossedFaces`; S slip walls,
N no-slip walls; `keep` is the previous behaviour):

| t/h | S off | S keep | S zero | N off | N keep | N zero |
| --- | --- | --- | --- | --- | --- | --- |
| 0.43 | 0.769 | 1.0000 | 1.0000 | 0.442 | 0.946 | 1.0000 |
| 0.5 | 0.889 | 1.0000 | 1.0000 | 0.568 | 0.954 | 1.0000 |
| 1 | 1.0000 | 1.0000 | 1.0000 | 0.964 | 0.981 | 1.0000 |
| 1.7 | 1.0000 | 1.0000 | 1.0000 | 0.980 | 0.987 | 0.9985 |
| 3 | 1.0000 | 1.0000 | 1.0000 | 0.989 | 0.992 | 1.0000 |

The no-slip deficit of 1-5 % with `keep` was the leakage through the blocked
faces, not a wall-layer effect. With `zero` the instantaneous Q/(wA) of the
no-slip plate of 1.7 cells still varies between 0.992 and 1.000.

Static plate under a pressure drop (inlet `p` = 0.01, outlet 0, plate at rest,
S slip or N no-slip walls, `crossedFaces`): leakage fraction at 6 s, the outlet
flux divided by that of the duct without a plate at the same time (the flow
is still accelerating):

| t/h | N off | N keep | N zero | S off | S keep | S zero |
| --- | --- | --- | --- | --- | --- | --- |
| 0.43 | 47 % | 5.1 % | 1e-8 | 33 % | 2.3 % | 2e-8 |
| 1 | 7.4 % | 4.5 % | 3e-8 | 3.5 % | 2.0 % | 2e-8 |
| 3 | 2.4 % | 2.0 % | 1e-9 | 1.1 % | 0.9 % | 1e-9 |

With `zero` the outlet flux is at the level of the solver tolerance (at most
2e-10 m3/s over the run, against 2.6e-5 m3/s without a plate, no-slip, and
6e-5 slip). Dependence, 0.43 cells, `zero`, leakage fraction at 6 s: time step
halved or doubled 2e-8 and 2e-8 (slip), 8e-8 and 1e-10 (no-slip); mesh
refined to h = 5 mm (same absolute thickness ratio of 0.43 cells) 2e-9 and
3e-9; thickness 0.1 cells 3e-10 and 8e-11. With `keep` a plate of 0.43 cells
on the h = 5 mm mesh, which contains a mesh point but no cell centre, leaks
99.6 %. For a static plate straddling a processor face the outlet flux is at most
2e-12 m3/s serial and 2e-16 m3/s in 4 ranks.

Moving mesh and deforming body (research-branch ALE test cases A-G, plates of
three cells, copied for this study; maximum deviation of the lab-frame
velocity from the exact uniform value, divided by the plate speed; same value
with the switch off, `keep`, `zero`): fixed mesh and moving plate (B)
1.5e-5, 3.1e-5, 1.5e-4; plate crossing the duct for 10 s (F) the same;
mesh moving with the plate (D) 0, 0, 0; mesh moving, plate at rest (C) 0, 0,
0; mesh and plate moving (G) 1.8e-5, 2.4e-5, 2.4e-5 and at half the cell size
5.8e-6, 1.1e-5, 1.8e-4; deforming plate (E) 2.0e-4, 2.5e-4, 6.8e-5 and at
half the cell size 9.2e-5, 9.8e-5, 2.9e-5. The volume flux is 1.000 of the
exact one in the fixed-mesh cases and the geometric conservation error is
5e-13 to 3e-8 of the reference flux in all of them. The deviations with `zero`
remain small but are larger than with `off` in three cases (B, F, G at h/2),
a few 1e-4 of the plate speed.

Results below are from the research branch the code comes from (the
crossing code of the first commit, not re-run for this pull request, and with
the earlier `keep` diffusivity):

- Thin shells of 0.3 cells thickness (a plate tilted by 30 degrees and a
  parabolic arc of 4.5 cm sag, closed prisms, translated along x at w, slip
  walls), time-averaged Q/(wA) over 2-8 s: tilted 0.72 (off) to 0.89
  (`crossedFaces`), arc 0.75 to 0.94. The shortfall is a start-up transient:
  Q/(wA) of the tilted shell reaches 0.995 at 10 s (0.999 for the arc). With
  no-slip walls the values are 0.40 to 0.57 and 0.42 to 0.66. These shells
  are not within 2 % of their swept volume.

Regression with the switch off (default), same compiler and OpenFOAM as the
base `ba563ba6a`, serial, 1 s: `translatingCylinderInChannel`,
`staticCylinderInChannel`, `oscillatingCylinderInChannel` and
`bendingBeamInChannel` (`MESH_LEVEL=1`) give fields (`U`, `p`, `phi`, ...) and
`postProcessing` files that are bit-identical to those of the unmodified base,
also with the pressure-equation changes.
A 4-rank run (simple `(4 1 1)`; the plate crosses the processor faces at
`x` = 0.5 and 1 m during 2-8 s) gives the same time-averaged Q/(wA) as the
serial run to the five digits shown (slip 1.00000, no-slip 1.00000); the
instantaneous outlet flux differs by at most 8e-4 of `w A` (1.7e-3 and 7.8e-3
of `w A` for the no-slip plates of 0.43 and 1.7 cells), which is the level of
the solver tolerance in parallel.

### Known limitations and open questions

- Faces are blocked entirely or not at all, not by a partial aperture; the
  flux changes in steps of order `h^2 w` as the set of blocked faces changes.
  A hole or gap narrower than a cell, or a plate edge that ends inside a
  face, is not represented properly.
- Oblique and curved thin shells transmit less than their swept volume during
  the start-up (0.3 cells thickness, above, with the earlier diffusivity).
- Faces entirely inside a body between two cells that are not entirely solid
  are blocked too (`zero`), which is the intended meaning for a thin body, and
  not tested for thick bodies (use `off` for those).
- Free edges, gaps, holes, leaflet tips and tapered leaflets are not treated:
  a face is blocked only if the segment between two cell centres crosses a
  triangle or the face is inside the body, so a slit narrower than a cell
  between two blocked faces is closed by the staircase.
- The force estimators (momentum exchange, surface traction) with the switch
  on, deforming bodies and contact were not examined.
- The cost of the crossing search (one `findLine` per candidate face at each
  update, and a `regionSplit` per pressure corrector) was not measured on
  large meshes; meshes with skewed processor faces were not tested.
- Tested only for prescribed rigid translation at constant velocity (and the
  deforming plate of the ALE tests), laminar flow and a mesh with `h` = 10 mm
  (5 mm for some checks).
- The "faces blocked" count in the log adds the coupled faces on both sides.

## Provenance

This library is a rewrite, for prescribed motion and as a finite volume
option, of the immersed boundary library and `pimpleHFDIBFoam` solver
contributed to cardiacFoam in
[solids4foam/cardiacFoam#20](https://github.com/solids4foam/cardiacFoam/pull/20)
by Sairam Pamulaparthi Venkata (UCD). That code is based on:

- **openHFDIB-DEM**, <https://github.com/techMathGroup/openHFDIB-DEM>,
  developed mostly by members of the techMathGroup of the Institute of
  Thermomechanics of the Czech Academy of Sciences and the Monolith group of
  the University of Chemistry and Technology, Prague; its main contributors
  are Martin Isoz, Martin Kotouč Šourek, Ondřej Studeník and Petr Kočí. The
  imported version is the library in `src/HFDIBDEM` and the
  `pimpleHFDIBFoam` solver at commit
  [`b4ea321`](https://github.com/techMathGroup/openHFDIB-DEM/commit/b4ea321)
  (2024-08-19), on the `main` and `Ver2.6` branches.
- **openHFDIB**, <https://github.com/fmuni/openHFDIB>, by Federico Municchi,
  from which openHFDIB-DEM is derived.

The lineage of the code is:

1. openHFDIB-DEM `b4ea321`, imported into
   [xenosim-erc/immersedBoundaryRigidMotion](https://github.com/xenosim-erc/immersedBoundaryRigidMotion)
   by Sairam Pamulaparthi Venkata;
2. the discrete element (contact, virtual mesh) parts removed and the code
   ported to OpenFOAM v2512 by Philip Cardiff (February 2026, branch
   `sairam-leftHeart-IB-Verf`);
3. prescribed motion laws (sinusoidal translation, bending, heart valves) and
   benchmark cases added by Sairam Pamulaparthi Venkata (June 2026,
   `xenosim-erc/openHFDIB-Benchmarks-v2412`, and
   solids4foam/cardiacFoam#20 commit `251e887`);
4. rewritten here for solids4foam by Philip Cardiff: the discrete element
   method, free bodies, body addition, adaptive mesh refinement and restart
   files are removed; bodies are moved identically on every processor at the
   current time; and the solver is replaced by the finite volume option.

The occupancy (`immersedBody::addOccupancy`) follows the openHFDIB-DEM
`nonConvexBody`, and the forcing and force calculation
(`immersedBoundaryForce`, `immersedBody::force`) follow the openHFDIB-DEM
`pimpleHFDIBFoam` solver and `immersedBody`. The sinusoidal translation,
the `quadraticBend` and `customProfileBend` deformations, and the
`valveSliceAxis` and `valveMorph` valve motions and their time laws follow
the versions of Sairam Pamulaparthi Venkata in solids4foam/cardiacFoam#20,
except the closed-to-open valve morph, which is from his later working copy
of xenosim-erc/immersedBoundaryRigidMotion
(`src/HFDIBDEM/geomModels/stlBased/stlBased.C`, June 2026, not committed). The openHFDIB-DEM
interpolation of the velocity at the immersed boundary (`lineInt`,
`leastSquares`), which the benchmark cases do not use, is not included.

With the `incremental` method, the `vertexFraction` occupancy and the same
settings, `immersedBoundaryForce` reproduces the forces of `pimpleHFDIBFoam`
for the static cylinder in a channel to six significant figures. The
`penalty` method, the `signedDistance` occupancy and the `volumeFraction`
weighting are new in solids4foam.

## Licence

openHFDIB is distributed under the GNU Lesser General Public License version
3, and the openHFDIB-DEM repository under the GNU General Public License
version 3 (its source files state the GNU Lesser General Public License
version 3). The Lesser General Public License allows the code to be
distributed under the GNU General Public License version 3, under which this
library is distributed as part of solids4foam.

## References

- Municchi, F. openHFDIB. <https://doi.org/10.5281/zenodo.3940236>
- Isoz, M., Kotouč Šourek, M., Studeník, O., Kočí, P. (2022). Hybrid
  fictitious domain-immersed boundary solver coupled with discrete element
  method for simulations of flows laden with arbitrarily-shaped particles.
  Computers & Fluids, 244, 105538.
  <https://doi.org/10.1016/j.compfluid.2022.105538>
- Studeník, O., Isoz, M., Kotouč Šourek, M., Kočí, P. (2024). OpenHFDIB-DEM:
  An extension to OpenFOAM for CFD-DEM simulations with arbitrary particle
  shapes. SoftwareX, 27. <https://doi.org/10.1016/j.softx.2024.101871>
- Uhlmann, M. (2005). An immersed boundary method with direct forcing for the
  simulation of particulate flows. Journal of Computational Physics, 209,
  448-476.
