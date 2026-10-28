# Closed-loop circulation with a pipe surrogate: `pipeClosedLoopMonolithic`

---

Prepared by David Kelly

---

## Tutorial Aims

- Solve the closed-loop pipe circuit of
  [`pipeClosedLoop`](../pipeClosedLoop/README.md) as a single system, with no
  CFD, by replacing the 3D pipe by a lumped resistance and inertance.
- Provide a reference for the CFD-coupled case, which has no partitioned
  coupling error.

---

## Case Overview

### Problem description

The circuit is the one of `pipeClosedLoop`: a pipe whose outlet drains into the
arterial side, the heart (atrium, mitral valve, ventricle, aortic valve), and
the return of the blood to the pipe inlet.

```text
C1 -R1-L1- C2 -R2- C3 -R3- atrium -mitral- ventricle -aortic(+R4)- C4 -pipe- C1
```

In `pipeClosedLoop` the pipe is a 3D CFD domain coupled to the circuit through
its two end patches. Here the pipe is replaced by its rigid, laminar
Poiseuille surrogate, a resistor `Rpipe` in series with an inductor `Lpipe`,
so the whole loop is one `circulationNetwork`:

- `Rpipe = 8 mu L / (pi r^4)`, with a viscosity `mu` of 0.004 Pa s, a pipe
  length `L` of 0.3 m and a radius `r` of 0.02 m;
- `Lpipe = rho L / A`, with a density `rho` of 1060 kg/m3 and the pipe
  cross-sectional area `A`.

The circuit is in `constant/circulationProperties`, in SI units, with a
cardiac period of 0.8 s. It is the same circuit as in `pipeClosedLoop`, with
the two coupling ports replaced by `Rpipe` and `Lpipe`.

### Solver settings

- `endTime` is 2.4 s, three cardiac cycles.
- `deltaT` is 0.1 ms. This is only the **output interval**: the ODE solver
  takes whatever internal steps its tolerances need.
- The ODE solver is `RKDP45` with tight tolerances. The hard switching of
  the valves makes it take many small steps.

---

## Running the Case

The case is run with `circulationFoam`, which is built along with solids4foam.
No mesh is needed:

```bash
./Allrun
```

This runs `circulationFoam -check`, which describes the circuit and checks the
set-up, and then solves it. The results are written to
`postProcessing/circulation/0/circulation.dat`: a row every `deltaT`, with the
pressure, volume and flow of every element. The case is cleaned with
`./Allclean`.

---

## Expected Results

The coupled CFD solution of `pipeClosedLoop` differs from this solution by the
modelling error of the surrogate: the resolved pipe is not exactly its rigid
laminar surrogate. At the coarse time steps used in the tutorials, the
coupling error dominates this difference.
