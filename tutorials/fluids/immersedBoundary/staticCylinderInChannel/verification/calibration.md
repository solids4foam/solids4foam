# Calibration of the time step and the ranks

solids4foam `f33b244a2` with OpenFOAM v2412, on one node of 2 x 96-core AMD
EPYC 9684X under SLURM, October 2026. The data are in
`reference/mpi_scaling.csv`.

## Fixed time step: 0.001 s

Two `cutLink` pilots on `MESH_LEVEL=4` (288 640 cells, 80 cells across the
cylinder) from rest to `t = 3 s`, with 4 ranks:

```bash
./Allverify generate --study calibration --methods cutLink --levels 4 \
    --tag dt0p001 --deltaT 0.001 --end-time 3 --write-interval 3 --ranks 4
# and likewise --tag dt0p00125 --deltaT 0.00125
```

| deltaT (s) | max Co, developed flow | Cd at t = 3 s |
| ---------- | ---------------------- | ------------- |
| 0.001      | 0.420                  | 5.58194       |
| 0.00125    | 0.525                  | 5.58192       |

The maximum Courant number of the developed flow on Level 4 is
`420 deltaT` (1/s), so that `Co <= 0.5` requires `deltaT <= 0.00119 s`;
`0.001 s` is the largest round time step below it that divides the end and
write times. It was chosen from this criterion only, independently of the
drag, and is used for all the formulations and levels (`max Co` 0.05, 0.105,
0.21 and 0.42 on levels 1-4). The legacy incremental formulation reaches
`Co = 0.52-0.54` in its first time steps on Level 4 only (0.42 after).

## Ranks: 1, 1, 2 and 8 on levels 1-4

`cutLink`, Level 4: the same 300 time steps from the developed state of the
`0.001 s` pilot with 1, 2, 4 and 8 ranks in one 8-CPU allocation:

```bash
sbatch --nodes=1 --ntasks=8 --cpus-per-task=1 --mem=16G --time=01:30:00 \
    --export=ALL --output=work/calibration/scaling-%j.out \
    scripts/mpi_scaling.sh work/calibration/cutLink-L4-dt0p001 300 "1 2 4 8"
```

| Ranks | Cells per rank | Wall time (s) | Speedup | Efficiency |
| ----- | -------------- | ------------- | ------- | ---------- |
| 1     | 288 640        | 364.8         | 1.00    | 1.00       |
| 2     | 144 320        | 172.0         | 2.12    | 1.06       |
| 4     | 72 160         | 88.4          | 4.13    | 1.03       |
| 8     | 36 080         | 39.8          | 9.17    | 1.15       |

The scaling is linear to 8 ranks, and the drag is the same to `3e-7` with all
the rank counts. The cases use 1 rank on levels 1 and 2, 2 on Level 3 and 8 on
Level 4, the same for all the formulations, with the `simple` decomposition
in `x`.

## 8 or 16 ranks for the long legacy incremental runs

Legacy incremental (`couplingCoeff 0.2`, `nOuterCorrectors 5`), Level 4: the
same 200 time steps from its developed state at `t = 15 s` with 8 and 16
ranks (stage C):

```bash
sbatch --nodes=1 --ntasks=16 --cpus-per-task=1 --mem=32G --time=01:30:00 \
    --export=ALL --output=work/final-C/scaling-%j.out \
    scripts/mpi_scaling.sh work/final-C/base-incremental-L4-cc0p2-n5 200 "8 16"
```

| Ranks | Wall time (s) | s per step | CPU hours | Cd mean | Cd range | Cl mean | Cl range |
| ----- | ------------- | ---------- | --------- | ------- | -------- | ------- | -------- |
| 8     | 253.2         | 1.266      | 0.563     | 5.733872 | 3.2e-4  | 0.012273 | 6.6e-5  |
| 16    | 146.5         | 0.732      | 0.651     | 5.733769 | 1.3e-5  | 0.012295 | 2.5e-6  |

The speedup is 1.73 for 16% more CPU hours, but the limit cycle of the forces
of this case (8 ranks) decays with 16 ranks: its amplitude depends on the
decomposition, while the mean drag agrees to `2e-5`. The long Level 4 legacy
incremental runs therefore keep 8 ranks, so that all of them share the
decomposition of the other Level 4 cases.

## Process placement

The SLURM installation does not confine the jobs to their CPUs
(`TaskPlugin=task/none`): with `mpirun --bind-to core` the ranks of two
concurrent jobs were bound to the same cores and ran at half speed. The
cases use `mpirun --bind-to none`, with which every rank ran at about 100% of
a separate core.
