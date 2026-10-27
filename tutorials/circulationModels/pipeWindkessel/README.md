# pipeWindkessel: a CFD pipe coupled to a Windkessel

A 5 m, 0.5 m radius pipe (25,004 cells) solved with stock `pimpleFoam`. The
inlet prescribes the flow `Q(t) = sin(2 pi t)` m3/s. The outlet pressure comes
from a three-element Windkessel described in `constant/circulationProperties`
and solved by `circulation0D`. The outlet is a single boundary condition
plus a dictionary, whichever coupling scheme is used.

## Choosing the coupling

```c++
// constant/circulationProperties
coupling
{
    scheme      implicit;   // explicit | implicit
    order       2;          // 1 .. 4
    rho         1060;       // the CFD solves p/rho
}
```

```c++
// 0/p
outlet { type circulationPressure; value uniform 0; }
```

That is the whole interface: the eight combinations of `scheme` and `order`
are the eight coupling algorithms.

| | explicit |
| --- | --- |
| exchanges per step | one |
| order 1 | `Q(t_n)` held over the step |
| order 2 | line through `Q(t_n-1)`, `Q(t_n)`, extrapolated |
| order 3 / 4 | quadratic / cubic through earlier flows, extrapolated |

| | implicit |
| --- | --- |
| exchanges per step | one per outer corrector (new flux) |
| order 1 | the current iterate held over the step |
| order 2 | line from the start flow to the current iterate |
| order 3 / 4 | adds one / two earlier flows to order 2 |

Orders 3 and 4 want a higher-order CFD time scheme so that the CFD does not
dominate the error. The stock `backward` scheme is used here, which limits the
coupled solution to second order (the coupler prints a warning). `makeCase`
seeds the flow history with the exact earlier flows so that every order starts
at full order.

Coupling ports are bound to patches by name: the `outlet` port in the circuit
has `patch outlet;`. A second coupled patch would be a second port; nothing else
changes.

## Running

```bash
# one case: implicit, order 2, deltaT 0.05 (about a minute)
./Allrun
# -> runs/implicit/2/0.01
./makeCase implicit 2 0.01
cd runs/implicit/2/0.01 && blockMesh && pimpleFoam
# many cases
./runLadder "explicit implicit" "1 2" "0.1 0.05 0.02 0.01 0.005" 10
# error against the exact solution, observed orders
python3 analyse.py
```

## Exact solution

With `R = 10, Rch = 1, C = 0.01, Pc(0) = 0`, `dPc/dt = (Q - Pc/R)/C` gives

```text
Pc(t) = A R/(1+(wRC)^2) [ sin wt - wRC cos wt + wRC exp(-t/RC) ]
Po(t) = Pc(t) + Rch A sin wt
```

and the error is the RMSE of `Po` over the steps in `0 < t <= 1`.

## Results

All 40 cases (explicit and implicit, orders 1-4, five time steps each) were run
and compared with the exact solution above (`python3 analyse.py`). The pairwise
orders of accuracy of the RMSE of `Po` converge to 1, 2, 3 and 4. The orders 3
and 4 in this table were obtained with third- and fourth-order BDF time schemes
in the CFD, which are not yet part of solids4foam; with `backward` they are
limited by the CFD time error:

| dt | explicit 1 | explicit 2 | explicit 3 | explicit 4 |
| --- | ---: | ---: | ---: | ---: |
| 0.02 -> 0.01 | 1.029 | 2.024 | 3.037 | 4.020 |
| 0.01 -> 0.005 | 1.015 | 2.012 | 3.019 | 4.012 |

| dt | implicit 1 | implicit 2 | implicit 3 | implicit 4 |
| --- | ---: | ---: | ---: | ---: |
| 0.02 -> 0.01 | 0.982 | 2.001 | 3.012 | 4.004 |
| 0.01 -> 0.005 | 0.991 | 2.001 | 3.006 | 4.003 |

The circuit ODE solver is `Trapezoid` (1e-10) for orders 1 and 2 and `RKDP45`
(1e-12) for orders 3 and 4 (`makeCase` sets it; the tolerances are on the
compliance volume `V = C Pc`, so they are scaled by `C = 0.01`).

The accepted state of every step is logged in
`postProcessing/circulationCoupling/0/coupling.dat`.
