---
sort: 1
---

# `solids4Foam` solver

`solids4Foam` is the single solver application provided by solids4foam. It is
used for solid mechanics, fluid dynamics and fluid-solid interaction analyses
alike: the physics, the discretisation and the solution algorithm are all
selected at run time, so the same executable runs every case in the toolbox.

---

## Usage

```bash
solids4Foam                 # serial
mpirun -np 4 solids4Foam -parallel   # parallel
```

The type of analysis is set by `constant/physicsProperties`:

```c++
type    solid;    // solid | fluid | fluidSolidInteraction
```

Each choice then reads its own dictionary:

| `type` | Dictionary | Base class |
| ------ | ---------- | ---------- |
| `solid` | `constant/solidProperties` | `solidModel` |
| `fluid` | `constant/fluidProperties` | `fluidModel` |
| `fluidSolidInteraction` | `constant/fsiProperties` | `fluidSolidInterface` |

---

## The solver code

The solver is deliberately thin; it contains no details of the physics or the
discretisation. Its source is
[`applications/solvers/solids4Foam/solids4Foam.C`](https://github.com/solids4foam/solids4foam/blob/development/applications/solvers/solids4Foam/solids4Foam.C):

```c++
int main(int argc, char *argv[])
{
#   include "setRootCase.H"
#   include "createTime.H"
#   include "solids4FoamWriteHeader.H"

    // Create the general physics class
    autoPtr<physicsModel> physics = physicsModel::New(runTime);

    while (runTime.run())
    {
        // Update deltaT, if desired, before moving to the next step
        physics().setDeltaT(runTime);

        runTime++;

        if (physics().printInfo())
        {
            Info<< "Time = " << runTime.timeName() << nl << endl;
        }

        // Solve the mathematical model
        physics().evolve();

        // Let the physics model know the end of the time-step has been reached
        physics().updateTotalFields();

        if (runTime.outputTime())
        {
            physics().writeFields(runTime);
        }

        if (physics().printInfo())
        {
            Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
                << "  ClockTime = " << runTime.elapsedClockTime() << " s"
                << nl << endl;
        }
    }

    physics().end();

    Info<< nl << "End" << nl << endl;

    return(0);
}
```

A run-time selectable `physicsModel` object encapsulates the specifics, and
virtual functions such as `evolve()` tell that object to solve its governing
equations for the current time step. Adding a new solid model, fluid model or
coupling scheme therefore requires no change to the solver.

The `physicsModel` base class has three derived families:

- `solidModel`, see
  [solid models](https://www.solids4foam.com/documentation/solid-models/);
- `fluidModel`, see
  [fluid models](https://www.solids4foam.com/documentation/fluid-models/);
- `fluidSolidInterface`, which itself creates one `fluidModel` and one
  `solidModel`, see
  [fluid-solid interfaces](https://www.solids4foam.com/documentation/fluid-solid-interfaces/).

---

## Source

[`applications/solvers/solids4Foam`](https://github.com/solids4foam/solids4foam/tree/development/applications/solvers/solids4Foam)
