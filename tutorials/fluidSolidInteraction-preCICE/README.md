The preCICE coupling library (http://precice.org/) can provide FSI coupling
between the solids4foam solid solvers and the OpenFOAM fluid solvers (or any
other fluid solvers supported by preCICE).

Both preCICE and the preCICE OpenFOAM adapter need to be installed:
- https://precice.org/quickstart.html
- https://precice.org/adapter-openfoam-overview.html

You can then try out the preCICE tutorials (https://precice.org/tutorials.html)
that demonstrate solids4foam. For example, see:
- https://precice.org/tutorials-perpendicular-flap.html [pending pull request https://github.com/precice/tutorials/pull/286]
These tutorials can be downloaded using the included `installTutorials` script, i.e. `./installTutorials`.

Note: the OpenFOAM adapter needs to have the solids4foam update [https://github.com/precice/openfoam-adapter/pull/236].