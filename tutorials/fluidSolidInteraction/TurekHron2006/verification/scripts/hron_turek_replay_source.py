#!/usr/bin/env python3
"""Record the converged coupled flag trajectory of an FSI3 run (replay source).

Restarts a completed FSI3 IQN-ILS run of the mesh study (default: fluid 2x,
solid 2x, which ends at t = 7 s in its periodic state) for a further
`--duration` seconds with a `coded` function object that writes, at every time
step, the positions of the points of the fluid `plate` patch to
`processor<N>/postProcessing/plateMotion.dat` (local mesh point ids on the first
line). The solid is a total-Lagrangian St. Venant-Kirchhoff material, so the
restart does not need the incremental fields (`restart no`); the IQN-ILS
history is rebuilt within a few steps. The trajectory is later fitted and
replayed on fluid meshes by hron_turek_replay.py.

    python3 scripts/hron_turek_replay_source.py --run iqnils_mesh_2x --cores 8
"""

from __future__ import annotations

import argparse
import re
import shutil
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import hron_turek_verification as driver  # noqa: E402

CODE = r"""
        const fvMesh& m = mesh();
        const label pid = m.boundaryMesh().findPatchID("plate");
        const labelList& mp = m.boundaryMesh()[pid].meshPoints();
        const pointField& pts = m.points();
        static bool first = true;
        mkDir(m.time().path()/"postProcessing");
        std::ofstream os
        (
            (m.time().path()/"postProcessing/plateMotion.dat").c_str(),
            std::ios::app
        );
        os.precision(12);
        if (first)
        {
            os << "IDS";
            forAll(mp, i) { os << " " << mp[i]; }
            os << "\n";
            first = false;
        }
        os << m.time().value();
        forAll(mp, i) { os << " " << pts[mp[i]].x() << " " << pts[mp[i]].y(); }
        os << "\n";
"""


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--run", default="iqnils_mesh_2x")
    parser.add_argument("--cores", type=int, required=True)
    parser.add_argument("--duration", type=float, default=1.2)
    parser.add_argument("--name", default=None)
    args = parser.parse_args()
    src = driver.WORK_ROOT / args.run
    name = args.name or f"replaysrc_{args.run}"
    case = driver.WORK_ROOT / name
    if case.exists():
        shutil.rmtree(case)
    case.mkdir()
    for sub in ("constant", "system", "0"):
        shutil.copytree(src / sub, case / sub, symlinks=True)
    for processor in sorted(src.glob("processor*")):
        destination = case / processor.name
        destination.mkdir()
        shutil.copytree(processor / "constant", destination / "constant")
        shutil.copytree(processor / "7", destination / "7")
    # Number of subdomains of the dictionaries must match the ranks
    for dictionary in ("system/decomposeParDict", "system/fluid/decomposeParDict",
                       "system/solid/decomposeParDict"):
        driver.replace_entry(case / dictionary, "numberOfSubdomains", str(args.cores))
    solid = case / "constant/solid/solidProperties"
    text, count = re.subn(r'(Coeffs"\s*\{)', r"\g<1>\n    restart no;", solid.read_text(), count=1)
    if count != 1:
        driver.fail("Could not set restart no")
    solid.write_text(text)
    control = case / "system/controlDict"
    driver.replace_entry(control, "endTime", f"{7 + args.duration:.8g}")
    driver.replace_entry(control, "writeInterval", "100000000")
    functions = case / "system/functions"
    text = functions.read_text()
    block = ("\n   plateMotion\n   {\n       type            coded;\n"
             "       libs            ( utilityFunctionObjects );\n"
             "       name            plateMotion;\n       region          fluid;\n"
             "       codeInclude\n       #{\n           #include <fstream>\n       #};\n"
             "       codeExecute\n       #{" + CODE + "       #};\n   }\n")
    functions.write_text(text.rstrip().rstrip("}") + block + "}\n")
    log = case / "log.solids4Foam"
    with log.open("w") as handle:
        result = subprocess.run(
            ["mpirun", "-np", str(args.cores), "solids4Foam", "-parallel"],
            cwd=case, stdout=handle, stderr=subprocess.STDOUT, text=True)
    text = log.read_text(errors="replace")
    if result.returncode or not re.search(r"^End\s*$", text, re.MULTILINE):
        driver.fail(f"{name} did not complete; see {log}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
