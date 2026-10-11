#!/usr/bin/env python3
"""Replay restart case with per-face traction output on the flag (tip localisation).

Builds the same restart case as hron_turek_energy.py (identical motion, function
object and restart, so Q_in is reproduced) and appends a coded function object
`plateTraction` that writes, per MPI rank, a binary file traction_<rank>.bin:

  header   int64 nf, int64 ncoef; then per local face
           xmin xmax ymin ymax len(=|Sf|/depth) and the 2*ncoef replay-table
           coefficients of the face (mean of its vertices), all float64, with the
           coordinates taken from the undeformed points
  records  float64 t, then per face fpx fpy fvx fvy (N/m, same definition as the
           energy function object), for every step with t > --write-from

    python3 hron_turek_tip.py --source <replay> --restart 25 --end 33.2 \
        --name t_1x --work ~/ht_tip/work
"""
from __future__ import annotations
import argparse, subprocess, sys
from pathlib import Path

HERE = Path(__file__).resolve().parent

FO = r"""
    plateTraction
    {
        type            coded;
        libs            ( utilityFunctionObjects );
        name            plateTraction;
        codeInclude
        #{
            #include <fstream>
            #include <map>
            #include <vector>
            #include <cmath>
            #include "fvcGrad.H"
            #include "pointFields.H"
        #};
        codeExecute
        #{
            static bool started = false;
            const fvMesh& m = mesh();
            const label pid = m.boundaryMesh().findPatchID("plate");
            const polyPatch& pp = m.boundaryMesh()[pid];
            const scalar t = m.time().value();
            const std::string fname =
                (m.time().globalPath()/("traction_" + std::to_string(Pstream::myProcNo()) + ".bin")).c_str();
            if (!started)
            {
                std::ifstream in((m.time().globalPath()/"constant/plateMotion.tab").c_str());
                long n = 0, ncoef = 0;
                double omega = 0;
                in >> omega >> n >> ncoef;
                std::map<std::pair<long, long>, std::vector<double> > table;
                for (long i = 0; i < n; i++)
                {
                    long kx, ky;
                    in >> kx >> ky;
                    std::vector<double> c(2*ncoef);
                    for (long j = 0; j < 2*ncoef; j++) { in >> c[j]; }
                    table[std::make_pair(kx, ky)] = c;
                }
                const pointIOField p0
                (
                    IOobject
                    (
                        "points", m.time().constant(), m.dbDir()/polyMesh::meshSubDir,
                        m.time(), IOobject::MUST_READ, IOobject::NO_WRITE, false
                    )
                );
                std::ofstream os(fname.c_str(), std::ios::binary | std::ios::trunc);
                const int64_t nf = pp.size(), nc = ncoef;
                os.write(reinterpret_cast<const char*>(&nf), sizeof(nf));
                os.write(reinterpret_cast<const char*>(&nc), sizeof(nc));
                const vectorField& Sf0 = m.Sf().boundaryField()[pid];
                forAll(pp, fI)
                {
                    const face& f = pp[fI];
                    std::vector<double> c(2*ncoef, 0.0);
                    double xl = 1e30, xh = -1e30, yl = 1e30, yh = -1e30;
                    forAll(f, k)
                    {
                        const point& q = p0[f[k]];
                        xl = std::min(xl, double(q.x())); xh = std::max(xh, double(q.x()));
                        yl = std::min(yl, double(q.y())); yh = std::max(yh, double(q.y()));
                        auto it = table.find
                        (
                            std::make_pair(long(std::llround(q.x()*%(key)r)), long(std::llround(q.y()*%(key)r)))
                        );
                        if (it == table.end())
                        {
                            FatalErrorInFunction << "No replay row for " << q << exit(FatalError);
                        }
                        for (long j = 0; j < 2*ncoef; j++) { c[j] += it->second[j]/f.size(); }
                    }
                    const double len = mag(Sf0[fI])/%(depth)r;
                    const double hd[5] = {xl, xh, yl, yh, len};
                    os.write(reinterpret_cast<const char*>(hd), sizeof(hd));
                    os.write(reinterpret_cast<const char*>(c.data()), 2*ncoef*sizeof(double));
                }
                started = true;
            }
            if (t <= %(tw)r) { return true; }
            const volScalarField& p = m.lookupObject<volScalarField>("p");
            const volVectorField& U = m.lookupObject<volVectorField>("U");
            const tmp<volTensorField> tgU(fvc::grad(U));
            const tensorField& gUb = tgU().boundaryField()[pid];
            const symmTensorField devR(-%(nu)r*dev(twoSymm(gUb)));
            const scalarField& pb = p.boundaryField()[pid];
            const vectorField& Sf = m.Sf().boundaryField()[pid];
            std::vector<double> rec(1 + 4*pp.size());
            rec[0] = t;
            forAll(pb, i)
            {
                const vector fp = %(rho)r*pb[i]*Sf[i]/%(depth)r;
                const vector fv = %(rho)r*(Sf[i] & devR[i])/%(depth)r;
                rec[1 + 4*i] = fp.x(); rec[2 + 4*i] = fp.y();
                rec[3 + 4*i] = fv.x(); rec[4 + 4*i] = fv.y();
            }
            std::ofstream os(fname.c_str(), std::ios::binary | std::ios::app);
            os.write(reinterpret_cast<const char*>(rec.data()), rec.size()*sizeof(double));
        #};
    }
"""


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--source", required=True, type=Path)
    ap.add_argument("--restart", required=True)
    ap.add_argument("--end", required=True, type=float)
    ap.add_argument("--name", required=True)
    ap.add_argument("--work", type=Path, required=True)
    ap.add_argument("--write-from", type=float, required=True, help="first time written is > this")
    a = ap.parse_args()
    subprocess.run([sys.executable, str(HERE / "hron_turek_energy.py"), "--source", str(a.source),
                    "--restart", a.restart, "--end", str(a.end), "--name", a.name,
                    "--work", str(a.work)], check=True)
    case = (a.work / a.name).resolve()
    fn = case / "system/functions"
    s = fn.read_text()
    i = s.rindex("}")
    fn.write_text(s[:i] + FO % dict(key=1e5, nu=1.0e-3, rho=1000.0, depth=0.015, tw=a.write_from) + s[i:])
    with open(case / "energy_setup.txt", "a") as fh:
        fh.write(f"traction_write_from {a.write_from}\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
