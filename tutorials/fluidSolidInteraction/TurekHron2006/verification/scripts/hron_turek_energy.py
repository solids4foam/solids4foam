#!/usr/bin/env python3
"""Energy diagnostics of the FSI3 trajectory replays (hron_turek_replay.py).

Restarts a finished replay case (prescribed, periodic flag motion fitted to the
coupled 2x FSI3 run) from one of its written times, in a copy under --work, with
a coded function object that logs at every step (line prefix ENRG):

  t, pressure and viscous power on the flag (W/m), pressure and viscous force on
  the flag (N/m, to check against forcesPlate), the discrete work increments
  sum_f f.(d^n - d^(n-1)) of each harmonic h = 1..4 and direction x, y of the
  replayed motion (J/m), and the generalised forces Q_a = sum_f f_c a_hc(f),
  Q_b = sum_f f_c b_hc(f) on the cos and sin shapes of h = 1, 2 (J/m; see
  energy_analysis.py), total and pressure-only for h = 1, y.

Force f = rho (p Sf + Sf . devReff), devReff = -nu dev(twoSymm(grad U)), as in
the forces function object (laminar). Optionally the motion amplitude is scaled
by --amplitude-factor and the frequency by --frequency-factor, both blended in
smoothly over --blend s from the restart time (phase kept continuous), to
measure dW/dA and dW/domega at fixed mesh.

    python3 hron_turek_energy.py --source <replay case> --restart 32 --end 33.2 \
        --name e_2x_dt0.0005 [--amplitude-factor 1.05]
"""
from __future__ import annotations
import argparse, re, shutil, subprocess, sys
from pathlib import Path

RHO, NU, DEPTH, KEY = 1000.0, 1.0e-3, 0.015, 1e5

MOD = r"""
        const scalar kAmp = %(k)r, rFreq = %(r)r, tB = %(tb)r, tR = %(tr)r;
        auto sm = [](scalar x) -> scalar
        {
            x = min(max(x, 0.0), 1.0);
            return x*x*(3.0 - 2.0*x);
        };
        auto Kf = [&](scalar tt) -> scalar { return 1.0 + (kAmp - 1.0)*sm((tt - tB)/tR); };
        auto Tf = [&](scalar tt) -> scalar
        {
            const scalar x = (tt - tB)/tR;
            const scalar I = x <= 0 ? 0.0 : (x <= 1 ? x*x*x - 0.5*x*x*x*x : 0.5 + (x - 1.0));
            return (tt - t0) + (rFreq - 1.0)*tR*I;
        };
"""

FO = r"""
    plateEnergy
    {
        type            coded;
        libs            ( utilityFunctionObjects );
        name            plateEnergy;
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
            static bool loaded = false;
            static std::vector<std::vector<double> > fc;
            static scalar omega = 0;
            static label ncoef = 0;
            const fvMesh& m = mesh();
            const label pid = m.boundaryMesh().findPatchID("plate");
            const polyPatch& pp = m.boundaryMesh()[pid];
            const scalar t = m.time().value();
            const scalar dt = m.time().deltaTValue();
            const scalar t0 = %(t0)r;
            %(mod)s
            auto s = [&](scalar tt) -> scalar
            {
                scalar r = min(max((tt - t0)/%(ramp)r, 0.0), 1.0);
                return r*r*(3.0 - 2.0*r);
            };
            if (!loaded)
            {
                std::ifstream in((m.time().globalPath()/"constant/plateMotion.tab").c_str());
                long n = 0;
                in >> omega >> n >> ncoef;
                std::map<std::pair<long, long>, std::vector<double> > table;
                for (long i = 0; i < n; i++)
                {
                    long kx, ky;
                    in >> kx >> ky;
                    std::vector<double> c(2*ncoef);
                    for (label j = 0; j < 2*ncoef; j++) { in >> c[j]; }
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
                fc.resize(pp.size());
                forAll(pp, fI)
                {
                    const face& f = pp[fI];
                    std::vector<double> c(2*ncoef, 0.0);
                    forAll(f, k)
                    {
                        const point& q = p0[f[k]];
                        auto it = table.find
                        (
                            std::make_pair(long(std::llround(q.x()*%(key)r)), long(std::llround(q.y()*%(key)r)))
                        );
                        if (it == table.end())
                        {
                            FatalErrorInFunction << "No replay row for " << q << exit(FatalError);
                        }
                        for (label j = 0; j < 2*ncoef; j++) { c[j] += it->second[j]/f.size(); }
                    }
                    fc[fI] = c;
                }
                loaded = true;
            }
            const volScalarField& p = m.lookupObject<volScalarField>("p");
            const volVectorField& U = m.lookupObject<volVectorField>("U");
            const tmp<volTensorField> tgU(fvc::grad(U));
            const tensorField& gUb = tgU().boundaryField()[pid];
            const symmTensorField devR(-%(nu)r*dev(twoSymm(gUb)));
            const scalarField& pb = p.boundaryField()[pid];
            const vectorField& Ub = U.boundaryField()[pid];
            const vectorField& Sf = m.Sf().boundaryField()[pid];
            const label H = ncoef/2;
            // basis values at t and t - dt for each harmonic: K s cos/sin(h w T)
            const scalar Tn = max(Tf(t), 0.0), To = max(Tf(t - dt), 0.0);
            const scalar gn = s(t)*Kf(t), go = s(t - dt)*Kf(t - dt);
            List<scalar> dW(2*H + 2, 0.0);   // h=1..H x,y and the mean x,y
            List<scalar> Q(10, 0.0);
            scalar Pp = 0, Pv = 0;
            vector Fp(Zero), Fv(Zero);
            forAll(pb, i)
            {
                const vector fp = %(rho)r*pb[i]*Sf[i]/%(depth)r;
                const vector fv = %(rho)r*(Sf[i] & devR[i])/%(depth)r;
                const vector f = fp + fv;
                Pp += fp & Ub[i];
                Pv += fv & Ub[i];
                Fp += fp;
                Fv += fv;
                const std::vector<double>& c = fc[i];
                for (label comp = 0; comp < 2; comp++)
                {
                    const scalar fcmp = f[comp];
                    dW[2*H + comp] += fcmp*c[comp*ncoef]*(gn - go);
                    for (label h = 1; h <= H; h++)
                    {
                        const scalar a = c[comp*ncoef + 2*h - 1], b = c[comp*ncoef + 2*h];
                        const scalar dn = gn*(a*std::cos(h*omega*Tn) + b*std::sin(h*omega*Tn));
                        const scalar d0 = go*(a*std::cos(h*omega*To) + b*std::sin(h*omega*To));
                        dW[2*(h - 1) + comp] += fcmp*(dn - d0);
                        if (h <= 2)
                        {
                            Q[4*(h - 1) + 2*comp] += fcmp*a;
                            Q[4*(h - 1) + 2*comp + 1] += fcmp*b;
                        }
                    }
                }
                Q[8] += fp[1]*c[ncoef + 1];
                Q[9] += fp[1]*c[ncoef + 2];
            }
            reduce(Pp, sumOp<scalar>());
            reduce(Pv, sumOp<scalar>());
            reduce(Fp, sumOp<vector>());
            reduce(Fv, sumOp<vector>());
            forAll(dW, k) { reduce(dW[k], sumOp<scalar>()); }
            forAll(Q, k) { reduce(Q[k], sumOp<scalar>()); }
            if (Pstream::master())
            {
                std::ofstream os((m.time().globalPath()/"energy.dat").c_str(), std::ios::app);
                os.precision(12);
                os << t << " " << Tn << " " << gn << " " << Pp << " " << Pv << " "
                   << Fp.x() << " " << Fp.y() << " " << Fv.x() << " " << Fv.y();
                forAll(dW, k) { os << " " << dW[k]; }
                forAll(Q, k) { os << " " << Q[k]; }
                os << "\n";
            }
        #};
    }
"""


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--source", required=True, type=Path)
    ap.add_argument("--restart", required=True)
    ap.add_argument("--end", required=True, type=float)
    ap.add_argument("--name", required=True)
    ap.add_argument("--work", type=Path, default=Path.home() / "ht_energy/work")
    ap.add_argument("--amplitude-factor", type=float, default=1.0)
    ap.add_argument("--frequency-factor", type=float, default=1.0)
    ap.add_argument("--blend", type=float, default=0.3)
    ap.add_argument("--setup-only", action="store_true")
    a = ap.parse_args()
    src, case = a.source.resolve(), (a.work / a.name).resolve()
    if case.exists():
        shutil.rmtree(case)
    case.mkdir(parents=True)
    procs = sorted(src.glob("processor*"))
    for d in ("constant", "system"):
        shutil.copytree(src / d, case / d)
    if procs:
        for pd in procs:
            (case / pd.name).mkdir()
            shutil.copytree(pd / "constant", case / pd.name / "constant")
            shutil.copytree(pd / a.restart, case / pd.name / a.restart)
        tdirs = [case / pd.name / a.restart for pd in procs]
    else:
        shutil.copytree(src / a.restart, case / a.restart)
        tdirs = [case / a.restart]
    tb = float(a.restart)
    bc_old = "            forAll(pts, i) { ref[i] = pts[i]; }"
    bc_new = ("            {\n"
              "                const polyMesh& pm = this->patch().boundaryMesh().mesh()();\n"
              "                const pointIOField p0(IOobject(\"points\", pm.time().constant(),"
              " pm.dbDir()/polyMesh::meshSubDir, pm.time(), IOobject::MUST_READ,"
              " IOobject::NO_WRITE, false));\n"
              "                const labelList& mp = this->patch().meshPoints();\n"
              "                forAll(pts, i) { ref[i] = p0[mp[i]]; }\n"
              "            }")
    m = re.search(r"const scalar t0 = ([-\d.eE+]+), ramp = ([-\d.eE+]+);", (tdirs[0] / "pointMotionU").read_text())
    t0, ramp = float(m.group(1)), float(m.group(2))
    mod = MOD % dict(k=a.amplitude_factor, r=a.frequency_factor, tb=tb, tr=a.blend)
    for td in tdirs:
        f = td / "pointMotionU"
        s = f.read_text()
        assert s.count(bc_old) == 1, f
        s = s.replace(bc_old, bc_new)
        anchor = "        auto s = [&](scalar tt) -> scalar"
        assert s.count(anchor) == 1
        s = s.replace(anchor, mod + anchor)
        old_t = "const scalar tn = max(t - t0, 0.0), to = max(t - dt - t0, 0.0);"
        assert s.count(old_t) == 1
        s = s.replace(old_t, "const scalar tn = max(Tf(t), 0.0), to = max(Tf(t - dt), 0.0);")
        s = s.replace("d1[comp] = s(t)*a;", "d1[comp] = s(t)*Kf(t)*a;")
        s = s.replace("d0[comp] = s(t - dt)*b;", "d0[comp] = s(t - dt)*Kf(t - dt)*b;")
        assert "Kf(t - dt)*b" in s
        f.write_text(s)
    fo = FO % dict(t0=t0, ramp=ramp, key=KEY, nu=NU, rho=RHO, depth=DEPTH,
                   mod=mod.replace("\n", "\n    "))
    fn = case / "system/functions"
    s = fn.read_text()
    i = s.rindex("}")
    fn.write_text(s[:i] + fo + s[i:])
    cd = case / "system/controlDict"
    s = cd.read_text()
    s = re.sub(r"(?m)^startFrom .*$", "startFrom       latestTime;", s)
    s = re.sub(r"(?m)^endTime .*$", f"endTime         {a.end:g};", s)
    s = re.sub(r"(?m)^writeInterval .*$", "writeInterval   100000000;", s)
    cd.write_text(s)
    (case / "energy_setup.txt").write_text(
        f"source {src}\nrestart {a.restart}\nend {a.end}\namplitude_factor {a.amplitude_factor}\n"
        f"frequency_factor {a.frequency_factor}\nblend {a.blend}\nt0 {t0}\nramp {ramp}\nnprocs {len(procs) or 1}\n")
    print(case, "nprocs", len(procs) or 1)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
