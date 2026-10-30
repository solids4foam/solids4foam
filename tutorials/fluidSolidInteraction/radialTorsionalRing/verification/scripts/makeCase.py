#!/usr/bin/env python3
"""Portable case generator for the radialTorsionalRing tutorial.

Writes every dictionary of the coupled radial-torsional ring case (full 2-D annulus,
solid A <= R <= B, fluid B <= r <= c, conformal interface at r = B): blockMeshDicts for
the three verification meshes, drive tables, coupling and solver settings, the coded
boundary conditions and the rtDiagFluid/rtDiagSolid conservation function objects.
The text of the dictionaries is identical to what the study's accepted runs used
(see the tutorial README, "Provenance").

Meshes (ntheta x nrs x nrf):  M1 64x8x16,  M2 128x16x32,  M3 256x32x64.
Configurations:               torsion (theta0=0.1, P0=0), radial (theta0=0, P0=-3000 Pa),
                              combined (theta0=0.1, P0=-3000 Pa).
Per-mesh outerCorrAbsoluteTolerance (explicit setting of the accepted runs): 1e-12 for
M1/M2, 1e-10 for M3; override with --absTol.
nNonOrthogonalCorrectors: 2 (the verification setting; 0 fails on M3 torsion near 1.81 s).

Usage (needs numpy, scipy for the drive tables and the reference):
  makeCase.py CASE_DIR [--mesh M1|M2|M3] [--config torsion|radial|combined]
                       [--np N] [--tend 2] [--nNonOrth 2] [--absTol X] [--no-reference]
Only the standard library, numpy and scipy are used; no machine paths.
"""
import sys as _sys
_sys.dont_write_bytecode = True
import argparse
import dataclasses
import json
import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import rt_reference as rt  # noqa: E402

MESHES = {"M1": (64, 8, 16), "M2": (128, 16, 32), "M3": (256, 32, 64)}
CONFIGS = {"torsion": (0.1, 0.0), "radial": (0.0, -3000.0), "combined": (0.1, -3000.0)}
ABSTOL = {"M1": 1e-12, "M2": 1e-12, "M3": 1e-10}


HEADER = """FoamFile
{{
    version     2.0;
    format      ascii;
    class       {cls};
    object      {obj};
}}
"""


def w(path, cls, obj, body):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as fh:
        fh.write(HEADER.format(cls=cls, obj=obj) + body)


def block_mesh(p, nth, nr, dz=0.1):
    """Four 90-degree blocks; vertices on the circles; arc edges."""
    A, B = p.A, p.B
    verts = []
    for z in (0.0, dz):
        for k in range(4):
            th = k * math.pi / 2
            verts.append((A * math.cos(th), A * math.sin(th), z))
            verts.append((B * math.cos(th), B * math.sin(th), z))
    # vertex index: (z, k, radial) -> 8*z + 2*k + r
    vi = lambda z, k, r: 8 * z + 2 * (k % 4) + r  # noqa: E731
    blocks, edges = [], []
    for k in range(4):
        h = [vi(0, k, 0), vi(0, k, 1), vi(0, k + 1, 1), vi(0, k + 1, 0),
             vi(1, k, 0), vi(1, k, 1), vi(1, k + 1, 1), vi(1, k + 1, 0)]
        blocks.append(f"    hex ({' '.join(map(str, h))}) ({nr} {nth // 4} 1) simpleGrading (1 1 1)")
        thm = (k + 0.5) * math.pi / 2
        for z, zz in ((0, 0.0), (1, dz)):
            for r, RR in ((0, A), (1, B)):
                edges.append(f"    arc {vi(z, k, r)} {vi(z, k + 1, r)} ({RR * math.cos(thm):.15g} {RR * math.sin(thm):.15g} {zz})")
    inner = [f"            ({vi(0, k, 0)} {vi(1, k, 0)} {vi(1, k + 1, 0)} {vi(0, k + 1, 0)})" for k in range(4)]
    outer = [f"            ({vi(0, k, 1)} {vi(0, k + 1, 1)} {vi(1, k + 1, 1)} {vi(1, k, 1)})" for k in range(4)]
    fb = [f"            ({vi(0, k, 0)} {vi(0, k + 1, 0)} {vi(0, k + 1, 1)} {vi(0, k, 1)})" for k in range(4)] + \
         [f"            ({vi(1, k, 0)} {vi(1, k, 1)} {vi(1, k + 1, 1)} {vi(1, k + 1, 0)})" for k in range(4)]
    vtxt = "\n".join(f"    ({x:.15g} {y:.15g} {z:.15g})" for x, y, z in verts)
    return f"""convertToMeters 1;
vertices
(
{vtxt}
);
blocks
(
{chr(10).join(blocks)}
);
edges
(
{chr(10).join(edges)}
);
boundary
(
    inner
    {{
        type patch;
        faces
        (
{chr(10).join(inner)}
        );
    }}
    outer
    {{
        type patch;
        faces
        (
{chr(10).join(outer)}
        );
    }}
    frontAndBack
    {{
        type empty;
        faces
        (
{chr(10).join(fb)}
        );
    }}
);
"""


def table(rows):
    return "(\n" + "\n".join(f"    ({t:.12g} {v:.17g})" for t, v in rows) + "\n)\n"


FLUID_DIAG = r'''
    rtDiag
    {
        type            coded;
        libs            (utilityFunctionObjects);
        region          fluid;
        name            rtDiagFluid;
        executeControl  timeStep;
        writeControl    timeStep;
        codeInclude
        #{
            #include "fvc.H"
            #include <fstream>
        #};
        codeExecute
        #{
            const scalar rho = RHO;
            const scalar nu = NU;
            const fvMesh& m = refCast<const fvMesh>(obr_);
            const volVectorField& U = m.lookupObject<volVectorField>("U");
            const volScalarField& p = m.lookupObject<volScalarField>("p");
            const surfaceScalarField& phi = m.lookupObject<surfaceScalarField>("phi");
            const volTensorField gradU(fvc::grad(U));
            const label iI = m.boundaryMesh().findPatchID("interface");
            const label iO = m.boundaryMesh().findPatchID("outer");
            // in-plane unit vectors from the x,y components of the face centre
            scalar A = 0, Rm = 0, pm = 0, srr = 0, snv = 0, srt = 0, Tq = 0, Pw = 0, wB = 0, ubn = 0;
            {
                const vectorField& Sf = m.Sf().boundaryField()[iI];
                const scalarField& mS = m.magSf().boundaryField()[iI];
                const vectorField& Cf = m.Cf().boundaryField()[iI];
                const tensorField& gb = gradU.boundaryField()[iI];
                const scalarField& pb = p.boundaryField()[iI];
                const vectorField& Ub = U.boundaryField()[iI];
                forAll(Cf, i)
                {
                    const symmTensor dRR(-rho*nu*dev(twoSymm(gb[i])));
                    const vector fP(rho*Sf[i]*pb[i]);
                    const vector fV(Sf[i] & dRR);
                    const scalar rxy = Foam::sqrt(sqr(Cf[i].x()) + sqr(Cf[i].y()));
                    const vector er(Cf[i].x()/rxy, Cf[i].y()/rxy, 0);
                    const vector et(-er.y(), er.x(), 0);
                    A += mS[i];
                    Rm += rxy*mS[i];
                    pm += rho*pb[i]*mS[i];
                    srr += (er & (fP + fV));
                    snv += (er & fV);
                    srt += (et & fV);
                    Tq += (Cf[i] ^ (fP + fV)).z();
                    Pw += ((fP + fV) & Ub[i]);
                    wB += (et & Ub[i])*mS[i];
                    ubn += (er & Ub[i])*mS[i];
                }
            }
            scalar Hf = 0, KEf = 0, Df = 0;
            forAll(U, c)
            {
                const vector Uc = U[c];
                Hf += rho*((m.C()[c] ^ Uc).z())*m.V()[c];
                KEf += 0.5*rho*magSqr(Uc)*m.V()[c];
                const symmTensor Dc(symm(gradU[c]));
                Df += 2.0*rho*nu*(Dc && Dc)*m.V()[c];
            }
            scalar Fc = 0, Tc = 0, Qo = 0, Pc = 0;
            {
                const vectorField& Sf = m.Sf().boundaryField()[iO];
                const vectorField& Cf = m.Cf().boundaryField()[iO];
                const tensorField& gb = gradU.boundaryField()[iO];
                const vectorField& Ub = U.boundaryField()[iO];
                const scalarField& pb = p.boundaryField()[iO];
                const scalarField& ph = phi.boundaryField()[iO];
                forAll(Cf, i)
                {
                    const symmTensor dRR(-rho*nu*dev(twoSymm(gb[i])));
                    const vector fOn(-rho*pb[i]*Sf[i] - (Sf[i] & dRR));   // force by exterior on fluid
                    Fc += rho*((Cf[i] ^ Ub[i]).z())*ph[i];
                    Tc += (Cf[i] ^ (-(Sf[i] & dRR))).z();
                    Qo += ph[i];
                    Pc += (fOn & Ub[i]) - 0.5*rho*magSqr(Ub[i])*ph[i];   // power in, minus KE convected out
                }
            }
            // parallel: global sums (serial: no-op)
            reduce(A, sumOp<scalar>()); reduce(Rm, sumOp<scalar>()); reduce(pm, sumOp<scalar>());
            reduce(srr, sumOp<scalar>()); reduce(snv, sumOp<scalar>()); reduce(srt, sumOp<scalar>());
            reduce(Tq, sumOp<scalar>()); reduce(Pw, sumOp<scalar>()); reduce(wB, sumOp<scalar>());
            reduce(ubn, sumOp<scalar>()); reduce(Hf, sumOp<scalar>()); reduce(KEf, sumOp<scalar>());
            reduce(Df, sumOp<scalar>()); reduce(Fc, sumOp<scalar>()); reduce(Tc, sumOp<scalar>());
            reduce(Qo, sumOp<scalar>()); reduce(Pc, sumOp<scalar>());
            if (!Pstream::master()) return true;
            std::ofstream os((m.time().globalPath()/"rtDiagFluid.dat").c_str(), std::ios::app);
            os.precision(16);
            os << m.time().value() << ' ' << A << ' ' << Rm/A << ' ' << pm/A << ' ' << srr/A << ' '
               << snv/A << ' ' << srt/A << ' ' << Tq << ' ' << Pw << ' ' << Hf << ' ' << KEf << ' '
               << Df << ' ' << Fc << ' ' << Tc << ' ' << Qo << ' ' << Pc << ' ' << wB/A << ' ' << ubn/A
               << std::endl;
        #};
    }
'''

SOLID_DIAG = r'''
    rtSolid
    {
        type            coded;
        libs            (utilityFunctionObjects);
        region          solid;
        name            rtDiagSolid;
        executeControl  timeStep;
        writeControl    timeStep;
        codeInclude
        #{
            #include "fvc.H"
            #include <fstream>
        #};
        codeExecute
        #{
            const scalar rho = RHOS;
            const scalar mu = MUS;
            const scalar kappa = KAPPAS;
            const fvMesh& m = refCast<const fvMesh>(obr_);
            const volVectorField& D = m.lookupObject<volVectorField>("D");
            const scalar dt = m.time().deltaTValue();
            const volVectorField& D0 = D.oldTime();
            const volVectorField& D00 = D.oldTime().oldTime();
            const volTensorField gD(fvc::grad(D));
            // interface (R = B) and inner (R = A) patches: radius, rotation, n = 2 Fourier modes
            scalar out[12] = {0,0,0,0,0,0,0,0,0,0,0,0};
            scalar hi[4] = {0,0,0,0};
            const char* pn[2] = {"interface", "inner"};
            for (label k = 0; k < 2; ++k)
            {
                const label pI = m.boundaryMesh().findPatchID(pn[k]);
                const vectorField& Cf = m.Cf().boundaryField()[pI];
                const scalarField& mS = m.magSf().boundaryField()[pI];
                const vectorField& Db = D.boundaryField()[pI];
                scalar A = 0, ur = 0, al = 0, c2 = 0, s2 = 0, R0 = 0, c3 = 0, s3 = 0, c4 = 0, s4 = 0;
                forAll(Cf, i)
                {
                    const scalar X = Cf[i].x(), Y = Cf[i].y();
                    const scalar R = Foam::sqrt(X*X + Y*Y);
                    const scalar th = Foam::atan2(Y, X);
                    const scalar x = X + Db[i].x(), y = Y + Db[i].y();
                    const scalar r = Foam::sqrt(x*x + y*y);
                    scalar da = Foam::atan2(y, x) - th;
                    if (da > constant::mathematical::pi) da -= constant::mathematical::twoPi;
                    if (da < -constant::mathematical::pi) da += constant::mathematical::twoPi;
                    A += mS[i]; ur += (r - R)*mS[i]; al += da*mS[i]; R0 += R*mS[i];
                    c2 += (r - R)*Foam::cos(2*th)*mS[i]; s2 += (r - R)*Foam::sin(2*th)*mS[i];
                    c3 += (r - R)*Foam::cos(3*th)*mS[i]; s3 += (r - R)*Foam::sin(3*th)*mS[i];
                    c4 += (r - R)*Foam::cos(4*th)*mS[i]; s4 += (r - R)*Foam::sin(4*th)*mS[i];
                }
                reduce(A, sumOp<scalar>()); reduce(ur, sumOp<scalar>()); reduce(al, sumOp<scalar>());
                reduce(c2, sumOp<scalar>()); reduce(s2, sumOp<scalar>()); reduce(R0, sumOp<scalar>());
                reduce(c3, sumOp<scalar>()); reduce(s3, sumOp<scalar>()); reduce(c4, sumOp<scalar>()); reduce(s4, sumOp<scalar>());
                if (k == 0) { hi[0] = 2*c3/A; hi[1] = 2*s3/A; hi[2] = 2*c4/A; hi[3] = 2*s4/A; }
                out[6*k + 0] = R0/A; out[6*k + 1] = ur/A; out[6*k + 2] = al/A;
                out[6*k + 3] = 2*c2/A; out[6*k + 4] = 2*s2/A; out[6*k + 5] = A;
            }
            scalar Hs = 0, KEs = 0, Ws = 0, lmin = GREAT, Jmin = GREAT, Jmax = 0;
            forAll(D, c)
            {
                const vector V((1.5*D[c] - 2.0*D0[c] + 0.5*D00[c])/dt);
                const vector x(m.C()[c] + D[c]);
                Hs += rho*((x ^ V).z())*m.V()[c];
                KEs += 0.5*rho*magSqr(V)*m.V()[c];
                tensor F(I + gD[c].T());
                F.zz() = 1.0; F.xz() = 0; F.yz() = 0; F.zx() = 0; F.zy() = 0;
                const scalar J = det(F);
                const symmTensor C(symm(F.T() & F));
                const scalar trC = tr(C);
                Ws += (0.5*mu*(Foam::pow(J, -2.0/3.0)*trC - 3.0)
                     + 0.25*kappa*(J*J - 1.0 - 2.0*Foam::log(J)))*m.V()[c];
                // in-plane principal stretches of C
                const scalar a = C.xx(), b = C.yy(), d = C.xy();
                const scalar lmn = 0.5*(a + b) - Foam::sqrt(0.25*sqr(a - b) + d*d);
                lmin = min(lmin, Foam::sqrt(max(lmn, 0.0)));
                Jmin = min(Jmin, J); Jmax = max(Jmax, J);
            }
            reduce(Hs, sumOp<scalar>()); reduce(KEs, sumOp<scalar>()); reduce(Ws, sumOp<scalar>());
            reduce(lmin, minOp<scalar>()); reduce(Jmin, minOp<scalar>()); reduce(Jmax, maxOp<scalar>());
            if (!Pstream::master()) return true;
            std::ofstream os((m.time().globalPath()/"rtDiagSolid.dat").c_str(), std::ios::app);
            os.precision(16);
            os << m.time().value();
            for (label k = 0; k < 12; ++k) os << ' ' << out[k];
            os << ' ' << Hs << ' ' << KEs << ' ' << Ws << ' ' << lmin << ' ' << Jmin << ' ' << Jmax
               << ' ' << hi[0] << ' ' << hi[1] << ' ' << hi[2] << ' ' << hi[3] << std::endl;
        #};
    }
'''


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("case")
    ap.add_argument("--mesh", choices=sorted(MESHES), default=None,
                    help="M1 64x8x16, M2 128x16x32, M3 256x32x64 (ntheta x nrs x nrf); default M1")
    ap.add_argument("--config", choices=sorted(CONFIGS), default=None,
                    help="torsion (theta0=0.1,P0=0), radial (theta0=0,P0=-3000), combined; default combined")
    ap.add_argument("--nNonOrth", type=int, default=2, help="PIMPLE nNonOrthogonalCorrectors (verification setting 2)")
    ap.add_argument("--no-reference", dest="reference", action="store_false",
                    help="do not write reference.json (the drive tables are always written)")
    ap.add_argument("--theta0", type=float, default=None)
    ap.add_argument("--Omega", type=float, default=3.0)
    ap.add_argument("--P0", type=float, default=None)
    ap.add_argument("--Tramp", type=float, default=1.0)
    ap.add_argument("--ntheta", type=int, default=None)
    ap.add_argument("--nrs", type=int, default=None)
    ap.add_argument("--nrf", type=int, default=None)
    ap.add_argument("--dt", type=float, default=0.0025)
    ap.add_argument("--tend", type=float, default=2.0)
    ap.add_argument("--write", type=float, default=0.1)
    ap.add_argument("--seed", type=float, default=0.0,
                    help="n=2 outer-pressure pulse amplitude [Pa], C-inf bump on 0.1..0.5 s (stability screen)")
    ap.add_argument("--seedModes", default="2", help="comma list of lobe numbers seeded together, e.g. 2,3,4")
    ap.add_argument("--outerTol", type=float, default=1e-6)
    ap.add_argument("--absTol", type=float, default=None,
                    help="outerCorrAbsoluteTolerance; default per mesh: 1e-12 (M1, M2), 1e-10 (M3)")
    ap.add_argument("--relax", type=float, default=0.05, help="IQN-ILS initial relaxation factor")
    ap.add_argument("--reuse", type=int, default=5, help="IQN-ILS couplingReuse (previous steps whose secants are reused)")
    ap.add_argument("--np", type=int, default=4, help="subdomains written to decomposeParDict (used only if run in parallel)")
    a = ap.parse_args()
    a.mesh = a.mesh or "M1"
    a.config = a.config or "combined"
    for k, v in zip(("ntheta", "nrs", "nrf"), MESHES[a.mesh]):
        if getattr(a, k) is None:
            setattr(a, k, v)
    if a.theta0 is None:
        a.theta0 = CONFIGS[a.config][0]
    if a.P0 is None:
        a.P0 = CONFIGS[a.config][1]
    if a.absTol is None:
        a.absTol = ABSTOL[a.mesh]
    p = rt.Params(theta0=a.theta0, Omega=a.Omega, P0=a.P0, T_ramp=a.Tramp, ramp="Cinf")
    C = a.case
    n = int(round(a.tend / a.dt))
    tt = np.arange(n + 3) * a.dt
    os.makedirs(f"{C}/constant", exist_ok=True)
    if a.reference:
        import make_reference
        make_reference.write_reference(C, p, a.ntheta, a.dt, a.tend, vars(a))
    with open(f"{C}/constant/timeVsAngle", "w") as fh:
        fh.write(table([(t, math.degrees(rt.rotation_drive(p, t)[0])) for t in tt]))
    with open(f"{C}/constant/pcTable", "w") as fh:
        fh.write(table([(t, rt.pressure_drive(p, t) / p.rho_f) for t in tt]))
    # ------------------------------------------------------------- meshes
    bs = block_mesh(p, a.ntheta, a.nrs).replace("outer", "interface")
    bs = bs.replace("interface\n    {\n        type patch;", "interface\n    {\n        type wall;")
    w(f"{C}/system/solid/blockMeshDict", "dictionary", "blockMeshDict", bs)
    q = dataclasses.replace(p, A=p.B, B=p.c)
    bf = block_mesh(q, a.ntheta, a.nrf).replace("inner", "interface")
    bf = bf.replace("interface\n    {\n        type patch;", "interface\n    {\n        type wall;")
    assert bs.count("type wall;") == 1 and bf.count("type wall;") == 1
    w(f"{C}/system/fluid/blockMeshDict", "dictionary", "blockMeshDict", bf)
    # ------------------------------------------------------------- global
    w(f"{C}/constant/g", "uniformDimensionedVectorField", "g", "dimensions [0 1 -2 0 0 0 0];\nvalue (0 0 0);\n")
    w(f"{C}/constant/physicsProperties", "dictionary", "physicsProperties", "type    fluidSolidInteraction;\n")
    w(f"{C}/constant/fsiProperties", "dictionary", "fsiProperties", f"""fluidSolidInterface    IQNILS;
IQNILSCoeffs
{{
    solidPatch interface;
    fluidPatch interface;
    predictor            yes;
    predictSolid         yes;
    relaxationFactor     {a.relax};
    outerCorrTolerance   {a.outerTol};
    // absolute interface point-mismatch tolerance [m]: the C-inf ramp starts with signals at
    // round-off level, where the relative residual (normalised by the largest residual seen so
    // far) cannot reach 1e-6 (pilot R1 attempt 1 stalled at t = 0.035 s); 1e-12 m is 1e-11 of the
    // expected interface displacement
    outerCorrAbsoluteTolerance {a.absTol};
    nOuterCorr           100;
    couplingReuse        {a.reuse};
    coupled              yes;
    interfaceTransferMethod directMap;
    writeResidualsToFile yes;
}}
""")
    diag = (FLUID_DIAG.replace("RHO", str(p.rho_f)).replace("NU", str(p.nu_f))
            + SOLID_DIAG.replace("RHOS", str(p.rho_s)).replace("MUS", repr(p.G)).replace("KAPPAS", repr(p.kappa)))
    w(f"{C}/system/controlDict", "dictionary", "controlDict", f"""application     solids4Foam;
startFrom       startTime;
startTime       0;
stopAt          endTime;
endTime         {a.tend};
deltaT          {a.dt};
writeControl    adjustableRunTime;
writeInterval   {a.write};
purgeWrite      0;
writeFormat     ascii;
writePrecision  14;
writeCompression off;
timeFormat      general;
timePrecision   10;
runTimeModifiable no;
functions
{{
{diag}
}}
""")
    # ------------------------------------------------------------- solid
    w(f"{C}/constant/solid/dynamicMeshDict", "dictionary", "dynamicMeshDict", "dynamicFvMesh   staticFvMesh;\n")
    w(f"{C}/constant/solid/g", "uniformDimensionedVectorField", "g", "dimensions [0 1 -2 0 0 0 0];\nvalue (0 0 0);\n")
    w(f"{C}/constant/solid/solidProperties", "dictionary", "solidProperties", """solidModel     nonLinearGeometryTotalLagrangianTotalDisplacement;
nonLinearGeometryTotalLagrangianTotalDisplacementCoeffs
{
    solutionAlgorithm PETScSNES;
    predictor yes;
    stabilisation
    {
        momentum
        {
            type        diffStencilLaplacian;
            scaleFactor 0.5;
        }
    }
}
""")
    w(f"{C}/constant/solid/mechanicalProperties", "dictionary", "mechanicalProperties", f"""planeStress     no;
mechanical
(
    ring
    {{
        type            neoHookeanElastic;
        rho             rho [1 -3 0 0 0 0 0] {p.rho_s};
        E               E [1 -1 -2 0 0 0 0] {p.E};
        nu              nu [0 0 0 0 0 0 0] {p.nu_s};
    }}
);
""")
    w(f"{C}/system/solid/fvSchemes", "dictionary", "fvSchemes", """d2dt2Schemes { default backward; }
ddtSchemes { default backward; }
gradSchemes { default leastSquaresS4f; }
divSchemes { default Gauss linear; }
laplacianSchemes { default Gauss linear corrected; }
snGradSchemes { default corrected; }
interpolationSchemes { default linear; }
""")
    w(f"{C}/system/solid/fvSolution", "dictionary", "fvSolution", """solvers
{
    D
    {
        solver petsc;
        options
        {
            snes_type newtonls;
            snes_linesearch_type basic;
            snes_rtol "1e-10";
            snes_atol "1e-11";
            snes_stol "1e-12";
            snes_max_it "200";
            snes_converged_reason;
            snes_mf;
            snes_mf_operator;
            ksp_type lgmres;
            ksp_gmres_restart "100";
            snes_lag_preconditioner_persists true;
            snes_lag_preconditioner "-2";
            pc_type bjacobi;
            sub_pc_type lu;
        }
    }
}
relaxationFactors {}
""")
    w(f"{C}/0/solid/D", "volVectorField", "D", """dimensions      [0 1 0 0 0 0 0];
internalField   uniform (0 0 0);
boundaryField
{
    inner
    {
        type            fixedRotation;
        rotationAngleSeries
        {
            file            "$FOAM_CASE/constant/timeVsAngle";
            outOfBounds     clamp;
        }
        rotationAxis    (0 0 1);
        rotationOrigin  (0 0 0);
        value           uniform (0 0 0);
    }
    interface
    {
        type            solidTraction;
        traction        uniform (0 0 0);
        pressure        uniform 0;
        value           uniform (0 0 0);
    }
    frontAndBack
    {
        type            empty;
    }
}
""")
    # ------------------------------------------------------------- fluid
    w(f"{C}/constant/fluid/fluidProperties", "dictionary", "fluidProperties", "fluidModel pimpleFluid;\npimpleFluidCoeffs {}\n")
    w(f"{C}/constant/fluid/transportProperties", "dictionary", "transportProperties",
      f"transportModel Newtonian;\nnu nu [0 2 -1 0 0 0 0] {p.nu_f};\nrho rho [1 -3 0 0 0 0 0] {p.rho_f};\n")
    for nm in ("turbulenceProperties", "momentumTransport"):
        w(f"{C}/constant/fluid/{nm}", "dictionary", nm, "simulationType laminar;\n")
    w(f"{C}/constant/fluid/g", "uniformDimensionedVectorField", "g", "dimensions [0 1 -2 0 0 0 0];\nvalue (0 0 0);\n")
    w(f"{C}/constant/fluid/dynamicMeshDict", "dictionary", "dynamicMeshDict", """dynamicFvMesh   dynamicMotionSolverFvMesh;
motionSolverLibs (fvMotionSolvers);
motionSolver    displacementLaplacian;
displacementLaplacianCoeffs
{
    diffusivity     quadratic inverseDistance (interface);
}
diffusivity     quadratic inverseDistance (interface);
""")
    w(f"{C}/system/fluid/fvSchemes", "dictionary", "fvSchemes", """ddtSchemes { default backward; }
gradSchemes { default leastSquares; }
divSchemes
{
    default         none;
    div(phi,U)      Gauss linear;
    div((nuEff*dev2(T(grad(U))))) Gauss linear;
    div((nuEff*dev(T(grad(U))))) Gauss linear;
}
laplacianSchemes { default Gauss linear corrected; }
interpolationSchemes { default linear; }
snGradSchemes { default corrected; }
fluxRequired { default no; p; }
""")
    w(f"{C}/system/fluid/fvSolution", "dictionary", "fvSolution", """solvers
{
    "p|pFinal"
    {
        solver          GAMG;
        tolerance       1e-12;
        relTol          0;
        minIter         1;
        maxIter         1000;
        smoother        GaussSeidel;
    }
    "cellDisplacement|cellDisplacementFinal"
    {
        solver          GAMG;
        tolerance       1e-13;
        relTol          0;
        minIter         1;
        maxIter         200;
        smoother        GaussSeidel;
    }
    "U|UFinal"
    {
        solver          PBiCG;
        preconditioner  DILU;
        tolerance       1e-13;
        relTol          0;
        minIter         1;
    }
}
PIMPLE
{
    fluxConsistentPatches (outer);
    momentumPredictor   on;
    nOuterCorrectors    2;
    nCorrectors         3;
    nNonOrthogonalCorrectors NNONORTH;
}
""".replace("NNONORTH", str(a.nNonOrth)))
    w(f"{C}/0/fluid/U", "volVectorField", "U", """dimensions      [0 1 -1 0 0 0 0];
internalField   uniform (0 0 0);
boundaryField
{
    interface
    {
        type            newMovingWallVelocity;
        value           uniform (0 0 0);
    }
    outer
    {
        type            codedMixed;
        refValue        uniform (0 0 0);
        refGradient     uniform (0 0 0);
        valueFraction   uniform 0;
        value           uniform (0 0 0);
        name            openRadialSwirlFsi;
        code
        #{
            const vectorField nf(this->patch().nf());
            const vectorField Cf(this->patch().Cf());
            scalarField rf(Cf.size());
            forAll(Cf, i) { rf[i] = Foam::sqrt(sqr(Cf[i].x()) + sqr(Cf[i].y())); }
            const vectorField Up(*this);
            this->refGrad() = (Up - 2.0*nf*(nf & Up))/rf;
            this->refValue() = Up;
            this->valueFraction() = 0.0;
        #};
    }
    frontAndBack
    {
        type            empty;
    }
}
""")
    seedP = a.seed / p.rho_f
    w(f"{C}/0/fluid/p", "volScalarField", "p", f"""dimensions      [0 2 -2 0 0 0 0];
internalField   uniform 0;
boundaryField
{{
    interface
    {{
        type            zeroGradient;
    }}
    outer
    {{
        type            codedFixedValue;
        value           uniform 0;
        name            outerPressureSeeded;
        codeInclude
        #{{
            #include "interpolationTable.H"
        #}};
        code
        #{{
            // p_c(t)/rho from the reference table; optional n = 2 pulse (stability screen):
            // seed*cos(2 theta)*B(t), B = C-inf bump on [0.1, 0.5] s, peak 1 at 0.3 s
            const scalar t = this->db().time().value();
            interpolationTable<scalar> pTab(fileName("$FOAM_CASE/constant/pcTable").expand());
            const scalar pc = pTab(t);
            const scalar seed = {seedP!r};
            const labelList modes({{{",".join(str(int(x)) for x in a.seedModes.split(","))}}});
            scalar bump = 0;
            const scalar tau = (t - 0.1)/0.4;
            if (tau > 0 && tau < 1)
            {{
                bump = Foam::exp(1.0 - 1.0/(1.0 - sqr(2*tau - 1)));
            }}
            const vectorField Cf(this->patch().Cf());
            scalarField pv(Cf.size());
            forAll(Cf, i)
            {{
                const scalar th = Foam::atan2(Cf[i].y(), Cf[i].x());
                scalar pert = 0;
                forAll(modes, j) {{ pert += Foam::cos(modes[j]*th); }}
                pv[i] = pc + seed*bump*pert;
            }}
            this->operator==(pv);
        #}};
    }}
    frontAndBack
    {{
        type            empty;
    }}
}}
""")
    w(f"{C}/0/fluid/pointDisplacement", "pointVectorField", "pointDisplacement", """dimensions      [0 1 0 0 0 0 0];
internalField   uniform (0 0 0);
boundaryField
{
    interface
    {
        type            fixedValue;
        value           uniform (0 0 0);
    }
    outer
    {
        type            fixedValue;
        value           uniform (0 0 0);
    }
    frontAndBack
    {
        type            empty;
    }
}
""")
    for reg in ("fluid", "solid"):
        w(f"{C}/system/{reg}/decomposeParDict", "dictionary", "decomposeParDict",
          f"numberOfSubdomains {a.np};\nmethod          scotch;\n")
    # informational record of the generator inputs; the solver dictionaries are authoritative
    with open(f"{C}/caseSetup.json", "w") as fh:
        json.dump(dict(note="generator inputs; the solver dictionaries in constant/ and system/ are authoritative",
                       mesh=a.mesh, config=a.config, ntheta=a.ntheta, nrs=a.nrs, nrf=a.nrf, dt=a.dt, tend=a.tend,
                       theta0=a.theta0, P0=a.P0, Omega=a.Omega, Tramp=a.Tramp, outerCorrTolerance=a.outerTol,
                       outerCorrAbsoluteTolerance=a.absTol, relaxationFactor=a.relax, couplingReuse=a.reuse,
                       nNonOrthogonalCorrectors=a.nNonOrth, np=a.np, params=repr(p)), fh, indent=1)
        fh.write("\n")
    print("written", C)


if __name__ == "__main__":
    main()
