#!/usr/bin/env python3
"""Family-consistency validation, run BEFORE any spatial-convergence or mixed-term assessment.

WHAT IS CHECKED (everything is read from the actual solver dictionaries of each case, never from
case_settings.txt, reference.json args or command.txt; dictionaries are compared SEMANTICALLY: comments, whitespace,
number formats and the order of entries do not matter):
  physics        E, nu, rho of the solid; nu, rho of the fluid; deltaT, endTime; PIMPLE nOuterCorrectors,
                 nCorrectors, nNonOrthogonalCorrectors, momentumPredictor, fluxConsistentPatches;
  drives         constant/timeVsAngle and constant/pcTable recomputed with rt_reference.py (theta0, P0, Omega, ramp);
  mesh           for BOTH regions: ntheta x nrs x nrf and the refinement ratio (2); all vertex coordinates (inner radius A,
                 interface radius B, outer radius c, depth), the arc points, the block topology, the grading (1 1 1) and the
                 patch definitions of the blockMeshDicts;
  coupling       fluidSolidInterface type (IQN-ILS), interfaceTransferMethod, predictors, outerCorrTolerance,
                 relaxationFactor, couplingReuse, nOuterCorr, outerCorrAbsoluteTolerance, allowUnconvergedCoupling,
                 requireAllResidualMeasures;
  solid          solidModel, solution algorithm, stabilisation type and scaleFactor, material law type, planeStress;
  fluid          fluidModel, laminar, dynamicMeshDict (motion solver, diffusivity);
  controlDict    every top-level entry, before AND after the functions block (deltaT, adjustTimeStep, writeControl,
                 startFrom, ...), except endTime and writeInterval; the parameters of every function object (type,
                 libs, region, name, execute/write control) and the normalised code of rtDiagFluid and rtDiagSolid
                 against the documented current and older versions;
  boundary       0/fluid/U, p, pointDisplacement and 0/solid/D: patch types and the normalised text of the coded
                 conditions (case paths and comments ignored), the seed of the outer pressure condition;
  numerics       linear-solver settings, fvSchemes, the solid fvSolution;
  safety net     a semantic fingerprint (canonical form of the parsed dictionary) of every small dictionary of constant/ and
                 system/ against the documented benchmark set and against the other meshes.
NOT checked here: the solver binary, OpenFOAM/PETSc versions, the points of the generated polyMesh (the blockMeshDicts that
generate them are), the log files (layer A does that), and function-object code other than the two documented versions.

Statuses of the family (shared severity order in status.py):
  INCOMPLETE     a run or dictionary is missing, empty, unreadable or garbled;
  INVALID-INPUT  an inconsistent setting between the meshes, or a setting that differs from the documented benchmark
                 set (the bands then do not apply);
  CONDITIONAL    a documented-equivalent difference only: absolute coupling tolerance 1e-12 m (M1, M2) against 1e-10 m (M3),
                 the older diagnostics function objects, the older seeded-pressure code with seed 0 (identical
                 boundary values), different rank counts and scotch partitions;
  PASS           none of these differences.
A CONDITIONAL family can never reach a layer-B PASS (verify.py caps its claims at CONDITIONAL).
"""
import sys as _sys
_sys.dont_write_bytecode = True
import hashlib
import math
import os
import re

import numpy as np

import rt_reference as rt
import status as ST

MESHES = ("M1", "M2", "M3")
SEVERITY = ST.ORDER

DICTS = [  # small dictionaries fingerprinted as normalised text
    "constant/physicsProperties", "constant/g",
    "constant/fluid/fluidProperties", "constant/fluid/transportProperties", "constant/fluid/momentumTransport",
    "constant/fluid/turbulenceProperties", "constant/fluid/dynamicMeshDict", "constant/fluid/g",
    "constant/solid/solidProperties", "constant/solid/mechanicalProperties", "constant/solid/dynamicMeshDict", "constant/solid/g",
    "system/fluid/fvSchemes", "system/fluid/fvSolution", "system/solid/fvSchemes", "system/solid/fvSolution",
    "constant/fsiProperties",
    "0/fluid/U", "0/fluid/p", "0/fluid/pointDisplacement", "0/solid/D",
]


def _read(path):
    with open(path, errors="replace") as fh:
        txt = fh.read()
    if not txt.strip():
        raise ValueError(f"{os.path.basename(path)} is empty")
    return txt


def _num(txt, key, cast=float):
    m = re.search(r"^\s*" + re.escape(key) + r"\s+([^;\s]+)\s*;", txt, re.M)
    if not m:
        return None
    v = m.group(1)
    if cast is int:
        x = float(v)
        if x != int(x):
            raise ValueError(f"{key} is not an integer: {v}")
        return int(x)
    return cast(v)


def _tab(path):
    rows = re.findall(r"\(\s*([-+0-9.eE]+)\s+([-+0-9.eE]+)\s*\)", _read(path))
    a = np.array([[float(x), float(y)] for x, y in rows])
    if a.ndim != 2 or len(a) < 3:
        raise ValueError(f"{os.path.basename(path)} holds no usable table")
    return a


def _hexes(path):
    return [tuple(int(x) for x in m) for m in re.findall(r"hex \([^)]*\) \((\d+) (\d+) (\d+)\)", _read(path))]


def strip_comments(t):
    t = re.sub(r"/\*.*?\*/", " ", t, flags=re.S)
    return re.sub(r"//[^\n]*", " ", t)


def split_code(t):
    """Return (text with every #{ ... #} block replaced by a token, [normalised code blocks])."""
    blocks = []

    def sub(m):
        blocks.append(" ".join(strip_comments(m.group(1)).split()))
        return " #CODE# "
    t = re.sub(r"#\{(.*?)#\}", sub, t, flags=re.S)
    return t, blocks


def sha(s):
    return hashlib.sha256(s.encode()).hexdigest()[:16]


# ---------------------------------------------------------------- semantic dictionary parser
_NUM = re.compile(r"^[-+]?(\d+\.?\d*|\.\d+)([eE][-+]?\d+)?$")


def _canon_token(t):
    return format(float(t), ".12g") if _NUM.match(t) else t


def tokenize(text):
    """Tokens of an OpenFOAM dictionary; #{ ... #} code blocks become one ('CODE', hash) token."""
    text = strip_comments(text)
    text = re.sub(r"FoamFile\s*\{[^}]*\}", " ", text, count=1)
    out, i, n = [], 0, len(text)
    while i < n:
        c = text[i]
        if c.isspace():
            i += 1
        elif text.startswith("#{", i):
            j = text.find("#}", i)
            if j < 0:
                raise ValueError("unterminated #{ code block")
            out.append(("CODE", sha(" ".join(text[i + 2:j].split()))))
            i = j + 2
        elif c in "{};()":
            out.append(c)
            i += 1
        elif c == '"':
            j = text.find('"', i + 1)
            if j < 0:
                raise ValueError("unterminated string")
            out.append(text[i:j + 1])
            i = j + 1
        else:
            j = i
            while j < n and not text[j].isspace() and text[j] not in '{};()"' and not text.startswith("#{", j):
                j += 1
            out.append(_canon_token(text[i:j]))
            i = j
    return out


def parse_entries(tok, pos=0, closing=False, mask=()):
    """Parse entries until the matching '}' (or the end). Returns (list of (key, canonical value string), new pos)."""
    entries = []
    while pos < len(tok):
        t = tok[pos]
        if t == "}":
            if not closing:
                raise ValueError("unbalanced '}'")
            return entries, pos + 1
        if isinstance(t, tuple) or t in "{;()":
            raise ValueError(f"unexpected token {t!r}")
        key = t
        pos += 1
        if pos >= len(tok):
            raise ValueError(f"entry {key} has no value")
        if tok[pos] == "{":
            sub, pos = parse_entries(tok, pos + 1, True, mask)
            entries.append((key, "{" + " ".join(sorted(f"{k} {v}" for k, v in sub)) + "}"))
            continue
        vals, depth = [], 0
        while True:
            if pos >= len(tok):
                raise ValueError(f"entry {key} is not terminated by ';'")
            u = tok[pos]
            if u == ";" and depth == 0:
                pos += 1
                break
            if u == "(":
                depth += 1
            elif u == ")":
                depth -= 1
                if depth < 0:
                    raise ValueError("unbalanced ')'")
            if u == "{":
                sub, pos = parse_entries(tok, pos + 1, True, mask)
                vals.append("{" + " ".join(sorted(f"{k} {v}" for k, v in sub)) + "}")
                continue
            vals.append(("CODE#" + u[1]) if isinstance(u, tuple) else u)
            pos += 1
        entries.append((key, "X" if key in mask else " ".join(vals)))
    if closing:
        raise ValueError("missing '}'")
    return entries, pos


def parse_dict(text, mask=()):
    tok = tokenize(text)
    entries, _ = parse_entries(tok, 0, False, mask)
    return entries


def canon(entries, mask=(), code_mask=False):
    """Canonical string: entries sorted; values of keys in `mask` replaced by X; code hashes masked on request."""
    out = []
    for k, v in entries:
        if k in mask:
            v = "X"
        if code_mask:
            v = re.sub(r"CODE#[0-9a-f]+", "CODE", v)
        out.append(f"{k} {v}")
    return " | ".join(sorted(out))


def patch_types(text):
    """[(patch, type)] of the boundaryField of a field file (first-level `type` of every patch), in canonical order."""
    tok = tokenize(text)
    out, i = [], 0
    while i < len(tok) - 1:
        if tok[i] == "boundaryField" and tok[i + 1] == "{":
            i += 2
            depth = 1
            patch, pdepth = None, 0
            while i < len(tok) and depth > 0:
                t = tok[i]
                if t == "{":
                    depth += 1
                elif t == "}":
                    depth -= 1
                elif depth == 1 and not isinstance(t, tuple) and i + 1 < len(tok) and tok[i + 1] == "{":
                    patch = t
                elif depth == 2 and t == "type" and patch is not None and i + 1 < len(tok):
                    out.append((patch, tok[i + 1]))
                    patch = None
                i += 1
            break
        i += 1
    return sorted(out)


def code_hashes_of(entries):
    return [h for _, v in entries for h in re.findall(r"CODE#([0-9a-f]+)", v)]


def mesh_info(path):
    """Semantic content of a blockMeshDict: vertices, arcs, block topology and grading, patches, cell counts."""
    entries = dict(parse_dict(_read(path)))
    for k in ("vertices", "blocks", "edges", "boundary"):
        if k not in entries:
            raise ValueError(f"blockMeshDict has no {k}")
    vtx = [tuple(float(x) for x in m) for m in re.findall(r"\(\s*([-+0-9.e]+) ([-+0-9.e]+) ([-+0-9.e]+)\s*\)", entries["vertices"])]
    blocks = re.findall(r"hex \(((?: \d+)+) \) \( (\d+) (\d+) (\d+) \) simpleGrading \( ([^)]*)\)", entries["blocks"].replace("( ", "( ").replace("(", "( ").replace("  ", " "))
    if len(vtx) < 8 or not blocks:
        raise ValueError("blockMeshDict vertices or blocks could not be parsed")
    edges = re.findall(r"arc (\d+) (\d+) \( ([-+0-9.e]+) ([-+0-9.e]+) ([-+0-9.e]+) \)", entries["edges"].replace("(", "( ").replace("  ", " "))
    radii = sorted({round(math.hypot(v[0], v[1]), 9) for v in vtx})
    arcs = sorted({round(math.hypot(float(e[2]), float(e[3])), 9) for e in edges})
    topo = [(b[0].split(), b[4].split()) for b in blocks]
    geom = canon([("vertices", entries["vertices"]), ("edges", entries["edges"]), ("boundary", entries["boundary"]),
                  ("convertToMeters", dict(parse_dict(_read(path))).get("convertToMeters", "1")),
                  ("topology", repr(topo))])
    return dict(vertices=vtx, radii=radii, z=sorted({round(v[2], 12) for v in vtx}), arc_radii=arcs, grading=[tuple(b[4].split()) for b in blocks],
                counts=[tuple(int(x) for x in b[1:4]) for b in blocks], geometry_fp=sha(geom), boundary_fp=sha(entries["boundary"]))


def parse_tree(tok, pos=0, closing=False):
    """Nested dict of a dictionary: values are canonical strings, sub-dictionaries are dicts (last entry of a repeated key wins)."""
    tree = {}
    while pos < len(tok):
        t = tok[pos]
        if t == "}":
            if not closing:
                raise ValueError("unbalanced '}'")
            return tree, pos + 1
        if isinstance(t, tuple) or t in "{;()":
            raise ValueError(f"unexpected token {t!r}")
        key = t
        pos += 1
        if pos >= len(tok):
            raise ValueError(f"entry {key} has no value")
        if tok[pos] == "{":
            tree[key], pos = parse_tree(tok, pos + 1, True)
            continue
        vals, depth = [], 0
        while True:
            if pos >= len(tok):
                raise ValueError(f"entry {key} is not terminated by ';'")
            u = tok[pos]
            if u == ";" and depth == 0:
                pos += 1
                break
            depth += 1 if u == "(" else (-1 if u == ")" else 0)
            if depth < 0:
                raise ValueError("unbalanced ')'")
            if u == "{":
                sub, pos = parse_tree(tok, pos + 1, True)
                vals.append("{" + " ".join(f"{k} {v}" for k, v in sorted(sub.items(), key=lambda kv: kv[0]) if isinstance(v, str)) + "}")
                continue
            vals.append(("CODE#" + u[1]) if isinstance(u, tuple) else u)
            pos += 1
        tree[key] = " ".join(vals)
    if closing:
        raise ValueError("missing '}'")
    return tree, pos


def tree_of(path):
    return parse_tree(tokenize(_read(path)))[0]


def _get(d, *keys):
    """Value at a key path of a tree; raises KeyError with the path if absent."""
    x = d
    for k in keys:
        if not isinstance(x, dict) or k not in x:
            raise KeyError("/".join(keys))
        x = x[k]
    return x


def _first_dict(d):
    for v in d.values():
        if isinstance(v, dict):
            return v
    raise KeyError("no coefficient dictionary")


def _f(x):
    return float(x.split()[-1])


def _i(x):
    v = float(x.split()[-1])
    if v != int(v):
        raise ValueError(f"not an integer: {x}")
    return int(v)


def read_case(d):
    """Settings of one case from its dictionaries (parsed, so formatting does not matter).
    Raises (OSError, ValueError, KeyError, AttributeError) if a dictionary is missing, empty, unreadable or garbled."""
    j = lambda *p: os.path.join(d, *p)  # noqa: E731
    mech = _read(j("constant", "solid", "mechanicalProperties"))
    tr = tree_of(j("constant", "fluid", "transportProperties"))
    ct = tree_of(j("system", "controlDict"))
    ctl = _read(j("system", "controlDict"))
    fsi = tree_of(j("constant", "fsiProperties"))
    fvt = tree_of(j("system", "fluid", "fvSolution"))
    spr = tree_of(j("constant", "solid", "solidProperties"))
    dmd = tree_of(j("constant", "fluid", "dynamicMeshDict"))
    flp = tree_of(j("constant", "fluid", "fluidProperties"))
    mi_s, mi_f = mesh_info(j("system", "solid", "blockMeshDict")), mesh_info(j("system", "fluid", "blockMeshDict"))
    hs, hf = [(c[0], c[1], c[2]) for c in mi_s["counts"]], [(c[0], c[1], c[2]) for c in mi_f["counts"]]
    if len(hs) != 4 or len(hf) != 4:
        raise ValueError("blockMeshDict does not hold four blocks")
    mech_c = " ".join(f"{k} {v}" for k, v in parse_dict(mech))
    sol = _get(fvt, "solvers")
    pim = _get(fvt, "PIMPLE")
    fc = _first_dict(fsi)
    algo = _get(fsi, "fluidSolidInterface")
    sm = _get(spr, "solidModel")
    sc = spr.get(sm + "Coeffs") or _first_dict(spr)
    g = lambda pat, txt: re.search(pat, txt).group(1)  # noqa: E731
    pS, uS = _get(sol, '"p|pFinal"'), _get(sol, '"U|UFinal"')
    dcp = os.path.join(d, "system", "fluid", "decomposeParDict")
    s = dict(
        E=float(g(r"\bE E \[[^\]]*\] (\S+)", mech_c)), nu_s=float(g(r"\bnu nu \[[^\]]*\] (\S+)", mech_c)), rho_s=float(g(r"\brho rho \[[^\]]*\] (\S+)", mech_c)),
        nu_f=_f(_get(tr, "nu")), rho_f=_f(_get(tr, "rho")),
        deltaT=_f(_get(ct, "deltaT")), endTime=_f(_get(ct, "endTime")),
        nNonOrthogonalCorrectors=_i(_get(pim, "nNonOrthogonalCorrectors")), nOuterCorrectors=_i(_get(pim, "nOuterCorrectors")), nCorrectors=_i(_get(pim, "nCorrectors")),
        outerCorrTolerance=_f(_get(fc, "outerCorrTolerance")), outerCorrAbsoluteTolerance=_f(_get(fc, "outerCorrAbsoluteTolerance")),
        relaxationFactor=_f(_get(fc, "relaxationFactor")), couplingReuse=_i(_get(fc, "couplingReuse")), nOuterCorr=_i(_get(fc, "nOuterCorr")),
        allowUnconvergedCoupling=("allowUnconvergedCoupling" in fc or "allowUnconvergedCoupling" in fsi),
        requireAllResidualMeasures=("requireAllResidualMeasures" in fc or "requireAllResidualMeasures" in fsi),
        ntheta=sum(h[1] for h in hf), nrs=hs[0][0], nrf=hf[0][0], ntheta_solid=sum(h[1] for h in hs),
        p_solver=" ".join(f"{k} {v}" for k, v in sorted(pS.items())), U_solver=" ".join(f"{k} {v}" for k, v in sorted(uS.items())),
        fluidSolidInterface=algo, interfaceTransferMethod=_get(fc, "interfaceTransferMethod"),
        predictor=_get(fc, "predictor"), predictSolid=_get(fc, "predictSolid"), coupled=_get(fc, "coupled"),
        solidModel=sm, solutionAlgorithm=_get(sc, "solutionAlgorithm"), solid_predictor=_get(sc, "predictor"),
        stabilisation_type=_get(sc, "stabilisation", "momentum", "type"), stabilisation_scaleFactor=_f(_get(sc, "stabilisation", "momentum", "scaleFactor")),
        law=g(r"type (\w+Elastic\w*)", mech_c), planeStress=g(r"planeStress (\w+)", mech_c),
        momentumPredictor=_get(pim, "momentumPredictor"),
        fluxConsistentPatches=" ".join(_get(pim, "fluxConsistentPatches").replace("(", " ").replace(")", " ").split()),
        fluidModel=_get(flp, "fluidModel"), simulationType=_get(tree_of(j("constant", "fluid", "momentumTransport")), "simulationType"),
        dynamicFvMesh=_get(dmd, "dynamicFvMesh"), motionSolver=_get(dmd, "motionSolver"), diffusivity=_get(dmd, "diffusivity"),
        diagnostics=("current" if ("reduce(" in ctl and "hi[4]" in ctl) else ("reduce, no n=3,4 columns" if "reduce(" in ctl else "legacy serial-only")),
        decompose=_i(_get(tree_of(dcp), "numberOfSubdomains")) if os.path.isfile(dcp) else None,
    )
    for k, v in s.items():
        if v is None and k != "decompose":
            raise ValueError(f"missing dictionary entry {k}")
    # semantic fingerprints, patch types and coded blocks
    bc, code, fp = {}, {}, {}
    for rel in DICTS:
        txt = _read(j(*rel.split("/")))
        entries = parse_dict(txt, mask=("outerCorrAbsoluteTolerance",) if rel == "constant/fsiProperties" else ())
        fp[rel] = sha(canon(entries, code_mask=True))
        hs_ = code_hashes_of(entries)
        if hs_:
            code[rel] = hs_
        if rel.startswith("0/"):
            bc[rel] = patch_types(txt)
    ce = parse_dict(ctl, mask=())
    cd = dict(ce)
    s["controlDict_entries"] = {k: v for k, v in ce if k not in ("functions", "endTime", "writeInterval")}
    s["controlDict_fp"] = sha(canon(ce, mask=("endTime", "writeInterval"), code_mask=True))
    s["diag_code_hashes"] = sorted(code_hashes_of(ce))
    s["mesh_solid"], s["mesh_fluid"] = {k: v for k, v in mi_s.items() if k != "vertices"}, {k: v for k, v in mi_f.items() if k != "vertices"}
    pcode = _read(j("0", "fluid", "p"))
    m = re.search(r"const scalar seed = ([^;]+);", pcode)
    s["p_seed"] = float(m.group(1)) if m else None
    s["fingerprints"], s["code_hashes"], s["bc_types"] = fp, code, bc
    s["angle_table"] = _tab(j("constant", "timeVsAngle"))
    s["pc_table"] = _tab(j("constant", "pcTable"))
    return s


def expected_tables(cfg, thr, times):
    D = thr["family_consistency"]["drives"]["value"]
    th, P0 = D[cfg]["theta0"], D[cfg]["P0"]
    p = rt.Params(theta0=th, Omega=D["Omega"], P0=P0, T_ramp=D["T_ramp"], ramp="Cinf")
    ang = np.array([math.degrees(rt.rotation_drive(p, t)[0]) for t in times])
    pc = np.array([rt.pressure_drive(p, t) / p.rho_f for t in times])
    return ang, pc


def worst_status(findings):
    return ST.worst(f["status"] for f in findings)


def check_family(work, cfg, thr, runs=None):
    """Returns dict(status, findings=[{check, status, detail}], settings={mesh: summary}). Never raises."""
    try:
        return _check_family(work, cfg, thr)
    except Exception as e:  # a garbled input must give INCOMPLETE, never a traceback
        return dict(status="INCOMPLETE", findings=[dict(check="family_consistency", status="INCOMPLETE",
                                                        detail=f"{type(e).__name__}: {e}")], settings={})


def _check_family(work, cfg, thr):
    F = thr["family_consistency"]
    findings = []

    def add(name, sev, detail):
        findings.append(dict(check=name, status=sev, detail=detail))

    S = {}
    for m in MESHES:
        d = os.path.join(work, f"{cfg}_{m}")
        if not os.path.isdir(d):
            add(f"run_{m}", "INCOMPLETE", f"case {cfg}_{m} missing")
            continue
        try:
            S[m] = read_case(d)
        except Exception as e:      # missing, empty, unreadable or garbled dictionary
            add(f"dictionaries_{m}", "INCOMPLETE", f"unreadable or garbled dictionaries of {cfg}_{m}: {type(e).__name__}: {e}")
    if len(S) < 3:
        return dict(status="INCOMPLETE", findings=findings, settings={})
    ph = F["physics"]["value"]
    cp = F["coupling"]["value"]
    for k, want in ph.items():
        bad = {m: S[m][k] for m in MESHES if abs(S[m][k] - want) > 1e-12 * max(1.0, abs(want))}
        add(f"physics_{k}", "INVALID-INPUT" if bad else "PASS", f"{k} expected {want:g}" + (f"; differs: {bad}" if bad else "; all meshes equal"))
    for k, want in cp.items():
        bad = {m: S[m][k] for m in MESHES if S[m][k] != want}
        add(f"coupling_{k}", "INVALID-INPUT" if bad else "PASS", f"{k} expected {want:g}" + (f"; differs: {bad}" if bad else "; all meshes equal"))
    for k in ("allowUnconvergedCoupling", "requireAllResidualMeasures"):
        bad = {m: S[m][k] for m in MESHES if S[m][k]}
        add(f"coupling_{k}", "INVALID-INPUT" if bad else "PASS", f"{k} must be unset; set in: {sorted(bad) if bad else 'none'}")
    doc = F["documented_absTol_by_mesh"]["value"]
    at = {m: S[m]["outerCorrAbsoluteTolerance"] for m in MESHES}
    if all(at[m] == doc[m] for m in MESHES):
        add("coupling_outerCorrAbsoluteTolerance", "CONDITIONAL",
            f"{at}: the documented per-mesh pattern (1e-12 m on M1 and M2, 1e-10 m on M3); its effect is bounded on M2 only, not on M3")
    elif len(set(at.values())) == 1:
        add("coupling_outerCorrAbsoluteTolerance", "PASS", f"equal on all meshes: {at}")
    else:
        add("coupling_outerCorrAbsoluteTolerance", "INVALID-INPUT", f"undocumented absolute tolerances {at}")
    # named settings: equal on all meshes AND equal to the documented benchmark set
    doc_set = F["documented_benchmark_settings"]["value"]
    for k, want in doc_set.items():
        vals = {m: S[m][k] for m in MESHES}
        if len(set(map(str, vals.values()))) > 1:
            add(f"setting_{k}", "INVALID-INPUT", f"differs between the meshes: {vals}")
        elif str(vals["M1"]) != str(want):
            add(f"setting_{k}", "INVALID-INPUT", f"{vals['M1']} on all meshes, documented benchmark value {want} (not the benchmark configuration; the bands do not apply)")
        else:
            add(f"setting_{k}", "PASS", f"{want} on all meshes")
    for k in ("p_solver", "U_solver"):
        vals = {m: S[m][k] for m in MESHES}
        ok = len(set(vals.values())) == 1 and vals["M1"] == F["documented_solver_text"]["value"][k]
        add(f"numerics_{k}", "PASS" if ok else "INVALID-INPUT", "identical on all meshes and equal to the documented set" if ok else f"differs: {vals}")
    # boundary conditions: patch types, coded blocks, seed
    bct = {m: S[m]["bc_types"] for m in MESHES}
    exp_bc = F["documented_bc_types"]["value"]
    same = all(bct[m] == bct["M1"] for m in MESHES)
    ok = same and {k: [list(x) for x in v] for k, v in bct["M1"].items()} == exp_bc
    add("boundary_conditions_types", "PASS" if ok else "INVALID-INPUT",
        "patch types of U, p, pointDisplacement and D equal on all meshes and documented" if ok else f"patch types differ from the meshes or the documented set: {bct if not same else bct['M1']}")
    cc = F["documented_code_hashes"]["value"]
    for rel in ("0/fluid/U", "0/fluid/p"):
        hs_ = {m: S[m]["code_hashes"].get(rel) for m in MESHES}
        allowed = cc[rel]
        if all(h == allowed["current"] for h in hs_.values()):
            add(f"coded_condition_{rel}", "PASS", "normalised code identical on all meshes and documented")
        elif rel == "0/fluid/p" and all(h in (allowed["current"], allowed.get("legacy")) for h in hs_.values()) \
                and all(S[m]["p_seed"] == 0.0 for m in MESHES):
            add(f"coded_condition_{rel}", "CONDITIONAL", "older seeded-pressure code on some meshes with seed 0: identical boundary values (documented equivalent)")
        else:
            add(f"coded_condition_{rel}", "INVALID-INPUT", f"coded condition differs from the documented set or between the meshes (seed {[S[m]['p_seed'] for m in MESHES]})")
    # fingerprint safety net: every small dictionary
    docfp = F["documented_fingerprints"]["value"]
    skip = {"0/fluid/p", "0/fluid/U"}   # judged above through patch types and code hashes
    bad_between, bad_doc = [], []
    for rel in DICTS:
        if rel in skip:
            continue
        v = {m: S[m]["fingerprints"][rel] for m in MESHES}
        if len(set(v.values())) > 1:
            bad_between.append(rel)
        elif v["M1"] != docfp.get(rel):
            bad_doc.append(rel)
    # the p and U files without their code: compare text with code masked
    for rel in sorted(skip):
        v = {m: S[m]["fingerprints"][rel] for m in MESHES}
        if len(set(v.values())) > 1:
            bad_between.append(rel)
        elif v["M1"] != docfp.get(rel):
            bad_doc.append(rel)
    add("dictionary_fingerprints", "INVALID-INPUT" if (bad_between or bad_doc) else "PASS",
        "normalised text of all small dictionaries equal on all meshes and to the documented benchmark set" if not (bad_between or bad_doc)
        else f"differ between the meshes: {bad_between}; differ from the documented set: {bad_doc}")
    # controlDict: every top-level entry (before and after `functions`) except endTime and writeInterval; function-object parameters
    docc = F["documented_controlDict"]["value"]
    bad_ce = {m: {k: (S[m]["controlDict_entries"].get(k), docc["entries"].get(k)) for k in set(S[m]["controlDict_entries"]) | set(docc["entries"])
                  if S[m]["controlDict_entries"].get(k) != docc["entries"].get(k)} for m in MESHES}
    bad_ce = {m: v for m, v in bad_ce.items() if v}
    add("controlDict_entries", "INVALID-INPUT" if bad_ce else "PASS",
        "all top-level controlDict entries (before and after the functions block) equal the documented set, except endTime and writeInterval" if not bad_ce
        else f"controlDict entries differ from the documented set (found, documented): {bad_ce}")
    bad_cfp = [m for m in MESHES if S[m]["controlDict_fp"] != docc["fingerprint"]]
    add("controlDict_semantic", "INVALID-INPUT" if bad_cfp else "PASS",
        "canonical controlDict (including the function-object parameters, code masked) equals the documented set" if not bad_cfp
        else f"canonical controlDict differs from the documented set on {bad_cfp}")
    dh = {m: S[m]["diag_code_hashes"] for m in MESHES}
    known = (docc["diag_code_current"], docc["diag_code_legacy"])
    if any(h not in known for h in dh.values()):
        add("diagnostics_function_objects", "INVALID-INPUT", f"function-object code of {[m for m in MESHES if dh[m] not in known]} is neither the documented current nor the documented older version")
    else:
        lab = {m: ("current" if dh[m] == docc["diag_code_current"] else "legacy serial-only") for m in MESHES}
        add("diagnostics_function_objects", "PASS" if len(set(lab.values())) == 1 else "CONDITIONAL",
            f"{lab}" + ("" if len(set(lab.values())) == 1 else ": different versions; the shared columns used here are the same quantities"))
    # mesh geometry and grading, both regions
    G = F["documented_geometry"]["value"]
    for reg, key, want in (("solid", "mesh_solid", [G["A"], G["B"]]), ("fluid", "mesh_fluid", [G["B"], G["c"]])):
        bad = {m: S[m][key]["radii"] for m in MESHES if S[m][key]["radii"] != [round(x, 9) for x in want]}
        badz = {m: S[m][key]["z"] for m in MESHES if S[m][key]["z"] != [0.0, G["depth"]]}
        bada = {m: S[m][key]["arc_radii"] for m in MESHES if S[m][key]["arc_radii"] != [round(x, 9) for x in want]}
        add(f"geometry_{reg}", "INVALID-INPUT" if (bad or badz or bada) else "PASS",
            f"{reg} region radii {want} (inner/outer), depth {G['depth']}: vertices and arc points " +
            ("as documented on all meshes" if not (bad or badz or bada) else f"differ: radii {bad}, z {badz}, arc radii {bada}"))
        gr = {m: set(S[m][key]["grading"]) for m in MESHES}
        add(f"grading_{reg}", "PASS" if all(g == {("1", "1", "1")} for g in gr.values()) else "INVALID-INPUT",
            "simpleGrading (1 1 1) in every block on all meshes" if all(g == {("1", "1", "1")} for g in gr.values()) else f"grading differs: {gr}")
        fps = {m: (S[m][key]["geometry_fp"], S[m][key]["boundary_fp"]) for m in MESHES}
        okfp = len(set(fps.values())) == 1 and fps["M1"] == tuple(F["documented_mesh_fingerprints"]["value"][reg])
        add(f"mesh_geometry_{reg}", "PASS" if okfp else "INVALID-INPUT",
            "vertices, arcs, block topology and patch definitions equal on all meshes and to the documented set" if okfp else f"geometry or patch definitions differ: {fps}")
    rk = {m: S[m]["decompose"] for m in MESHES}
    add("ranks_decomposition", "PASS" if len(set(rk.values())) == 1 else "CONDITIONAL",
        f"stored numberOfSubdomains {rk}; scotch partitions are not reproducible (noise allowance, numerical_noise section)")
    # mesh family and refinement ratio
    exp = F["expected_mesh_ntheta_nrs_nrf"]["value"]
    got = {m: [S[m]["ntheta"], S[m]["nrs"], S[m]["nrf"]] for m in MESHES}
    ok = all(got[m] == exp[m] for m in MESHES) and all(S[m]["ntheta_solid"] == S[m]["ntheta"] for m in MESHES)
    rr = F["expected_refinement_ratio"]["value"]
    ratios = [[got["M2"][i] / got["M1"][i], got["M3"][i] / got["M2"][i]] for i in range(3)]
    ok = ok and all(x == rr for r in ratios for x in r)
    add("mesh_family", "PASS" if ok else "INVALID-INPUT", f"ntheta x nrs x nrf {got}; refinement ratios per direction {ratios} (expected {rr})")
    # drives
    tol = F["table_relative_tolerance"]["value"]
    for m in MESHES:
        a_t, p_t = S[m]["angle_table"], S[m]["pc_table"]
        if a_t.shape != p_t.shape or not np.allclose(a_t[:, 0], p_t[:, 0]):
            add(f"drives_{m}", "INVALID-INPUT", f"timeVsAngle and pcTable differ in length or time grid ({a_t.shape} against {p_t.shape}): truncated or inconsistent series")
            continue
        t = a_t[:, 0]
        ang, pc = expected_tables(cfg, thr, t)
        e1 = float(np.max(np.abs(a_t[:, 1] - ang)) / max(1e-300, np.max(np.abs(ang)))) if np.max(np.abs(ang)) > 0 else float(np.max(np.abs(a_t[:, 1])))
        e2 = float(np.max(np.abs(p_t[:, 1] - pc)) / max(1e-300, np.max(np.abs(pc)))) if np.max(np.abs(pc)) > 0 else float(np.max(np.abs(p_t[:, 1])))
        add(f"drives_{m}", "PASS" if (e1 <= tol and e2 <= tol) else "INVALID-INPUT",
            f"timeVsAngle and pcTable against rt_reference {cfg} (theta0, P0): relative difference {e1:.1e}, {e2:.1e}")
    tg = {m: (len(S[m]["angle_table"]), float(S[m]["angle_table"][-1, 0])) for m in MESHES}
    ng = int(round(ph["endTime"] / ph["deltaT"])) + 3
    add("time_grid_tables", "PASS" if (len(set(tg.values())) == 1 and tg["M1"][0] == ng) else "INVALID-INPUT", f"{tg} (expected {ng} rows)")
    summary = {m: {k: v for k, v in S[m].items() if k not in ("angle_table", "pc_table", "fingerprints", "code_hashes", "bc_types", "mesh_solid", "mesh_fluid", "controlDict_entries")} for m in MESHES}
    return dict(status=worst_status(findings), findings=findings, settings=summary)
