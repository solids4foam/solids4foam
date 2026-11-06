#!/usr/bin/env python3
"""Unit tests of the oscillatingCylinderInChannel verification driver, with
synthetic data, and consistency tests of the tracked reference files,
compact results, results.md and figures (no OpenFOAM, no SLURM: sbatch and
squeue are never called).

    python3 -m unittest discover -s scripts -p 'test_*.py'
"""

import csv
import importlib.util
import json
import math
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock

SCRIPT = Path(__file__).resolve().parent / "oscillating_cylinder_verification.py"
spec = importlib.util.spec_from_file_location("ocv", SCRIPT)
ocv = importlib.util.module_from_spec(spec)
spec.loader.exec_module(ocv)

try:
    import scipy.signal  # noqa: F401
    HAVE_SCIPY = True
except ImportError:
    HAVE_SCIPY = False

HEADER_CUTLINK = ("# Time\tforce (x y z)\ttorque (x y z)\tinertia of the "
                  "fluid inside (x y z)\tforce from the momentum exchange "
                  "(x y z)")
HEADER_FORCING = ("# Time\tforce (x y z)\ttorque (x y z)\tinertia of the "
                  "fluid inside (x y z)")


def force_line(t, *vectors):
    return f"{t}\t" + "\t".join(" ".join(str(c) for c in v)
                                for v in vectors)


class Sandbox(unittest.TestCase):
    """Runs each test with WORK_DIR, CASES_DIR, MANIFEST and OUTPUT_DIR in a
    temporary directory"""

    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp())
        self.patches = [
            mock.patch.object(ocv, "WORK_DIR", self.tmp / "work"),
            mock.patch.object(ocv, "CASES_DIR", self.tmp / "work/cases"),
            mock.patch.object(ocv, "MANIFEST", self.tmp / "work/manifest.json"),
            mock.patch.object(ocv, "OUTPUT_DIR", self.tmp / "postProcessing"),
        ]
        for p in self.patches:
            p.start()

    def tearDown(self):
        for p in self.patches:
            p.stop()
        shutil.rmtree(self.tmp)

    def write_force_file(self, case, start, header, lines):
        d = case / "postProcessing/immersedBoundary" / start
        d.mkdir(parents=True, exist_ok=True)
        (d / "cylinder.dat").write_text(
            "# Immersed body cylinder\n" + header + "\n" + "\n".join(lines)
            + "\n")


class ForceColumns(Sandbox):
    def test_header_cutlink(self):
        cols = ocv.parse_force_header(HEADER_CUTLINK)
        self.assertEqual(cols, {"force": 1, "torque": 4, "inertia": 7,
                                "momentumExchange": 10})

    def test_header_both_alternatives(self):
        cols = ocv.parse_force_header(
            HEADER_CUTLINK + "\tforce from the surface traction (x y z)")
        self.assertEqual(cols["surfaceTraction"], 13)

    def test_unknown_column_rejected(self):
        with self.assertRaises(RuntimeError):
            ocv.parse_force_header(HEADER_FORCING + "\tsomething (x y z)")

    def test_cutlink_estimates_never_added(self):
        cols = ocv.parse_force_header(HEADER_CUTLINK)
        rows = [[0.1, 1.0, 2.0, 0, 0, 0, 0, 0, 0, 0, 10.0, 20.0, 0]]
        f = ocv.hydrodynamic_forces("cutLink", cols, rows)
        self.assertEqual(f["surfaceTraction"], [(1.0, 2.0)])
        self.assertEqual(f["momentumExchange"], [(10.0, 20.0)])
        self.assertEqual(f["total"], [(1.0, 2.0)])  # not 11, 22

    def test_cutlink_nonzero_inertia_rejected(self):
        cols = ocv.parse_force_header(HEADER_CUTLINK)
        rows = [[0.1, 1.0, 2.0, 0, 0, 0, 0, 0.5, 0, 0, 10.0, 20.0, 0]]
        with self.assertRaises(RuntimeError):
            ocv.hydrodynamic_forces("cutLink", cols, rows)

    def test_forcing_plus_inertia(self):
        cols = ocv.parse_force_header(HEADER_FORCING)
        rows = [[0.1, 1.0, 2.0, 0.0, 9, 9, 9, 0.25, -0.5, 0.0]]
        for method in ("penalty", "incremental"):
            f = ocv.hydrodynamic_forces(method, cols, rows)
            self.assertEqual(f["forcing"], [(1.0, 2.0)])
            self.assertEqual(f["inertia"], [(0.25, -0.5)])
            self.assertEqual(f["total"], [(1.25, 1.5)])

    def test_truncated_line_dropped(self):
        case = self.tmp / "c"
        self.write_force_file(case, "0", HEADER_FORCING, [
            force_line(0.1, (1, 2, 0), (0, 0, 0), (0, 0, 0)),
            "0.2\t1 2 0\t0 0"])
        cols, rows = ocv.read_force_file(
            case / "postProcessing/immersedBoundary/0/cylinder.dat")
        self.assertEqual(len(rows), 1)

    def test_restart_segments(self):
        case = self.tmp / "c"
        self.write_force_file(case, "0", HEADER_FORCING, [
            force_line(t, (t, 0, 0), (0, 0, 0), (0, 0, 0))
            for t in (0.1, 0.2, 0.3, 0.4)])
        self.write_force_file(case, "0.3", HEADER_FORCING, [
            force_line(t, (10*t, 0, 0), (0, 0, 0), (0, 0, 0))
            for t in (0.3, 0.4, 0.5)])
        rows = ocv.force_rows(case)
        self.assertEqual([r[0] for r in rows], [0.1, 0.2, 0.3, 0.4, 0.5])
        self.assertEqual([r[1] for r in rows], [0.1, 0.2, 3.0, 4.0, 5.0])

    def test_restart_with_changed_columns_rejected(self):
        with self.assertRaises(RuntimeError):
            ocv.merge_segments([(0.0, {"force": 1}, []),
                                (1.0, {"force": 1, "inertia": 4}, [])])


class Scaling(unittest.TestCase):
    def test_parameters_from_tutorial(self):
        p = ocv.parameters(ocv.TUTORIAL_DIR)
        self.assertAlmostEqual(p["U_ref"], 2*math.pi*0.25/4, places=12)
        self.assertAlmostEqual(p["D"], 0.1, places=6)
        self.assertAlmostEqual(p["Lz"], 0.1, places=12)
        self.assertAlmostEqual(p["rho"], 1.0)
        self.assertAlmostEqual(p["nu"], 1e-3)
        expected = 2/(1.0*(2*math.pi*0.25/4)**2*0.1*0.1)
        self.assertAlmostEqual(p["coeff_scale"]/expected, 1.0, places=6)

    def test_body_centre(self):
        p = ocv.parameters(ocv.TUTORIAL_DIR)
        self.assertAlmostEqual(ocv.body_centre(p, 0.5),
                               1.1 + 0.25*math.sin(math.pi/4), places=6)
        self.assertAlmostEqual(ocv.body_centre(p, 0.5), 1.27678, places=4)


class Series(unittest.TestCase):
    def test_clean_series(self):
        t, y = ocv.clean_series([0.3, 0.1, 0.2, 0.2], [3, 1, 2, 4])
        self.assertEqual(t, [0.1, 0.2, 0.3])
        self.assertEqual(y, [1, 3, 3])

    def test_interpolation_and_range(self):
        t, y = [0.0, 1.0, 2.0], [0.0, 10.0, 0.0]
        self.assertAlmostEqual(ocv.interpolate(t, y, 0.25), 2.5)
        self.assertTrue(math.isnan(ocv.interpolate(t, y, 2.5)))
        self.assertTrue(math.isnan(ocv.interpolate(t, y, -0.1)))

    def test_rms_constant_offset(self):
        t = [i*0.01 for i in range(801)]
        ref = [math.sin(2*math.pi*x/4) for x in t]
        sim = [r + 0.1 for r in ref]
        e = ocv.rms_error(t, sim, t, ref)
        self.assertAlmostEqual(e["rms"], 0.1, places=10)
        self.assertAlmostEqual(e["t0"], 0.25)
        self.assertAlmostEqual(e["t1"], 7.5)
        # Mean of sin^2(pi t/2) over [0.25, 7.5], which is not a whole
        # number of periods
        w = math.pi/2
        mean = ((7.5/2 - math.sin(2*w*7.5)/(4*w))
                - (0.25/2 - math.sin(2*w*0.25)/(4*w)))/7.25
        self.assertAlmostEqual(e["ref_rms"], math.sqrt(mean), places=4)

    def test_rms_irregular_reference_time_weighted(self):
        # A reference sampled densely in one half must not dominate
        ref_t = [0.25 + 0.001*i for i in range(1000)] + [1.5, 3.0, 4.5, 6.0,
                                                         7.5]
        ref_y = [0.0]*len(ref_t)
        sim_t = [i*0.01 for i in range(801)]
        sim_y = [1.0 if x > 1.25 else 0.0 for x in sim_t]
        e = ocv.rms_error(sim_t, sim_y, ref_t, ref_y)
        # The step is between the simulation samples at 1.25 and 1.26
        self.assertAlmostEqual(e["rms"]**2, (7.5 - 1.255)/(7.5 - 0.25),
                               places=3)

    def test_rms_no_extrapolation(self):
        t = [i*0.01 for i in range(801)]
        ref_t = [x for x in t if x <= 7.0]
        e = ocv.rms_error(t, [1.0]*len(t), ref_t, [1.0]*len(ref_t))
        self.assertAlmostEqual(e["t1"], 7.0)
        sim_t = [x for x in t if x <= 5.0]
        e = ocv.rms_error(sim_t, [1.0]*len(sim_t), t, [1.0]*len(t))
        self.assertAlmostEqual(e["t1"], 5.0)

    def test_lag_diagnostic(self):
        t = [i*0.0025 for i in range(3201)]
        ref = [math.sin(2*math.pi*x/4) for x in t]
        sim = [math.sin(2*math.pi*(x - 0.02)/4) for x in t]
        self.assertAlmostEqual(ocv.best_lag(t, sim, t, ref), 0.02,
                               places=3)
        # The primary error is not shifted
        self.assertGreater(ocv.rms_error(t, sim, t, ref)["rms"], 0.02)

    def test_extrema(self):
        t = [i*0.01 for i in range(801)]
        y = [math.sin(2*math.pi*x/4) for x in t]
        tmax, ymax, tmin, ymin = ocv.extrema(t, y, 0.25, 4.0)
        self.assertAlmostEqual(tmax, 1.0, places=3)
        self.assertAlmostEqual(ymax, 1.0, places=4)
        self.assertAlmostEqual(tmin, 3.0, places=3)

    @unittest.skipUnless(HAVE_SCIPY, "scipy is not installed")
    def test_high_frequency_rms(self):
        t = [i*0.0025 for i in range(3201)]
        y = [3*math.sin(2*math.pi*x/4) + 0.01*math.sin(2*math.pi*40*x)
             for x in t]
        self.assertAlmostEqual(ocv.high_frequency_rms(t, y, 5.0),
                               0.01/math.sqrt(2), places=3)
        smooth = [3*math.sin(2*math.pi*x/4) for x in t]
        self.assertLess(ocv.high_frequency_rms(t, smooth, 5.0), 1e-4)
        # A large start-up spike outside the window does not leak into it
        spiky = [30.0 if x < 0.02 else v for x, v in zip(t, smooth)]
        self.assertLess(ocv.high_frequency_rms(t, spiky, 5.0), 1e-4)

    @unittest.skipUnless(HAVE_SCIPY, "scipy is not installed")
    def test_low_pass_diagnostic(self):
        t = [i*0.0025 for i in range(3201)]
        smooth = [3*math.sin(2*math.pi*x/4) for x in t]
        noisy = [v + 0.5*math.sin(2*math.pi*40*x) for x, v in zip(t, smooth)]
        tl, yl = ocv.low_pass(t, noisy, 5.0)
        self.assertGreaterEqual(tl[0], 0.25 - 1e-12)
        e = ocv.rms_error(tl, yl, t, smooth, window=(0.5, 7.25))
        self.assertLess(e["rms"], 1e-3)
        # The primary error keeps the noise
        self.assertAlmostEqual(ocv.rms_error(t, noisy, t, smooth)["rms"],
                               0.5/math.sqrt(2), places=2)

    @unittest.skipUnless(HAVE_SCIPY, "scipy is not installed")
    def test_high_frequency_rejects_irregular(self):
        with self.assertRaises(RuntimeError):
            ocv.high_frequency_rms([0.25 + 0.01*i + (0.003 if i % 2 else 0)
                                    for i in range(100)], [0.0]*100, 5.0)


class Generation(Sandbox):
    class A:
        runs = ["smoke"]
        ranks = 0
        force = False
        dry_run = False
        traction_quadrature = "greedy"

    def test_generate_settings(self):
        ocv.cmd_generate(self.A)
        case = ocv.CASES_DIR / "cutLink-L1-dt0p0025-t0p5"
        control = (case / "system/controlDict").read_text()
        self.assertRegex(control, r"(?m)^deltaT\s+0\.0025;")
        self.assertRegex(control, r"(?m)^adjustTimeStep\s+no;")
        self.assertRegex(control, r"(?m)^endTime\s+0\.5;")
        self.assertRegex(control, r"(?m)^timePrecision\s+10;")
        self.assertRegex((case / "system/blockMeshDict").read_text(),
                         r"(?m)^MESH_LEVEL\s+1;")
        fv = (case / "constant/fvOptions").read_text()
        self.assertRegex(fv, r"method\s+cutLink;")
        self.assertRegex(fv, r"forceEstimator\s+surfaceTraction;")
        # The historical quadrature is written explicitly
        self.assertRegex(fv, r"tractionQuadrature\s+greedy;")
        self.assertRegex(fv, r"amplitude\s+0\.25;")
        self.assertRegex(fv, r"period\s+4;")
        sol = (case / "system/fvSolution").read_text()
        self.assertRegex(sol, r"nOuterCorrectors\s+1;")
        self.assertRegex(sol, r"(?m)^\s+nCorrectors\s+3;")
        self.assertIn("OMP_NUM_THREADS=1", (case / "Allrun.verify").read_text())
        self.assertIn("./Allrun", (case / "Allrun.verify").read_text())
        manifest = json.loads(ocv.MANIFEST.read_text())
        self.assertIn("cutLink-L1-dt0p0025-t0p5", manifest["cases"])
        # The generated case gives the same physical parameters
        self.assertAlmostEqual(ocv.parameters(case)["coeff_scale"],
                               ocv.parameters(ocv.TUTORIAL_DIR)["coeff_scale"])

    def test_generate_methods(self):
        class B(self.A):
            runs = ["R1B", "R1C", "R2C"]
            ranks = 4
        ocv.cmd_generate(B)
        pen = ocv.CASES_DIR / "penalty-L1-dt0p0025"
        inc = ocv.CASES_DIR / "incremental-L2-dt0p00125"
        self.assertRegex((pen / "constant/fvOptions").read_text(),
                         r"occupancy\s+signedDistance;")
        fv = (inc / "constant/fvOptions").read_text()
        self.assertRegex(fv, r"couplingCoeff\s+0\.2;")
        self.assertRegex(fv, r"occupancy\s+vertexFraction;")
        self.assertRegex((inc / "system/fvSolution").read_text(),
                         r"nOuterCorrectors\s+5;")
        self.assertRegex((inc / "system/decomposeParDict").read_text(),
                         r"numberOfSubdomains\s+4;")
        self.assertIn("solids4Foam::runParallel",
                      (inc / "Allrun.verify").read_text())

    def test_generate_ghost_cell_and_level_4(self):
        class G(self.A):
            runs = ["R1D", "R5A"]
        ocv.cmd_generate(G)
        ghost = ocv.CASES_DIR / "ghostCell-L1-dt0p0025"
        fv = (ghost / "constant/fvOptions").read_text()
        self.assertRegex(fv, r"method\s+ghostCell;")
        self.assertRegex(fv, r"imageDistance\s+1\.5;")
        self.assertRegex(fv, r"imageInterpolation\s+cellPoint;")
        self.assertNotIn("tractionQuadrature", fv)
        self.assertNotIn("forceEstimator", fv)
        self.assertRegex((ghost / "system/fvSolution").read_text(),
                         r"nOuterCorrectors\s+1;")
        l4 = ocv.CASES_DIR / "cutLink-L4-dt0p0003125"
        self.assertEqual(json.loads((l4 / "case.json").read_text())["ranks"],
                         8)
        self.assertRegex((l4 / "system/blockMeshDict").read_text(),
                         r"(?m)^MESH_LEVEL\s+4;")

    def test_outer_corrector_runs(self):
        class N(self.A):
            runs = ["N5", "N10", "N20"]
        ocv.cmd_generate(N)
        for n, name in ((5, "incremental-L2-dt0p00125-t2p25-np2"),
                        (10, "incremental-L2-dt0p00125-t2p25-np2-nOuter10"),
                        (20, "incremental-L2-dt0p00125-t2p25-np2-nOuter20")):
            case = ocv.CASES_DIR / name
            self.assertRegex((case / "system/fvSolution").read_text(),
                             rf"nOuterCorrectors\s+{n};")
            meta = json.loads((case / "case.json").read_text())
            self.assertEqual((meta["nOuterCorrectors"], meta["end_time"],
                              meta["ranks"]), (n, 2.25, 2))

    def test_quadrature_is_explicit_not_read_from_source(self):
        # No solver source is searched: the default option is greedy and
        # the entry is written whatever the sources say
        self.assertFalse(hasattr(ocv, "SOLVER_SOURCE"))
        self.assertEqual(ocv.HISTORICAL_QUADRATURE, "greedy")
        self.assertEqual(ocv.quadrature_tag("cutLink", "greedy"), "")
        self.assertEqual(ocv.quadrature_tag("cutLink", "clustered"),
                         "clustered")
        self.assertEqual(ocv.quadrature_tag("penalty", "clustered"), "")
        self.assertIn("greedy", ocv.QUADRATURES)
        ocv.cmd_generate(self.A)
        meta = json.loads((ocv.CASES_DIR / "cutLink-L1-dt0p0025-t0p5"
                           / "case.json").read_text())
        self.assertEqual(meta["traction_quadrature"], "greedy")
        self.assertEqual(meta["fvOptions"]["tractionQuadrature"], "greedy")

    def test_default_option_is_greedy(self):
        # The default of the command line is the historical quadrature, not
        # that of the solver (clustered after solids4foam pull request #612)
        with mock.patch.object(ocv, "cmd_dry_run", return_value=0) as m:
            ocv.main(["dry-run"])
        self.assertEqual(m.call_args[0][0].traction_quadrature, "greedy")

    def test_other_quadrature_tags_case(self):
        class Q(self.A):
            traction_quadrature = "clustered"
        ocv.cmd_generate(Q)
        case = ocv.CASES_DIR / "cutLink-L1-dt0p0025-t0p5-clustered"
        self.assertRegex((case / "constant/fvOptions").read_text(),
                         r"tractionQuadrature\s+clustered;")

    def test_other_methods_have_no_quadrature(self):
        class P(self.A):
            runs = ["R1B"]
        ocv.cmd_generate(P)
        fv = (ocv.CASES_DIR / "penalty-L1-dt0p0025"
              / "constant/fvOptions").read_text()
        self.assertNotIn("tractionQuadrature", fv)

    def test_existing_case_protected(self):
        ocv.cmd_generate(self.A)
        with self.assertRaises(RuntimeError):
            ocv.cmd_generate(self.A)

    def test_submitted_case_protected(self):
        ocv.cmd_generate(self.A)
        m = json.loads(ocv.MANIFEST.read_text())
        m["cases"]["cutLink-L1-dt0p0025-t0p5"]["job_id"] = "123"
        ocv.save_manifest(m)

        class F(self.A):
            force = False
        with self.assertRaises(RuntimeError):
            ocv.cmd_generate(F)

    def test_run_with_fixed_ranks(self):
        class P(self.A):
            runs = ["O4"]
            ranks = 8
        ocv.cmd_generate(P)
        case = ocv.CASES_DIR / "cutLink-L2-dt0p00125-np1"
        self.assertEqual(json.loads((case / "case.json").read_text())["ranks"],
                         1)
        self.assertIn("./Allrun", (case / "Allrun.verify").read_text())

    def test_unknown_run(self):
        class U(self.A):
            runs = ["R9"]
        with self.assertRaises(RuntimeError):
            ocv.cmd_generate(U)


class Submission(Sandbox):
    class S:
        runs = ["R1A"]
        ranks = 0
        force = False
        env_script = "/private/env.sh"
        partition = "main"
        max_cpus = 96
        time = ""
        mem = ""
        sbatch_args = "--qos=normal"
        yes = True

    def setUp(self):
        super().setUp()

        class G:
            runs = ["R1A"]
            ranks = 0
            force = False
        ocv.cmd_generate(G)

    def test_job_script(self):
        calls = []

        def fake_run(cmd, **kw):
            calls.append(cmd)
            return mock.Mock(stdout="4242\n")
        with mock.patch.object(ocv, "queue_jobs", return_value={}), \
                mock.patch.object(ocv.shutil, "which", return_value="/x"), \
                mock.patch.object(ocv.subprocess, "run", side_effect=fake_run):
            ocv.cmd_submit(self.S)
        cmd = calls[0]
        self.assertEqual(cmd[0], "sbatch")
        for opt in ("--ntasks=1", "--cpus-per-task=1", "--partition=main",
                    "--qos=normal", "--mem=2G", "--time=06:00:00"):
            self.assertIn(opt, cmd)
        self.assertIn("source /private/env.sh", cmd[-1])
        self.assertIn("./Allrun.verify", cmd[-1])
        m = json.loads(ocv.MANIFEST.read_text())
        self.assertEqual(m["cases"]["cutLink-L1-dt0p0025"]["job_id"], "4242")

    def test_no_duplicate_submission(self):
        m = json.loads(ocv.MANIFEST.read_text())
        m["cases"]["cutLink-L1-dt0p0025"]["job_id"] = "77"
        ocv.save_manifest(m)
        with mock.patch.object(ocv, "queue_jobs",
                               return_value={"77": ("RUNNING", 1)}), \
                mock.patch.object(ocv.shutil, "which", return_value=None), \
                mock.patch.object(ocv.subprocess, "run",
                                  side_effect=AssertionError("submitted")):
            # Nothing to submit succeeds, with or without --yes, without
            # sbatch
            self.assertEqual(ocv.cmd_submit(self.S), 0)

            class NoYes(self.S):
                yes = False
            self.assertEqual(ocv.cmd_submit(NoYes), 0)

    def test_nonempty_submission_needs_sbatch(self):
        with mock.patch.object(ocv, "queue_jobs", return_value={}), \
                mock.patch.object(ocv.shutil, "which", return_value=None), \
                mock.patch.object(ocv.subprocess, "run",
                                  side_effect=AssertionError("submitted")):
            with self.assertRaisesRegex(RuntimeError, "sbatch is unavailable"):
                ocv.cmd_submit(self.S)

    def test_cpu_cap(self):
        with mock.patch.object(ocv, "queue_jobs",
                               return_value={"1": ("RUNNING", 96)}), \
                mock.patch.object(ocv.subprocess, "run",
                                  side_effect=AssertionError("submitted")):
            with self.assertRaises(RuntimeError):
                ocv.cmd_submit(self.S)

    def test_summary_without_yes(self):
        class N(self.S):
            yes = False
        with mock.patch.object(ocv, "queue_jobs", return_value={}), \
                mock.patch.object(ocv.subprocess, "run",
                                  side_effect=AssertionError("submitted")):
            self.assertEqual(ocv.cmd_submit(N), 0)


class Failures(Sandbox):
    def test_solver_problem(self):
        case = self.tmp / "c"
        case.mkdir()
        self.assertEqual(ocv.solver_problem(case), "no log.solids4Foam")
        (case / "log.solids4Foam").write_text("Time = 0.1\n")
        self.assertIn("did not complete", ocv.solver_problem(case))
        (case / "log.solids4Foam").write_text("--> FOAM FATAL ERROR\nEnd\n")
        self.assertIn("failed", ocv.solver_problem(case))
        (case / "log.solids4Foam").write_text("Time = 0.1\nEnd\n")
        self.assertEqual(ocv.solver_problem(case), "")
        (case / "log.solids4Foam").write_text("[1,0]<stdout>:End\n")
        self.assertEqual(ocv.solver_problem(case), "")
        (case / "log.solids4Foam").write_text(
            "trapFpe: Floating point exception trapping enabled "
            "(FOAM_SIGFPE).\nTime = 0.1\nEnd\n")
        self.assertEqual(ocv.solver_problem(case), "")
        (case / "log.solids4Foam").write_text(
            "Time = 0.1\n#0  Foam::error::printStack(Foam::Ostream&)\n"
            "#1  Foam::sigFpe::sigHandler(int)\n")
        self.assertIn("failed", ocv.solver_problem(case))

    def test_log_stats(self):
        case = self.tmp / "c"
        case.mkdir()
        (case / "log.solids4Foam").write_text(
            "Courant Number mean: 0.01 max: 0.3\nTime = 0.1\n"
            "time step continuity errors : sum local = 1e-9, global = 1e-12,"
            " cumulative = 2e-12\n"
            "[1,0]<stdout>:Courant Number mean: 0.02 max: 0.45\nTime = 0.2\n"
            "ExecutionTime = 12.5 s  ClockTime = 13 s\nEnd\n")
        s = ocv.log_stats(case)
        self.assertEqual(s["max_co"], 0.45)
        self.assertEqual(s["n_steps"], 2)
        self.assertEqual(s["execution_s"], 12.5)
        self.assertAlmostEqual(s["max_continuity_local"], 1e-9)

    def test_cleanup_scope(self):
        text = (ocv.VERIFICATION_DIR / "Allclean").read_text()
        command = [l for l in text.splitlines() if l.startswith("rm ")]
        self.assertEqual(command, ["rm -rf work postProcessing "
                                   "scripts/__pycache__"])
        self.assertIn('cd "${0%/*}" || exit 1', text)


RESULTS_DIR = ocv.REFERENCE_DIR / "results"


class Check(Sandbox):
    """./Allverify check against the tracked results"""

    def write_current(self, rows):
        ocv.OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
        keys = ["case", "integrity", "integrity_notes"] + \
            list(ocv.CHECK_METRICS)
        with (ocv.OUTPUT_DIR / "cases.csv").open("w", newline="") as h:
            w = csv.DictWriter(h, fieldnames=keys, extrasaction="ignore")
            w.writeheader()
            w.writerows(rows)

    def tracked(self, case):
        with (RESULTS_DIR / "cases.csv").open() as h:
            return next(dict(r, integrity_notes="")
                        for r in csv.DictReader(h) if r["case"] == case)

    def test_reproduced_case_passes(self):
        self.write_current([self.tracked("cutLink-L1-dt0p0025")])
        self.assertEqual(ocv.cmd_check(None), 0)

    def test_changed_lift_fails(self):
        r = self.tracked("cutLink-L1-dt0p0025")
        r["Cl_nrms_B"] = str(float(r["Cl_nrms_B"])*0.1)
        self.write_current([r])
        self.assertEqual(ocv.cmd_check(None), 1)

    def test_integrity_failure_fails(self):
        r = self.tracked("penalty-L1-dt0p0025")
        r["integrity"], r["integrity_notes"] = "FAIL", "no force file"
        self.write_current([r])
        self.assertEqual(ocv.cmd_check(None), 1)

    def scipy_absent(self):
        return mock.patch.object(ocv, "have_scipy", return_value=False)

    def test_nan_filtered_metric_skipped_without_scipy(self):
        r = self.tracked("cutLink-L1-dt0p0025")
        r["Cd_hf_rms_above_5Hz"] = "nan"
        self.write_current([r])
        with self.scipy_absent(), mock.patch("builtins.print") as out:
            self.assertEqual(ocv.cmd_check(None), 0)
        lines = "\n".join(str(c.args[0]) for c in out.call_args_list)
        self.assertIn("SKIP cutLink-L1-dt0p0025: Cd_hf_rms_above_5Hz", lines)
        self.assertIn("1 skipped (scipy)", lines)
        # The raw-force metrics were still compared
        self.assertIn("PASS cutLink-L1-dt0p0025: Cd_nrms_B", lines)

    def test_nan_filtered_metric_fails_with_scipy(self):
        r = self.tracked("cutLink-L1-dt0p0025")
        r["Cd_hf_rms_above_5Hz"] = "nan"
        self.write_current([r])
        with mock.patch.object(ocv, "have_scipy", return_value=True):
            self.assertEqual(ocv.cmd_check(None), 1)

    def test_nan_primary_metric_fails_even_without_scipy(self):
        for key in ("Cd_nrms_B", "Cl_nrms_B", "Cd_nrms_A", "Cl_nrms_A"):
            r = self.tracked("cutLink-L1-dt0p0025")
            r[key] = "nan"
            self.write_current([r])
            with self.scipy_absent():
                self.assertEqual(ocv.cmd_check(None), 1, key)

    def test_wrong_filtered_metric_still_fails_without_scipy(self):
        # Only a nonfinite value is skipped: a finite one is compared
        r = self.tracked("cutLink-L1-dt0p0025")
        r["Cd_hf_rms_above_5Hz"] = str(float(r["Cd_hf_rms_above_5Hz"])*3)
        self.write_current([r])
        with self.scipy_absent():
            self.assertEqual(ocv.cmd_check(None), 1)

    def test_integrity_failure_fails_without_scipy(self):
        r = self.tracked("penalty-L1-dt0p0025")
        r["integrity"], r["integrity_notes"] = "FAIL", "non-finite force"
        r["Cd_hf_rms_above_5Hz"] = "nan"
        self.write_current([r])
        with self.scipy_absent():
            self.assertEqual(ocv.cmd_check(None), 1)

    def test_filters_without_scipy_are_nan_not_numbers(self):
        t = [i*0.01 for i in range(1000)]
        y = [math.sin(i) for i in range(1000)]
        real = __import__
        def no_scipy(name, *a, **k):
            if name.startswith("scipy"):
                raise ImportError(name)
            return real(name, *a, **k)
        with mock.patch("builtins.__import__", side_effect=no_scipy):
            self.assertFalse(ocv.have_scipy())
            self.assertTrue(math.isnan(ocv.high_frequency_rms(t, y, 5.0)))
            self.assertEqual(ocv.low_pass(t, y, 5.0), ([], []))

    def test_untracked_case_not_checked(self):
        r = self.tracked("cutLink-L1-dt0p0025")
        r["case"] = "cutLink-L1-dt0p0025-clustered"
        self.write_current([r])
        self.assertEqual(ocv.cmd_check(None), 0)


class Provenance(Sandbox):
    def test_variant(self):
        self.assertEqual(ocv.openfoam_variant("v2512"), "OpenFOAM.com")
        self.assertEqual(ocv.openfoam_variant("12"), "openfoam.org")
        self.assertEqual(ocv.openfoam_variant("", ""), "")
        self.assertEqual(ocv.openfoam_variant("4.1", "extend"), "extend")

    def make_build(self, with_info=True):
        root = self.tmp / "build"
        exe = root / "platforms/linux64GccDPInt32Opt/bin/solids4Foam"
        exe.parent.mkdir(parents=True)
        exe.write_bytes(b"binary")
        if with_info:
            (root / "BUILD_INFO.txt").write_text(
                "SHA=0123456789abcdef\nWM_PROJECT_VERSION=v2512\n")
        return exe

    def test_build_with_record(self):
        exe = self.make_build()
        b = ocv.solver_build(str(exe))
        self.assertEqual(b["solver_build"], "0123456789abcdef")
        self.assertEqual(len(b["solver_binary_sha256"]), 64)

    def test_build_without_record_is_not_invented(self):
        exe = self.make_build(with_info=False)
        b = ocv.solver_build(str(exe))
        self.assertEqual(b["solver_build"], "")
        self.assertEqual(len(b["solver_binary_sha256"]), 64)
        with mock.patch.object(ocv.shutil, "which", return_value=None):
            self.assertEqual(ocv.solver_build(), {
                "solver_executable": "", "solver_binary_sha256": "",
                "solver_build": ""})

    def test_case_json_records_provenance_separately(self):
        env = {"WM_PROJECT_VERSION": "v2512"}
        with mock.patch.dict(ocv.os.environ, env), \
                mock.patch.object(ocv.shutil, "which", return_value=None):
            class A:
                runs = ["smoke"]
                ranks = 0
                force = False
                traction_quadrature = "greedy"
            ocv.cmd_generate(A)
        meta = json.loads((ocv.CASES_DIR / "cutLink-L1-dt0p0025-t0p5"
                           / "case.json").read_text())
        p = meta["provenance"]
        self.assertEqual(p["openfoam_version"], "v2512")
        self.assertEqual(p["openfoam_variant"], "OpenFOAM.com")
        self.assertEqual(p["solver_build"], "")
        self.assertEqual(p["tutorial_revision"], meta["git_sha"])
        self.assertEqual(meta["tracked_baseline"],
                         {"openfoam": "v2412", "solids4foam": "f33b244a2"})
        self.assertEqual(set(ocv.PROVENANCE_COLUMNS) - set(p), set())

    def test_check_reports_solver_of_the_run(self):
        r = next(dict(r, integrity_notes="", openfoam_version="v2512",
                      solver_build="abc123")
                 for r in csv.DictReader((RESULTS_DIR / "cases.csv").open())
                 if r["case"] == "cutLink-L1-dt0p0025")
        ocv.OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
        keys = ["case", "integrity", "integrity_notes", "openfoam_version",
                "solver_build"] + list(ocv.CHECK_METRICS)
        with (ocv.OUTPUT_DIR / "cases.csv").open("w", newline="") as h:
            w = csv.DictWriter(h, fieldnames=keys, extrasaction="ignore")
            w.writeheader()
            w.writerow(r)
        with mock.patch("builtins.print") as out:
            self.assertEqual(ocv.cmd_check(None), 0)
        text = "\n".join(str(c.args[0]) for c in out.call_args_list)
        self.assertIn("OpenFOAM v2512, solids4foam build abc123; tracked "
                      "results: OpenFOAM v2412, solids4foam f33b244a2", text)


class Export(Sandbox):
    """./Allverify export: which histories are written"""

    def setUp(self):
        super().setUp()
        patch = mock.patch.object(ocv, "REFERENCE_DIR", self.tmp / "reference")
        patch.start()
        self.addCleanup(patch.stop)
        ocv.OUTPUT_DIR.mkdir(parents=True)
        cases = [("cutLink-L2-dt0p00125", "R2A", 8.0, 2),
                 ("cutLink-L2-dt0p00125-np1", "O4", 8.0, 2),
                 ("cutLink-L1-dt0p0025-t0p5", "smoke", 0.5, 1)]
        with (ocv.OUTPUT_DIR / "cases.csv").open("w", newline="") as h:
            w = csv.writer(h)
            w.writerow(["case", "run", "method", "level", "deltaT",
                        "end_time", "integrity", "openfoam_version"])
            for name, run, end, level in cases:
                w.writerow([name, run, "cutLink", level, 0.00125, end,
                            "PASS", "v2512"])
                with (ocv.OUTPUT_DIR / f"history-{name}.csv").open(
                        "w", newline="") as hh:
                    hw = csv.writer(hh)
                    hw.writerow(["t", "Cd_total"])
                    for i in range(0, 1601):
                        hw.writerow([f"{i*0.005:.10g}", "1"])

    def exported(self):
        out = ocv.REFERENCE_DIR / "results"
        return sorted(p.name for p in out.iterdir())

    def test_histories_of_figure_runs_only(self):
        class A:
            full_resolution = False
        ocv.cmd_export(A)
        self.assertEqual(self.exported(), [
            "cases.csv", "history-cutLink-L2-dt0p00125.csv"])
        with (ocv.REFERENCE_DIR / "results/cases.csv").open() as h:
            rows = list(csv.DictReader(h))
        # All the cases keep their rows; provenance only if recorded
        self.assertEqual(len(rows), 3)
        self.assertEqual(rows[0]["openfoam_version"], "v2512")

    def test_full_resolution_is_optional(self):
        class A:
            full_resolution = True
        ocv.cmd_export(A)
        self.assertIn("zoom-cutLink-L2-dt0p00125.csv", self.exported())


def markdown_table(text: str, after: str) -> list:
    """Body rows of the first table after the text 'after' in results.md,
    as lists of cells"""
    rows = []
    for line in text.split(after, 1)[1].splitlines():
        if line.startswith("|"):
            rows.append([c.strip() for c in line.strip("|").split("|")])
        elif rows:
            break
    return rows[2:]


def rounding(cell: str) -> float:
    """Half a unit of the last digit of a table cell"""
    digits = cell.rstrip("%").split(".")
    return 0.5*10**-len(digits[1]) if len(digits) > 1 else 0.5


def percent(cell: str) -> tuple:
    """A percentage cell as a fraction, and its rounding"""
    return float(cell.rstrip("%"))/100, rounding(cell)/100 + 1e-12


class TrackedResults(unittest.TestCase):
    """The tracked references, compact results, results.md and figures
    agree"""

    @classmethod
    def setUpClass(cls):
        with (RESULTS_DIR / "cases.csv").open() as handle:
            cls.header = next(csv.reader(handle))
        with (RESULTS_DIR / "cases.csv").open() as handle:
            cls.rows = list(csv.DictReader(handle))
        cls.by_run = {r["run"]: r for r in cls.rows}
        cls.text = (ocv.VERIFICATION_DIR / "results.md").read_text()

    def test_references(self):
        refs = ocv.references()
        self.assertEqual(set(refs), {"A", "B", "digitised"})
        t, _ = refs["B"]["Cd"]
        self.assertEqual((len(t), t[0], t[-1]), (3200, 0.0025, 8.0))
        self.assertAlmostEqual(refs["A"]["Cd"][0][-1], 7.58)
        self.assertEqual(len(refs["digitised"]["Cd"][0]), 449)

    def test_cases_csv(self):
        self.assertEqual(self.header, list(ocv.EXPORT_COLUMNS))
        hashes = {}
        for r in self.rows:
            self.assertEqual(r["integrity"], "PASS")
            self.assertAlmostEqual(float(r["t_final"]), float(r["end_time"]))
            # The histories of the 8 s runs that the figures draw (not the
            # serial repeat O4 and the shorter runs, whose metrics are here)
            drawn = (float(r["end_time"]) == ocv.END_TIME
                     and r["run"] not in ocv.HISTORY_NOT_TRACKED)
            self.assertEqual(
                (RESULTS_DIR / f"history-{r['case']}.csv").exists(), drawn,
                r["case"])
            hashes.setdefault(r["level"], set()).add(r["mesh_hash"])
        self.assertTrue(all(len(h) == 1 for h in hashes.values()))

    def test_runs_table(self):
        rows = markdown_table(self.text, "## Runs")
        self.assertEqual({r[0] for r in rows}, set(self.by_run))
        for run, method, _, dt, end, ranks, co, cpu in rows:
            r = self.by_run[run]
            self.assertEqual(method, r["method"])
            self.assertEqual(float(dt), float(r["deltaT"]))
            self.assertEqual(float(end), float(r["end_time"]))
            self.assertEqual(int(ranks), int(r["ranks"]))
            self.assertAlmostEqual(float(co), float(r["max_co"]),
                                   delta=rounding(co) + 1e-12)
            self.assertAlmostEqual(float(cpu), float(r["cpu_hours"]),
                                   delta=rounding(cpu) + 1e-12)

    def test_agreement_tables(self):
        columns = ("Cd_nrms_A", "Cd_nrms_B", "Cd_nrms_B_lowpass5Hz",
                   "Cl_nrms_B", "Cl_nrms_B_lowpass5Hz")
        rows = markdown_table(self.text, "## Agreement with the references")
        full = {r["run"] for r in self.rows
                if float(r["end_time"]) == ocv.END_TIME and r["run"] != "O4"}
        self.assertEqual({r[0] for r in rows}, full)
        for row in rows:
            for cell, column in zip(row[1:], columns):
                value, tol = percent(cell)
                self.assertAlmostEqual(
                    value, float(self.by_run[row[0]][column]), delta=tol,
                    msg=(row[0], column))
        rows = markdown_table(self.text, "columns 11-13, from the same runs")
        self.assertEqual({r[0] for r in rows},
                         {r for r in full if self.by_run[r]["method"]
                          == "cutLink"})
        for run, cd, cl, hf in rows:
            r = self.by_run[run]
            for cell, column in ((cd, "Cd_ME_nrms_B"), (cl, "Cl_ME_nrms_B")):
                value, tol = percent(cell)
                self.assertAlmostEqual(value, float(r[column]), delta=tol,
                                       msg=(run, column))
            self.assertAlmostEqual(
                float(hf), float(r["Cd_ME_hf_rms_above_5Hz"]),
                delta=rounding(hf) + 1e-12, msg=run)

    def test_fluctuation_table(self):
        runs = {"cutLink": ("R1A", "R2A", "R3A", "R5A"),
                "ghostCell": ("R1D", "R2D", "R3D", "R5D"),
                "penalty": ("R1B", "R2B", "R3B", "R5B"),
                "legacy incremental": ("R1C", "R2C", "R3C", None)}
        rows = markdown_table(self.text, "Its RMS in Cd is:")
        self.assertEqual(len(rows), 4)
        for method, *cells in rows:
            for cell, run in zip(cells, runs[method]):
                if run is None:
                    continue
                value = float(self.by_run[run]["Cd_hf_rms_above_5Hz"])
                self.assertAlmostEqual(float(cell), value,
                                       delta=rounding(cell) + 1e-12,
                                       msg=(method, run))

    def test_figures_regenerate(self):
        try:
            import matplotlib  # noqa: F401
            import scipy  # noqa: F401
        except ImportError:
            self.skipTest("matplotlib or scipy is not installed")
        out = Path(tempfile.mkdtemp())
        self.addCleanup(shutil.rmtree, out, True)
        subprocess.run(
            [sys.executable, str(SCRIPT.parent / "plot_verification.py"),
             "--data", str(RESULTS_DIR), "--out", str(out), "--dpi", "40",
             "--png-only"], check=True, capture_output=True, timeout=300)
        # cell_crossing needs the full-resolution extracts, which are not
        # tracked (./Allverify export --full-resolution)
        self.assertFalse(list(RESULTS_DIR.glob("zoom-*.csv")))
        self.assertEqual(
            sorted(p.name for p in out.iterdir()),
            sorted(p.name for p in
                   (ocv.VERIFICATION_DIR / "images").glob("*.png")
                   if p.name != "cell_crossing.png"))


if __name__ == "__main__":
    unittest.main()
