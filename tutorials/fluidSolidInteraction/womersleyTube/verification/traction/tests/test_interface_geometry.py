"""womersley_interface_traction.py: undeformed positions, face-area weights and
full-coordinate face matching, on synthetic faces and synthetic .dat files."""
import math
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import womersley_interface_traction as W  # noqa: E402

L, R, PERIOD = 15.0, 1.0, 1.0
OMEGA = 2 * math.pi / PERIOD
K = 0.144 - 0.01j            # complex wave number, exp(-i k x)
UX = 7.3e-4                  # axial interface motion amplitude


def faces(xs, n_circ=1):
    """face centres on a cylinder of radius R: n_circ faces per axial station"""
    out, area = [], []
    for i, x in enumerate(xs):
        for j in range(n_circ):
            th = 2 * math.pi * (j + 0.5) / n_circ
            out.append((x, R * math.cos(th), R * math.sin(th)))
            area.append((xs[min(i + 1, len(xs) - 1)] - xs[max(i - 1, 0)]) / 2 / n_circ * 2 * math.pi * R)
    return np.array(out), np.array(area)


def write_dat(path, rows):
    np.savetxt(path, np.array(rows), fmt="%.15g")


def make_case(case: Path, xs, amp, n_circ=1, d0=UX, solid_shuffle=False):
    """Fluid/solid .dat files for P = Re[amp exp(i(w t - k x0))] carried by
    material faces. The fluid mesh moves by D(t) - D(0) with D(t) = UX cos(w t)
    and D(0) = d0: d0 = UX is a start from the periodic state (the mesh offset is
    zero at every whole period); d0 = 0 is a start from the undisplaced state
    (the mesh then oscillates about +UX, offset UX at a whole period); d0 = None
    keeps the fluid mesh fixed. The solid mesh never moves."""
    X0, area = faces(xs, n_circ)
    nf = len(X0)
    nvec = np.zeros((nf, 3))
    nvec[:, 1:] = X0[:, 1:] / R
    steps = 20
    ts = [(i + 1) * PERIOD / steps for i in range(4 * steps)]
    fl, so = [], []
    perm = np.random.RandomState(0).permutation(nf) if solid_shuffle else np.arange(nf)
    for t in ts:
        x = X0.copy()
        if d0 is not None:
            x[:, 0] += UX * math.cos(OMEGA * t) - d0   # mesh moves by D(t) - D(0)
        for f in range(nf):
            p = (amp * np.exp(1j * (OMEGA * t - K * X0[f, 0]))).real
            row = [t, f, *x[f], area[f], *nvec[f], p] + [0.0] * 12 + [0.0, 1.0]
            fl.append(row)
        for g in range(nf):
            f = perm[g]
            so.append([t, g, *X0[f], area[f], *nvec[f]] + [0.0] * 6)
    (case / "postProcessing/interfaceTraction").mkdir(parents=True)
    write_dat(case / "postProcessing/interfaceTraction/fluid.dat", fl)
    write_dat(case / "postProcessing/interfaceTraction/solid.dat", so)


EX = {"rho": 1.0, "period": PERIOD, "omega": OMEGA, "k": K, "L": L}


class Projection(unittest.TestCase):
    def test_projection_does_not_depend_on_the_fluid_mesh_motion(self):
        # corrected code: exact for any D(0), including a mesh that oscillates
        # about an offset position, and for a fixed fluid mesh
        xs = np.linspace(0.1, L - 0.1, 41)
        for d0 in (UX, 0.0, None):
            with self.subTest(d0=d0), tempfile.TemporaryDirectory() as tmp:
                case = Path(tmp)
                make_case(case, xs, 0.4 + 0.3j, d0=d0)
                res = W.analyse_case(case, EX, 2)
            self.assertLess(abs(res["P"] - (0.4 + 0.3j)), 5e-6)

    def test_fluid_window_mean_positions_carry_the_static_offset(self):
        # the alternative of projecting on the window-mean fluid centres (kept
        # as --fluid-mean-positions) is displaced by -D(0): phase k UX ~ 1e-4
        xs = np.linspace(0.1, L - 0.1, 41)
        with tempfile.TemporaryDirectory() as tmp:
            case = Path(tmp)
            make_case(case, xs, 0.4 + 0.3j, d0=UX)
            fluid_mean = W.analyse_case(case, EX, 2, fluid_mean_positions=True)
        self.assertGreater(abs(fluid_mean["P"] - (0.4 + 0.3j)), 2e-5)


def old_projection_P(case: Path, ex: dict, window: int) -> complex:
    """Faithful copy of what the analysis did at 7c9f1e2cf (the PR head) for
    the pressure series: positions xf = the deformed fluid face centre at the
    LAST sample, solid faces matched by nearest axial position with a 1e-6 L
    tolerance, weights from midpoints between the sorted axial positions. Only
    the lines that decide the projection are kept."""
    tf, F = W.per_time(W.read(case / "postProcessing/interfaceTraction/fluid.dat", 24), 24)
    ts, S = W.per_time(W.read(case / "postProcessing/interfaceTraction/solid.dat", 15), 15)
    t_end = tf[-1]
    t0 = t_end - window * ex["period"]
    xf = F[-1, :, 2]
    xs = S[-1, :, 2]
    match = np.array([np.argmin(np.abs(xs - x)) for x in xf])
    if np.max(np.abs(xs[match] - xf)) > 1e-6 * ex["L"]:
        raise RuntimeError("solid and fluid interface faces do not coincide")
    sel = (xf >= 0.25 * ex["L"]) & (xf <= 0.75 * ex["L"])
    order = np.argsort(xf)
    xs_sorted = xf[order]
    edges = np.concatenate([[xs_sorted[0] - 0.5 * (xs_sorted[1] - xs_sorted[0])],
                            0.5 * (xs_sorted[1:] + xs_sorted[:-1]),
                            [xs_sorted[-1] + 0.5 * (xs_sorted[-1] - xs_sorted[-2])]])
    dx = np.empty_like(xf)
    dx[order] = np.diff(edges)
    w = np.where(sel, dx, 0.0)
    proj = np.exp(1j * ex["k"] * xf)
    a, _ = W.harmonic(tf, ex["rho"] * F[:, :, 9], ex["omega"], t0, t_end)
    return complex(np.sum(a * proj * w) / np.sum(w))


class HistoricalBehaviour(unittest.TestCase):
    """What is and is not claimed about the code at 7c9f1e2cf.

    NOT claimed: that the old projection silently produced a ~1e-4 phase error
    in the PR's runs. The old code took the deformed fluid face centres at the
    LAST sample, a whole number of periods after the start. If the run starts
    from the periodic state (D(0) equal to the periodic value at that phase,
    d0 = UX here) the mesh offset D(t_end) - D(0) is ~0 there, and the old and
    the corrected code agree (test_old_code_is_fine_when_the_last_sample_is_
    undisplaced). That D(0) matches the periodic state in the real runs has not
    been checked (no stored output on this machine).

    Claimed: the old result depended on the fluid mesh position at one
    sample. The old safety net was a 1e-6 L match tolerance, which allows a
    phase error of only k x 1.5e-5 m ~ 2e-6; a larger offset (D(0) not equal
    to the periodic value: here up to UX = 7.3e-4 m, the size Philip quotes)
    makes the old code RAISE rather than silently be wrong, so it could not
    carry the 1e-4 error either; the corrected code needs no such coincidence.
    Also, the old nearest-axial-position matching and sorted-x midpoint weights
    break down with several faces per axial station.
    """

    def test_old_code_is_fine_when_the_last_sample_is_undisplaced(self):
        xs = np.linspace(0.1, L - 0.1, 41)
        with tempfile.TemporaryDirectory() as tmp:
            case = Path(tmp)
            make_case(case, xs, 0.4 + 0.3j, d0=UX)
            old = old_projection_P(case, EX, 2)
            new = W.analyse_case(case, EX, 2)["P"]
        self.assertLess(abs(old - (0.4 + 0.3j)), 5e-6)
        self.assertLess(abs(new - (0.4 + 0.3j)), 5e-6)

    def test_old_code_stops_when_the_fluid_mesh_is_displaced_at_the_last_sample(self):
        # D(0) = 0: the mesh oscillates about +UX; last sample displaced by
        # UX = 7.3e-4 m >> 1e-6 L = 1.5e-5 m. Old: raises. New: correct.
        xs = np.linspace(0.1, L - 0.1, 41)
        with tempfile.TemporaryDirectory() as tmp:
            case = Path(tmp)
            make_case(case, xs, 0.4 + 0.3j, d0=0.0)
            with self.assertRaises(RuntimeError):
                old_projection_P(case, EX, 2)
            new = W.analyse_case(case, EX, 2)["P"]
        self.assertLess(abs(new - (0.4 + 0.3j)), 5e-6)

    def test_old_code_phase_error_is_bounded_by_its_tolerance(self):
        # largest offset the old tolerance accepts: error <~ k * 1.5e-5 ~ 2e-6,
        # i.e. the old code could not produce a 1e-4 error without raising
        xs = np.linspace(0.1, L - 0.1, 41)
        with tempfile.TemporaryDirectory() as tmp:
            case = Path(tmp)
            make_case(case, xs, 0.4 + 0.3j, d0=UX - 1.4e-5)
            old = old_projection_P(case, EX, 2)
        err = abs(old - (0.4 + 0.3j))
        self.assertLess(err, 1e-5)
        self.assertGreater(err, 1e-7)     # the offset is visible, but small


class ProjectionGeometry(unittest.TestCase):
    def test_moving_solid_mesh_is_refused(self):
        xs = np.linspace(0.1, L - 0.1, 11)
        with tempfile.TemporaryDirectory() as tmp:
            case = Path(tmp)
            make_case(case, xs, 0.4 + 0.3j)
            f = case / "postProcessing/interfaceTraction/solid.dat"
            d = np.loadtxt(f)
            d[:, 2] += 1e-3 * np.sin(d[:, 0])
            write_dat(f, d)
            with self.assertRaises(RuntimeError):
                W.analyse_case(case, EX, 2)
            W.analyse_case(case, EX, 2, fluid_mean_positions=True)   # explicit opt-in

    def test_several_faces_per_station_and_shuffled_solid_order(self):
        xs = np.linspace(0.1, L - 0.1, 25)
        with tempfile.TemporaryDirectory() as tmp:
            case = Path(tmp)
            make_case(case, xs, 0.4 + 0.3j, n_circ=4, solid_shuffle=True)
            res = W.analyse_case(case, EX, 2)
        self.assertLess(abs(res["P"] - (0.4 + 0.3j)), 5e-6)
        self.assertEqual(res["nFaces"], int(np.sum((xs >= 0.25 * L) & (xs <= 0.75 * L))) * 4)

    def test_face_area_weights_nonuniform_axial_spacing(self):
        xs = np.cumsum(np.r_[0.05, np.linspace(0.05, 0.6, 40)])
        xs = xs * (L - 0.1) / xs[-1]
        # a field that is not a pure wave on the faces: weights matter; compare
        # the projection with the exact area-weighted sum
        X, area = faces(xs)
        a = np.exp(-1j * K * xs) * (1 + 0.02 * np.sin(xs))
        sel = (xs >= 0.25 * L) & (xs <= 0.75 * L)
        w = np.where(sel, area, 0.0)
        got = W.wave_projection(a, xs, w, K)
        want = np.sum(a * np.exp(1j * K * xs) * area * sel) / np.sum(area * sel)
        self.assertAlmostEqual(abs(got - want), 0, 14)
        unweighted = W.wave_projection(a, xs, sel.astype(float), K)
        self.assertGreater(abs(unweighted - want), 1e-6)   # equal weights would differ


class Matching(unittest.TestCase):
    def setUp(self):
        self.Xf, _ = faces(np.linspace(0.1, L - 0.1, 11), n_circ=3)

    def test_matches_shuffled_faces_with_small_offset(self):
        perm = np.random.RandomState(1).permutation(len(self.Xf))
        Xs = (self.Xf + [7e-4, 0, 0])[perm]
        m = W.match_faces(self.Xf, Xs)
        np.testing.assert_array_equal(Xs[m], self.Xf + [7e-4, 0, 0])

    def test_axial_position_alone_would_be_ambiguous(self):
        # same x, different circumferential position: x-only nearest match would
        # give duplicates; the full coordinates give a one-to-one match
        Xs = self.Xf[::-1].copy()
        m = W.match_faces(self.Xf, Xs)
        self.assertEqual(np.unique(m).size, len(self.Xf))

    def test_offset_larger_than_quarter_spacing_raises(self):
        spacing = np.linalg.norm(self.Xf[0] - self.Xf[1])
        with self.assertRaises(RuntimeError):
            W.match_faces(self.Xf, self.Xf + [0.3 * spacing, 0, 0])

    def test_duplicate_match_and_count_mismatch_raise(self):
        Xs = self.Xf.copy()
        Xs[1] = Xs[0] + [1e-9, 0, 0]
        with self.assertRaises(RuntimeError):
            W.match_faces(self.Xf, Xs)
        with self.assertRaises(RuntimeError):
            W.match_faces(self.Xf, self.Xf[:-1])

    def test_repeated_axial_stations_without_circumferential_info_raise_in_analysis(self):
        xs = np.linspace(0.1, L - 0.1, 11)
        with tempfile.TemporaryDirectory() as tmp:
            case = Path(tmp)
            make_case(case, xs, 0.4 + 0j, n_circ=2)
            # corrupt the solid file: collapse y,z so faces coincide in 3D
            f = case / "postProcessing/interfaceTraction/solid.dat"
            d = np.loadtxt(f)
            d[:, 3:5] = 0.0
            write_dat(f, d)
            with self.assertRaises(RuntimeError):
                W.analyse_case(case, EX, 2)


if __name__ == "__main__":
    unittest.main()
