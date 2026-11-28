"""read_samples/analyse: the profile-sample movement tolerance is an argument
(replaces the string patching and exec of the driver in amplitude_series.py).

Run: python3 -m unittest discover -s tests  (from verification/traction)
"""
import sys
import tempfile
import types
import unittest
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "scripts"))
import womersley_tube_verification as V  # noqa: E402


class ReadSamplesTolerance(unittest.TestCase):
    def make_case(self, root, move):
        period, omega = 1.0, 2 * 3.141592653589793
        exact = types.SimpleNamespace(period=period, omega=omega)
        times = V.period_times(0.0, period, V.SAMPLES_PER_PERIOD)
        coord = [0.1 * i for i in range(5)]
        for name in ("profile", "axial"):
            for j, t in enumerate(times):
                d = root / V.SAMPLES.format(name) / f"{t:.12g}"
                d.mkdir(parents=True, exist_ok=True)
                # sample 2 moves by `move` x spacing at odd samples
                rows = [f"{c + (move * 0.1 if (i == 2 and j % 2) else 0):.9g} 1 2 3 4"
                        for i, c in enumerate(coord)]
                (d / f"{name}_p_U.xy").write_text("\n".join(rows) + "\n")
        return exact

    def read(self, move, **kw):
        with tempfile.TemporaryDirectory() as tmp:
            exact = self.make_case(Path(tmp), move)
            return V.read_samples(Path(tmp), exact, 1, 1,
                                  {"profile": 5, "axial": 5}, **kw)

    def test_default_rejects_large_movement(self):
        with self.assertRaises(RuntimeError):
            self.read(0.2)

    def test_relaxed_tolerance_accepts_it(self):
        self.assertEqual(len(self.read(0.2, profile_tolerance=0.5)["profile"]["coord"]), 5)

    def test_default_accepts_small_movement(self):
        self.read(0.01)

    def test_analyse_signature_exposes_argument(self):
        import inspect
        self.assertEqual(inspect.signature(V.analyse).parameters["profile_tolerance"].default, 0.05)


if __name__ == "__main__":
    unittest.main()
