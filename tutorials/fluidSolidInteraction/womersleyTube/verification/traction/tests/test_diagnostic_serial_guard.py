"""Static check (no OpenFOAM needed) that the coded function object stops in a
parallel run before it opens any output file. The compiled behaviour has to be
confirmed in OpenFOAM (see the README)."""
import re
import unittest
from pathlib import Path

SRC = (Path(__file__).resolve().parents[1] / "scripts/interfaceTractionDiagnostic").read_text()


class SerialGuard(unittest.TestCase):
    def test_guard_first_in_codeExecute(self):
        body = SRC[SRC.index("codeExecute"):]
        guard = re.search(r"if\s*\(\s*Pstream::parRun\(\)\s*\)\s*\{[^}]*FatalErrorInFunction[^}]*"
                          r"exit\(FatalError\)", body, re.S)
        self.assertIsNotNone(guard)
        for token in ("OFstream(", "mkDir(", "fluidFile_.reset", "mesh()"):
            self.assertLess(guard.start(), body.index(token), token)

    def test_only_one_codeExecute(self):
        self.assertEqual(SRC.count("codeExecute"), 1)


if __name__ == "__main__":
    unittest.main()
