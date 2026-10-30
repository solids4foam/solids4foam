#!/usr/bin/env python3
"""One ordered definition of the status severities, shared by acceptance.py, family.py and verify.py.

Order, from least to most severe:
  0  PASS
  1  INFO                       informational item (no claim)
  2  NOT-VERIFIED               not checked here: never a PASS
  3  CONDITIONAL, DISCLOSED, ABOVE-FLOOR   documented difference, disclosed limitation, quantity above the disclosed floor
  4  UNKNOWN                    an unrecognised status string: non-PASS, handled safely, blocks the exit code
  5  INVALID-INPUT              inputs inconsistent: the assessment is meaningless
  6  INCOMPLETE                 missing, truncated or unreadable input
  7  FAIL                       a criterion is violated
worst() of statuses never returns PASS unless every status is PASS; worst() of nothing is NOT-VERIFIED (nothing was
evaluated, so no PASS). An unrecognised status never raises: it is reported as UNKNOWN.
"""
import sys as _sys
_sys.dont_write_bytecode = True

ORDER = {"PASS": 0, "INFO": 1, "NOT-VERIFIED": 2, "CONDITIONAL": 3, "DISCLOSED": 3, "ABOVE-FLOOR": 3,
         "UNKNOWN": 4, "INVALID-INPUT": 5, "INCOMPLETE": 6, "FAIL": 7}
BLOCKING = ("UNKNOWN", "INVALID-INPUT", "INCOMPLETE", "FAIL")     # layer A: exit code 1
NEUTRAL = ("INFO", "NOT-VERIFIED")                                # reported, but they do not decide the status of a run


def norm(status):
    """Recognised status, else UNKNOWN."""
    return status if (isinstance(status, str) and status in ORDER) else "UNKNOWN"


def severity(status):
    return ORDER[norm(status)]


def worst(statuses):
    """The most severe status; empty -> NOT-VERIFIED; PASS only if all are PASS."""
    out = None
    for s in statuses:
        s = norm(s)
        if out is None or ORDER[s] > ORDER[out]:
            out = s
    return "NOT-VERIFIED" if out is None else out


def decide(statuses):
    """Status of a run from its checks: INFO and NOT-VERIFIED items are reported but do not decide;
    nothing left to decide on -> NOT-VERIFIED (never PASS)."""
    return worst([norm(s) for s in statuses if norm(s) not in NEUTRAL])


def cap_pass(status, ceiling):
    """A PASS is lowered to `ceiling` (a family with documented differences cannot give PASS)."""
    return ceiling if (status == "PASS" and severity(ceiling) > 0) else status
