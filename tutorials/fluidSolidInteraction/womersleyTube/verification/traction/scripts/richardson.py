"""Observed order and Richardson extrapolation for complex errors on three
levels with refinement ratio 2.

For complex values Richardson is meaningful only if
r = (e1 - e2)/(e2 - e4) is close to a positive real: |r| = 2^p gives the
order, arg(r) measures the departure from an asymptotic (single-power,
same-phase) series. The order is accepted only if
  - e1 != e2 and e2 != e4 (finite, non-zero differences),
  - MIN_ORDER < p (a near-zero or negative p makes 1/(2^p - 1) blow up),
  - |arg r| <= MAX_ARG_DEG (otherwise the series is non-asymptotic).
Otherwise the observed-order limit falls back to the formal-order limit and
the result is flagged.
"""
import cmath
import math

MIN_ORDER = 0.3
MAX_ARG_DEG = 30.0


def observed_order(e1, e2, e4, min_order=MIN_ORDER, max_arg_deg=MAX_ARG_DEG):
    """Returns dict: order (nan if undefined), arg_deg (nan if undefined),
    valid (bool), reason ('' if valid)."""
    d1, d2 = e1 - e2, e2 - e4
    out = {"order": float("nan"), "arg_deg": float("nan"), "valid": False, "reason": ""}
    if not all(map(cmath.isfinite, (e1, e2, e4))):
        out["reason"] = "non-finite value"
        return out
    if abs(d1) == 0 or abs(d2) == 0:
        out["reason"] = "zero difference between levels"
        return out
    ratio = d1 / d2
    out["order"] = math.log(abs(ratio)) / math.log(2)
    out["arg_deg"] = math.degrees(cmath.phase(ratio))
    if not out["order"] > min_order:
        out["reason"] = f"observed order {out['order']:.2f} <= {min_order:g}"
    elif abs(out["arg_deg"]) > max_arg_deg:
        out["reason"] = f"non-asymptotic: |arg ratio| = {abs(out['arg_deg']):.0f} deg > {max_arg_deg:g}"
    else:
        out["valid"] = True
    return out


def richardson_limits(e1, e2, e4, formal_order=2.0):
    """(info, limit_formal, limit_observed). limit_observed equals limit_formal
    when the observed order is not valid (see observed_order)."""
    info = observed_order(e1, e2, e4)
    lim_formal = e4 + (e4 - e2) / (2.0 ** formal_order - 1.0)
    if info["valid"]:
        lim_obs = e4 + (e4 - e2) / (2.0 ** info["order"] - 1.0)
    else:
        lim_obs = lim_formal
    return info, lim_formal, lim_obs
