"""Example solver: probe the black box, then fit a hypothesised model.

Strategy (what an expert would do after exploratory probing):
  1. value_d takes only 1 or 2 and scales everything like z^2 or more; the output is
     negative, -> 0 as value_b -> 0 and falls like -sqrt(value_b) at low value_b:
     log10 of a mean ionic activity coefficient of a z:z salt at concentration value_b.
  2. For z = 1 the value_a dependence has the extended Debye-Hueckel form
     -A sqrt(I)/(1 + B a sqrt(I)) with a in Angstrom, EXCEPT below a temperature-dependent
     kink near 3.5-4 A, where the output drops extra: ion pairs switch on at Bjerrum's
     critical distance q = z^2 e^2 / (8 pi eps kT). For z = 2 (q ~ 14 A) pairing is always on,
     and the size of the drop matches K_A = 4 pi N_A l^3 Q(b) with the stoichiometric
     gamma = alpha * gamma_free.
  3. value_c (0..100) acts only through eps_r(T) * T: temperature in degC in water.
     eps_r(t) is recovered from probes in the pairing-free regime (z = 1, a = 8 A),
     where the output is a monotone function of eps_r, then fitted by a cubic in t.

This solver is illustrative, not guaranteed optimal. Standard library only.
"""
import math

PROBE_BUDGET = 200

E, KB, NA, EPS0 = 1.602176634e-19, 1.380649e-23, 6.02214076e23, 8.8541878188e-12
T_ZERO = 273.15
TEMPS = (0.0, 12.5, 25.0, 37.5, 50.0, 62.5, 75.0, 87.5, 100.0)

_coef = None


def _electro(t_c, z, eps_r):
    ekt = EPS0 * eps_r * KB * (t_c + T_ZERO)
    length = z * z * E * E / (4 * math.pi * ekt)
    b_dh = math.sqrt(2 * NA * 1000.0 * E * E / ekt)
    a_dh = E * E * b_dh / (8 * math.pi * ekt) / math.log(10)
    return length, a_dh, b_dh


def _q(b, n=2000):
    if b <= 2:
        return 0.0
    h = (b - 2) / n
    f = lambda x: math.exp(x) / x ** 4
    s = f(2) + f(b) + sum((4 if i % 2 else 2) * f(2 + i * h) for i in range(1, n))
    return s * h / 3


def _model(a_ang, c, t_c, z, eps_r):
    length, a_dh, b_dh = _electro(t_c, z, eps_r)
    a = a_ang * 1e-10
    k_a = 4 * math.pi * NA * length ** 3 * _q(length / a) * 1000.0

    def lg(al):
        s = math.sqrt(al * c * z * z)
        return -a_dh * z * z * s / (1 + b_dh * a * s)

    lo, hi = 0.0, 1.0
    for _ in range(100):
        mid = 0.5 * (lo + hi)
        if k_a * mid * mid * c * (10 ** lg(mid)) ** 2 - (1 - mid) > 0:
            hi = mid
        else:
            lo = mid
    al = 0.5 * (lo + hi)
    return math.log10(al) + lg(al)


def _invert(target, f, lo, hi):
    """Bisection for x with f(x) = target, f monotone increasing."""
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        if f(mid) < target:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def _lstsq(xs, ys, deg):
    """Polynomial least squares via normal equations (Gaussian elimination)."""
    n = deg + 1
    m = [[sum(x ** (i + j) for x in xs) for j in range(n)] for i in range(n)]
    v = [sum(y * x ** i for x, y in zip(xs, ys)) for i in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            f = m[j][i] / m[i][i]
            m[j] = [mj - f * mi for mj, mi in zip(m[j], m[i])]
            v[j] -= f * v[i]
    x = [0.0] * n
    for i in reversed(range(n)):
        x[i] = (v[i] - sum(m[i][j] * x[j] for j in range(i + 1, n))) / m[i][i]
    return x


def _eps(t_c):
    return sum(c * t_c ** i for i, c in enumerate(_coef))


def fit(query):
    global _coef

    def q(a, b, c, d):
        return query({"value_a": a, "value_b": b, "value_c": c, "value_d": d})["output"]

    # Pairing-free, most activity-sensitive probe: z = 1, a = 8 A (> q at any T), 0.1 M.
    # A larger eps_r means a weaker correction, i.e. a larger (less negative) output.
    eps = []
    for t in TEMPS:
        y = q(8.0, 0.1, t, 1)
        eps.append(_invert(y, lambda e, t=t: _model(8.0, 0.1, t, 1, e), 40.0, 100.0))
    _coef = _lstsq(TEMPS, eps, 3)

    # Check the pairing hypothesis on a strongly associated probe (not used for fitting).
    y = q(4.0, 0.01, 25.0, 2)
    pred = _model(4.0, 0.01, 25.0, 2, _eps(25.0))
    if abs(y - pred) > 0.002:
        raise RuntimeError(f"Bjerrum pairing hypothesis rejected: {y} vs {pred:.4f}")


def solve(inputs):
    if _coef is None:
        raise RuntimeError("call fit(query) first")
    a, b, c = (float(inputs[k]) for k in ("value_a", "value_b", "value_c"))
    d = int(inputs["value_d"])
    return {"output": round(_model(a, b, c, d, _eps(c)), 4)}
