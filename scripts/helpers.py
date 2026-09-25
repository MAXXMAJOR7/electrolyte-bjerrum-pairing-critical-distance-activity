"""Shared utilities for task scripts and grader."""
import importlib.util
import json
import math
import os

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
GOLDEN_PATH = os.path.join(ROOT, "golden", "test_data.json")
ORACLE_PATH = os.path.join(ROOT, "oracle", "implement.py")

INPUT_KEYS = ("value_a", "value_b", "value_c", "value_d")
INPUT_RANGES = {
    "value_a": (2.5, 8.0),
    "value_b": (0.0001, 0.1),
    "value_c": (0.0, 100.0),
    "value_d": (1, 2),
}
TOLERANCE = 0.001  # absolute, log10 units (output spans ~ -1.2 .. -0.005)


def load_module(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def load_oracle():
    return load_module(ORACLE_PATH, "oracle_impl")


def load_golden(path=GOLDEN_PATH):
    with open(path, "r", encoding="utf-8") as fh:
        return json.load(fh)


def variant(a_ang, conc, t_c, z, pairing="bjerrum", stoichiometric=True, thermal=True, ion_size=True):
    """Oracle with individual ingredients swapped out, for twist-integrity checks and baselines.

    pairing="bjerrum"   -> oracle: K_A = 4 pi N_A l^3 Q(b), zero for b <= 2 (critical distance)
    pairing="none"      -> plain extended Debye-Hueckel, no ion pairs (K_A = 0)
    pairing="fuoss"     -> Fuoss (1958) K_A = 4 pi N_A a^3 e^b / 3: no critical-distance cut-off
    stoichiometric=False -> reports the free-ion gamma (forgets the alpha factor)
    thermal=False       -> all electrostatics (eps_r and T inside l, A, B) frozen at 25 degC
    ion_size=False      -> Debye-Hueckel limiting law (no 1 + B a sqrt(I) denominator)
    """
    o = load_oracle()
    z = int(z)
    a_m = a_ang * o.M_PER_ANGSTROM
    t_eval = t_c if thermal else 25.0
    length, a_dh, b_dh = o.electrostatics(t_eval, z)
    if pairing == "bjerrum":
        k_a = o.association_constant(a_m, length)
    elif pairing == "fuoss":
        k_a = 4 * math.pi * o.N_A * a_m ** 3 * math.exp(length / a_m) / 3 * o.L_PER_M3
    else:
        k_a = 0.0

    def lg(alpha):
        s = math.sqrt(alpha * conc * z * z)
        return -a_dh * z * z * s / (1 + (b_dh * a_m * s if ion_size else 0.0))

    lo, hi = 0.0, 1.0
    for _ in range(o.N_BISECT):
        mid = 0.5 * (lo + hi)
        g = 10 ** lg(mid)
        if k_a * mid * mid * conc * g * g - (1 - mid) > 0:
            hi = mid
        else:
            lo = mid
    alpha = 0.5 * (lo + hi)
    return round((math.log10(alpha) if stoichiometric else 0.0) + lg(alpha), 4)
