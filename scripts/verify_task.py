"""Automated GATE 2 checks. Exits non-zero on any failure."""
import ast
import json
import math
import random
import sys

from helpers import (GOLDEN_PATH, INPUT_KEYS, INPUT_RANGES, ORACLE_PATH, load_golden,
                     load_oracle, variant)

# Every numeric literal allowed in oracle/implement.py, with its justification.
ALLOWED_CONSTANTS = {
    1.602176634e-19: "elementary charge (exact, SI 2019)",
    1.380649e-23: "Boltzmann constant (exact, SI 2019)",
    6.02214076e23: "Avogadro constant (exact, SI 2019)",
    8.8541878188e-12: "vacuum permittivity (CODATA 2022)",
    273.15: "0 degC in K (exact)",
    1000.0: "L per m^3 (exact)",
    1e-10: "m per Angstrom (exact)",
    87.740: "Malmberg & Maryott 1956", 0.40008: "Malmberg & Maryott 1956",
    9.398e-4: "Malmberg & Maryott 1956", 1.410e-6: "Malmberg & Maryott 1956",
    2.0: "Bjerrum 1926: critical distance q = l/2 -> lower limit x = 2 of Q(b)",
    4: "structural: 4 pi (Coulomb / pair volume); x^-4 in Q(b); Simpson weight 4; output rounding",
    8: "Debye-Hueckel: 8 pi eps kT denominator (exact theory)",
    2: "structural: kappa^2 = 2 N_A e^2 I / eps kT; squares; Simpson weight 2",
    3: "structural: cubic term of eps_r(t); Simpson h/3",
    1: "structural: 1 + B a sqrt(I); 1 - alpha; alpha upper bound; z in (1, 2)",
    10: "log10 base",
    0.0: "structural: alpha lower bound; domain check",
    0: "structural: domain check",
    0.5: "bisection midpoint",
    2.5: "input-domain check (a lower bound, Angstrom)",
    8.0: "input-domain check (a upper bound, Angstrom)",
    0.0001: "input-domain check (c lower bound, mol/L)",
    0.1: "input-domain check (c upper bound, mol/L)",
    100: "input-domain check: upper temperature of Malmberg & Maryott fit (degC)",
    4000: "fixed Simpson panel count (determinism; quadrature error < 1e-9 relative)",
    200: "fixed bisection iteration count (determinism)",
}

failures = []


def check(name, cond, detail=""):
    print(f"[{'PASS' if cond else 'FAIL'}] {name}" + (f" -- {detail}" if detail else ""))
    if not cond:
        failures.append(name)


def numeric_literals(path):
    tree = ast.parse(open(path, encoding="utf-8").read())
    vals = []
    for node in ast.walk(tree):
        if isinstance(node, ast.UnaryOp) and isinstance(node.op, ast.USub) and isinstance(node.operand, ast.Constant):
            if isinstance(node.operand.value, (int, float)) and not isinstance(node.operand.value, bool):
                vals.append(-node.operand.value)
                node.operand.value = None
        elif isinstance(node, ast.Constant) and isinstance(node.value, (int, float)) and not isinstance(node.value, bool):
            vals.append(node.value)
    return vals


def sample(rng):
    lo, hi = (math.log10(x) for x in INPUT_RANGES["value_b"])
    return (rng.uniform(*INPUT_RANGES["value_a"]), 10 ** rng.uniform(lo, hi),
            rng.uniform(*INPUT_RANGES["value_c"]), rng.choice((1, 2)))


def main():
    o = load_oracle()
    rng = random.Random(7)
    pts = [sample(rng) for _ in range(300)]
    ref = {p: o.compute(*p) for p in pts}

    check("deterministic", all(o.compute(*p) == ref[p] for p in pts[:60]))
    out = o.oracle({"value_a": 4.0, "value_b": 0.01, "value_c": 25.0, "value_d": 2})
    check("output key is 'output' and scalar float",
          list(out) == ["output"] and isinstance(out["output"], float))
    check("4-decimal precision", all(round(v, 4) == v for v in ref.values()))

    unknown = [v for v in numeric_literals(ORACLE_PATH) if v not in ALLOWED_CONSTANTS]
    check("every oracle numeric literal is justified", not unknown, f"unjustified: {unknown}" if unknown else "")

    # Literature anchors
    eps25 = o.eps_r(25.0)
    check("eps_r(25 degC) reproduces Malmberg & Maryott 78.30", abs(eps25 - 78.30) < 0.01, f"{eps25:.3f}")
    length, a25, b25 = o.electrostatics(25.0, 1)
    check("A(25 degC) ~ 0.51 (mol/L)^-1/2", abs(a25 - 0.5115) < 1e-3, f"{a25:.5f}")
    check("B(25 degC) ~ 3.29 nm^-1 (mol/L)^-1/2", abs(b25 * 1e-9 - 3.29) < 0.01, f"{b25 * 1e-9:.4f}")
    check("Bjerrum q(1:1, 25 degC) ~ 3.57 A", abs(length / 2 * 1e10 - 3.57) < 0.02, f"{length / 2 * 1e10:.4f}")
    q_exact = o.bjerrum_q(10.0)
    # closed form cross-check: Q(b) = [Ei(b) - Ei(2)]/6 - e^x (1/(3x^3) + 1/(6x^2) + 1/(6x)) |_2^b
    ei = lambda x: 0.5772156649015329 + math.log(x) + sum(x ** k / (k * math.factorial(k)) for k in range(1, 120))
    prim = lambda x: ei(x) / 6 - math.exp(x) * (1 / (3 * x ** 3) + 1 / (6 * x ** 2) + 1 / (6 * x))
    check("Simpson Q(10) matches exponential-integral closed form", abs(q_exact / (prim(10) - prim(2)) - 1) < 1e-9,
          f"{q_exact:.8f} vs {prim(10) - prim(2):.8f}")

    check("helper variant(oracle settings) == oracle", all(variant(*p) == ref[p] for p in pts))

    # Twist integrity: Bjerrum pairing with a critical-distance cut-off
    t_none = max(abs(ref[p] - variant(*p, pairing="none")) for p in pts)
    t_fuoss = max(abs(ref[p] - variant(*p, pairing="fuoss")) for p in pts)
    check("pairing twist changes output (> 100x tolerance)", t_none > 0.1, f"max diff vs no pairing {t_none:.4f}")
    check("cut-off distinguishes Bjerrum from Fuoss (> 10x tolerance)", t_fuoss > 0.01, f"max diff {t_fuoss:.4f}")
    above = [p for p in pts if p[3] == 1 and p[0] >= o.electrostatics(p[2], 1)[0] / 2 * 1e10]
    check("a >= q  =>  pairing exactly neutral (bit-exact vs no-pairing variant)",
          len(above) > 50 and all(ref[p] == variant(*p, pairing="none") for p in above), f"{len(above)} points")
    on = o.compute(3.8, 0.1, 100.0, 1) - variant(3.8, 0.1, 100.0, 1, pairing="none")
    off = o.compute(3.8, 0.1, 0.0, 1) - variant(3.8, 0.1, 0.0, 1, pairing="none")
    check("temperature alone switches pairing on (a = 3.8 A, 1:1)", off == 0 and on < -0.002,
          f"pairing shift {off:.4f} at 0 degC vs {on:.4f} at 100 degC")
    st = max(abs(ref[p] - variant(*p, stoichiometric=False)) for p in pts)
    check("stoichiometric alpha factor matters", st > 0.1, f"max diff {st:.4f}")
    th = max(abs(ref[p] - variant(*p, thermal=False)) for p in pts)
    check("eps_r(T) law matters", th > 0.05, f"max diff {th:.4f}")
    dh = max(abs(ref[p] - variant(*p, ion_size=False)) for p in pts)
    check("ion-size term matters", dh > 0.05, f"max diff {dh:.4f}")
    check("output monotone decreasing in concentration (spot)",
          all(o.compute(p[0], 0.1, p[2], p[3]) < o.compute(p[0], 0.001, p[2], p[3]) for p in pts[:50]))

    g = load_golden()
    cats = [c["category"] for c in g["cases"]]
    check(">=2 discriminating edge cases", cats.count("discriminating") >= 2)
    check(">=2 control edge cases", cats.count("control") >= 2)
    check(">=1 boundary edge case", cats.count("boundary") >= 1)
    check("golden inputs neutral and in range", all(
        sorted(c["input"]) == sorted(INPUT_KEYS)
        and all(INPUT_RANGES[k][0] <= c["input"][k] <= INPUT_RANGES[k][1] for k in INPUT_KEYS)
        and isinstance(c["input"]["value_d"], int)
        for c in g["cases"]))
    check("golden expected outputs match oracle",
          all(o.oracle(c["input"]) == c["expected"] for c in g["cases"]), GOLDEN_PATH)

    print(json.dumps({"failures": failures}))
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
