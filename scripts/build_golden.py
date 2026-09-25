"""Regenerate golden/test_data.json from the oracle (deterministic, seeded)."""
import json
import math
import random

from helpers import GOLDEN_PATH, load_oracle

# (id, category, label, value_a, value_b, value_c, value_d)
EDGE_CASES = [
    ("edge_control_1", "control", "1:1, large a, 1e-4 M, 25 degC: pairing off (a > q), output ~ DH limiting law",
     5.0, 0.0001, 25.0, 1),
    ("edge_control_2", "control", "1:1, a = 6 A, 0.01 M, 25 degC: a > q, pure extended Debye-Hueckel",
     6.0, 0.01, 25.0, 1),
    ("edge_control_3", "control", "1:1, a = 4.5 A, 0.1 M, 25 degC: a > q, pairing exactly zero (Fuoss-type K_A would give -0.135)",
     4.5, 0.1, 25.0, 1),
    ("edge_discrim_1", "discriminating", "1:1, a = 2.5 A < q = 3.58 A, 0.1 M, 25 degC: Bjerrum pairs lower log gamma by 0.016 vs no pairing",
     2.5, 0.1, 25.0, 1),
    ("edge_discrim_2", "discriminating", "1:1, a = 3.8 A, 0.1 M, 0 degC: q = 3.49 A < a, no pairs",
     3.8, 0.1, 0.0, 1),
    ("edge_discrim_3", "discriminating", "Same salt at 100 degC: q = 4.02 A > a, heating switches pairing ON (eps_r T falls)",
     3.8, 0.1, 100.0, 1),
    ("edge_discrim_4", "discriminating", "2:2, a = 4 A, 0.01 M, 25 degC: b = 7.2, K_A ~ 265 L/mol, alpha = 0.67",
     4.0, 0.01, 25.0, 2),
    ("edge_discrim_5", "discriminating", "2:2, a = 4 A, 0.01 M, 100 degC: pairing grows with temperature (K_A 265 -> 509 L/mol, alpha = 0.58)",
     4.0, 0.01, 100.0, 2),
    ("edge_boundary_1", "boundary", "1:1 just above the critical distance (a = 3.60 A > q = 3.579 A): K_A = 0 exactly",
     3.6, 0.1, 25.0, 1),
    ("edge_boundary_2", "boundary", "1:1 just below the critical distance (a = 3.50 A): pairing switches on with a kink",
     3.5, 0.1, 25.0, 1),
    ("edge_boundary_3", "boundary", "Most associated corner: 2:2, a = 2.5 A, 0.1 M, 100 degC (b = 12.9, K_A ~ 5700 L/mol, alpha ~ 0.11)",
     2.5, 0.1, 100.0, 2),
    ("edge_boundary_4", "boundary", "Infinite-dilution corner: 1:1, 1e-4 M, 0 degC, a = 8 A: output -> 0",
     8.0, 0.0001, 0.0, 1),
]

N_RANDOM = 40
SEED = 49


def main():
    o = load_oracle()
    cases = []
    for cid, cat, label, a, b, c, d in EDGE_CASES:
        inp = {"value_a": a, "value_b": b, "value_c": c, "value_d": d}
        cases.append({"id": cid, "category": cat, "label": label, "input": inp, "expected": o.oracle(inp)})
    rng = random.Random(SEED)
    lo, hi = math.log10(0.0001), math.log10(0.1)
    for i in range(N_RANDOM):
        inp = {
            "value_a": round(rng.uniform(2.5, 8.0), 3),
            "value_b": round(10 ** rng.uniform(lo, hi), 6),
            "value_c": round(rng.uniform(0.0, 100.0), 2),
            "value_d": rng.choice((1, 2)),
        }
        cases.append({"id": f"rand_{i:02d}", "category": "random", "label": "uniform value_a, value_c; log-uniform value_b",
                      "input": inp, "expected": o.oracle(inp)})
    doc = {
        "schema": {
            "input": {
                "value_a": "float [2.5, 8.0]",
                "value_b": "float [0.0001, 0.1]",
                "value_c": "float [0, 100]",
                "value_d": "int {1, 2}",
            },
            "output": {"output": "float, 4 decimals"},
        },
        "seed": SEED,
        "cases": cases,
    }
    assert all(math.isfinite(c["expected"]["output"]) for c in cases)
    with open(GOLDEN_PATH, "w", encoding="utf-8") as fh:
        json.dump(doc, fh, indent=2)
        fh.write("\n")
    print(f"wrote {len(cases)} cases to {GOLDEN_PATH}")


if __name__ == "__main__":
    main()
