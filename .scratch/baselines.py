import sys, math, statistics
sys.path.insert(0, '.')
from helpers import *
o = load_oracle(); cases = load_golden()["cases"]
models = {
 "Extended DH only (no pairing)": dict(pairing="none"),
 "Fuoss 1958 K_A (no cut-off)": dict(pairing="fuoss"),
 "Bjerrum, free-ion gamma reported (no alpha)": dict(stoichiometric=False),
 "Electrostatics frozen at 25 degC": dict(thermal=False),
 "DH limiting law in place of extended": dict(ion_size=False),
}
for name, kw in models.items():
    errs = [abs(variant(*(c["input"][k] for k in INPUT_KEYS), **kw) - c["expected"]["output"]) for c in cases]
    print(f"| {name} | {sum(e <= TOLERANCE for e in errs)}/{len(cases)} | {statistics.median(errs):.4f} | {max(errs):.4f} |")
outs=[c["expected"]["output"] for c in cases]; print("range", min(outs), max(outs))
# alpha for labels
def alpha(a,c,t,z):
    l,A,B=o.electrostatics(t,z); am=a*1e-10; ka=o.association_constant(am,l)
    return ka, l/am, (10**(o.compute(a,c,t,z)))/(10**variant(a,c,t,z,stoichiometric=False))
for c in [(2.5,0.1,25,1),(3.8,0.1,100,1),(4.0,0.01,25,2),(4.0,0.01,100,2),(2.5,0.1,100,2),(3.5,0.1,25,1)]:
    print(c, ["%.4g"%x for x in alpha(*c)], "none", variant(*c,pairing='none'), "fuoss", variant(*c,pairing='fuoss'))
