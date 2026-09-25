# Blackbox Task 49 - Stoichiometric Activity Coefficient of a z:z Salt with Bjerrum Pairing Switched On at the Critical Distance

**Domain:** Chemistry - physical chemistry / electrolyte solution theory (Debye-Hückel theory, ion association, dielectric properties of water)
**Tier:** Low (floor: 2h, 1 genuine twist, 3 substantive steps)
**Honest counts:** 4 substantive steps (plus a trivial rounding step), 1 genuine twist

**Repo name:** `electrolyte-bjerrum-pairing-critical-distance-activity`
**GitHub description:** Reverse-engineer the measured mean ionic activity coefficient of 1:1 and 2:2 salts in water from 0 to 100 °C, where free ions follow extended Debye-Hückel theory and ion pairs form only inside Bjerrum's critical distance. Heating can switch pairing on. Requires electrolyte theory, the Bjerrum association integral, and the dielectric behaviour of water.

---

## Task description (for reviewers)

The hidden function maps three real inputs and one small integer to the log of a thermodynamic property of an aqueous electrolyte. The free-ion part is a textbook activity model. On top of it sits an association equilibrium whose constant is an integral with a hard lower cut-off, so it is exactly zero over part of the domain and turns on with a kink. The reported quantity is on the stoichiometric scale, so it combines the fraction of free ions with their activity. Expected solver knowledge: graduate electrolyte theory (Robinson & Stokes level), including the difference between Bjerrum and Fuoss association and why water's permittivity makes the Coulomb coupling grow as it heats.

---

## 1. Domain & Algorithm

Electrolyte solution chemistry. The output is **log₁₀ γ±** (stoichiometric) of a symmetric z:z salt (z = 1 or 2) at molar concentration c in water at t °C. Cation and anion have distance of closest approach a. The ions are split into free ions and Bjerrum ion pairs:

- Free ions follow the extended Debye-Hückel law with ion-size parameter a. Their ionic strength is I = αcz².
- Pairs form according to Bjerrum's association constant K_A(a, z, T), which is **zero whenever a ≥ q**, where q is Bjerrum's critical distance.
- α comes from mass action with activity-corrected free ions, solved exactly by bisection.
- All temperature dependence enters through the relative permittivity ε_r(t) of water (Malmberg-Maryott) and kT, using exact SI constants and CODATA 2022 ε₀.

## 2. Core Method

- **Relative permittivity of water.** Malmberg & Maryott (1956), *J. Res. NBS* 56, 1: ε_r = 87.740 − 0.40008 t + 9.398×10⁻⁴ t² − 1.410×10⁻⁶ t³ (t in °C, 0-100 °C).
- **Pair length and Bjerrum parameter.** l = z²e²/(4πε₀ε_r kT), b = l/a. The critical distance is q = l/2 (Bjerrum 1926, the minimum of the pair-distribution integrand).
- **Bjerrum association constant.** K_A = 4πN_A l³ Q(b), with Q(b) = ∫₂ᵇ x⁻⁴ eˣ dx for b > 2, and K_A = 0 for b ≤ 2 (a ≥ q). Here ×1000 converts m³ mol⁻¹ to L mol⁻¹. Source: Bjerrum 1926; Robinson & Stokes ch. 14. The oracle integrates Q(b) with 4000-panel composite Simpson. `verify_task.py` checks this against the exponential-integral closed form to 1e-9.
- **Extended Debye-Hückel for free ions.** log₁₀ γ_f = −A z² √I / (1 + B a √I), with A = e²B/(ln10 · 8πε₀ε_r kT) and B = (2N_A e² · 1000/(ε₀ε_r kT))^½ (Debye & Hückel 1923; Robinson & Stokes ch. 9).
- **Mass action.** K_A = (1 − α)/(α² c γ_f²). The pair is neutral, so γ_IP = 1.
- **Stoichiometric scale.** μ_salt = μ° + 2RT ln(c γ±) = μ° + 2RT ln(αc γ_f), which gives γ± = α γ_f. This is exact thermodynamics and has no free constant.

## 3. Twist (1 genuine)

**T1 - Ion pairing that exists only inside Bjerrum's critical distance.**

The textbook activity model for a strong electrolyte is extended Debye-Hückel. When association is added, the formula most often quoted in practice is Fuoss (1958), K_A = (4πN_A a³/3) e^b. It is smooth in a and never vanishes. Bjerrum's constant instead integrates the Boltzmann-weighted pair distribution only from a out to q = l/2. That gives three behaviours that a DH or Fuoss hypothesis cannot reproduce:

1. **Exact neutrality for a ≥ q.** For 1:1 salts, q is 3.49 Å at 0 °C and 4.02 Å at 100 °C. Above q the output equals the plain extended-DH value bit for bit (107/107 sampled points). Below q, pairing appears with a kink in ∂output/∂a. The integrand at x = 2 is e²/16, so the kink is in the slope, not the value.
2. **Heating switches pairing on.** ε_r T falls with temperature (78.3 × 298 K > 55.7 × 373 K), so l and q grow on heating. At a = 3.8 Å, 0.1 M 1:1, there is no pairing at 0 °C but pairing at 100 °C (K_A = 0.20 L mol⁻¹, −0.0041 extra). For 2:2 at a = 4 Å, K_A rises from 265 to 509 L mol⁻¹ between 25 and 100 °C. The intuition "ions dissociate more when hot" is wrong in water.
3. **Magnitude set by l³ Q(b), not a³ e^b.** For 2:2 salts, q ≈ 14 Å exceeds every a in the domain, so pairing is always on. Its size distinguishes Bjerrum from Fuoss by up to 0.20 in log γ±.

The stoichiometric factor α in γ± = αγ_f is required to fit the data, but it is standard thermodynamics, so it is **not** counted as a twist. The same goes for ε_r(T), which is pipeline rigour here.

I/O isolation:
- value_d = 1 with value_a ≥ 4.02 Å makes T1 exactly neutral at every temperature (control cases).
- value_a and value_c each move the output across the threshold independently (a sweep, or a temperature sweep at a = 3.8 Å).
- value_d = 2 puts T1 on everywhere, with its largest effect.

## 4. Multi-Step Pipeline

| Step | Computation | Inputs involved |
|---|---|---|
| S1 | ε_r(t) → l, b, A(T), B(T) | value_a, value_c, value_d |
| S2 | Bjerrum K_A = 4πN_A l³ Q(b) with the critical-distance cut-off (T1) | value_a, value_c, value_d |
| S3 | Solve K_A = (1 − α)/(α² c γ_f²) with γ_f from extended DH at I = αcz² (bisection; the root is unique, checked on a 4000-point log grid over 1500 random points) | all |
| S4 | γ± = α γ_f → log₁₀ | all |
| (S5) | round to 4 decimals | trivial, not counted |

Each step is substantive. S2 needs a special-function integral. S3 is a genuinely coupled non-linear solve: α sets I, I sets γ_f, and γ_f sets α. Dropping the coupled activity term moves 2:2 outputs by tenths.

## 5. Input Schema

```json
{
  "value_a": "float [2.5, 8.0]",     // distance of closest approach a, Angstrom
  "value_b": "float [0.0001, 0.1]",  // stoichiometric salt concentration, mol L^-1
  "value_c": "float [0, 100]",       // temperature, degC (validity range of Malmberg & Maryott)
  "value_d": "int {1, 2}"            // charge number z of the symmetric z:z salt
}
```

## 6. Output Schema

```json
{ "output": "float, 4 decimals" }    // log10 of the stoichiometric mean ionic activity coefficient; range ~ -1.40 to -0.005
```

Grading tolerance is 0.001 (absolute, log₁₀ units).

## 7. Hardness & Twists

These baselines are measured on the 52 golden cases (`.scratch/baselines.py`, tolerance 0.001):

| Model | Pass | Median abs err | Max abs err |
|---|---|---|---|
| Extended DH only (no pairing) | 24/52 | 0.0024 | 0.408 |
| Fuoss 1958 K_A (no cut-off) | 11/52 | 0.0128 | 0.201 |
| Bjerrum, but free-ion γ reported (α forgotten) | 24/52 | 0.0025 | 0.977 |
| Electrostatics frozen at 25 °C | 21/52 | 0.0027 | 0.232 |
| DH limiting law instead of extended | 9/52 | 0.0195 | 0.339 |
| Full model (example solver, 10 probes) | 52/52 | 6e-5 | - |

Why this is hard:
- **Plain DH fits half the domain perfectly.** Any 1:1 probe with a ≳ 4 Å is exactly extended DH, which invites the solver to stop there. The discrepancy appears only for small a, for 2:2 salts, or at high temperature.
- **The kink is the fingerprint.** A smooth association law (Fuoss, or any fitted exponential) fails both above q, where it predicts pairing that isn't there, and below it, where its magnitude is wrong. The solver must notice that the onset moves with temperature exactly as ε_r T does.
- **Counter-intuitive temperature direction.** Pairing increasing on heating is a known property of water, but only to someone who reasons through ε_r(T).
- **Linear or additive surrogates fail.** The output is non-separable: a enters both the DH denominator and b, and temperature enters A, B, and l. The integer z raises l by 4× and the DH term by z².
- **Cold read:** four neutral inputs, one of them an integer in {1, 2}, and a negative output near 0. "Some log correction" is guessable. The Bjerrum cut-off, the association law, and the stoichiometric α are not.

## 8. Constants & Citations

Every numeric literal in `oracle/implement.py` is whitelisted with its justification in `scripts/verify_task.py` (AST scan).

| Constant | Value | Status | Source |
|---|---|---|---|
| e | 1.602176634×10⁻¹⁹ C | EXACT | SI 2019 defining constant |
| k_B | 1.380649×10⁻²³ J K⁻¹ | EXACT | SI 2019 defining constant |
| N_A | 6.02214076×10²³ mol⁻¹ | EXACT | SI 2019 defining constant |
| ε₀ | 8.8541878188×10⁻¹² F m⁻¹ | CITED | CODATA 2022 (NIST) |
| 273.15, 1000, 1e-10 | - | EXACT | °C → K; L per m³; m per Å |
| 87.740, 0.40008, 9.398e-4, 1.410e-6 | - | CITED | Malmberg & Maryott 1956 |
| 2.0 (lower limit of Q) | - | CITED | Bjerrum 1926: q = l/2 ⇒ x = l/r = 2 |
| 4, 8, 2 (4π, 8π, κ² factor 2) | - | EXACT | Coulomb law, Debye-Hückel theory |
| 4000, 200 | - | STRUCTURAL | fixed Simpson panels / bisection steps (determinism) |
| 2.5, 8.0, 0.0001, 0.1, 100 | - | STRUCTURAL | input-domain check (100 °C = top of Malmberg-Maryott fit) |

### Citations fetched & verified (quote < 15 words, location)

| Item | Verbatim quote | Where fetched |
|---|---|---|
| ε_r(t) | "E = 87.740 - 0.40008t + 9.398(10^-4)t^2 - 1.410(10^-6)t^3" | Malmberg & Maryott 1956, via archive.org jresv56n1p1 / academia.edu listing |
| ε_r(25 °C) | "a value of 78.30 at 25° C" | same (oracle gives 78.303) |
| q | "q = \|z_i z_j\| e² / (8π ε₀ ε_r k T)" | LibreTexts, *Topics in Thermodynamics of Solutions*, 1.16.1 Ion Association |
| b | "b = \|z₊ z₋\| e² / (4π ε₀ ε_r k T a)" | same |
| K_A form | "(4π N / 10³) × (\|z₊ z₋\| e² / 4π ε₀ ε_r k T)³ × Q(b)" | same |
| Q(b) limits | "∫ e^y y^−4 dy" from 2 to b (eq. 216) | M. Dalal, *Textbook of Physical Chemistry* Vol. 1, §4.12, eq. 216 |
| Cut-off | "formation of ion-pair is feasible if a < q and infeasible if a > q" | same, after eq. 214 |
| Extended DH | "log₁₀γ± = −Az_j² √I/(1 + Ba₀√I)" | LibreTexts *Physical Chemistry* 25.6 (A = 0.51, B = 3.29 nm⁻¹ at 25 °C) |
| ε₀ | "8.854 187 8188 x 10⁻¹² F m⁻¹" | physics.nist.gov, CODATA 2022 |

**Citation-code match:** `MM_0..MM_3 = 87.740, 0.40008, 9.398e-4, 1.410e-6`; `BJERRUM_X_MIN = 2.0`; `association_constant = 4 π N_A l³ Q(l/a) × 1000`; `electrostatics()` gives A = 0.51156 and B = 3.2914 nm⁻¹ at 25 °C, and q(1:1, 25 °C) = 3.579 Å. `verify_task.py` asserts all four anchors plus the Q(b) closed form.

**Honest notes (CONCERN-level, documented):**
- **Primary texts not fetched.** Bjerrum (1926) and Robinson & Stokes (ch. 14) are not open-access. The formulas were verified in two independent secondary sources (LibreTexts, Dalal) that reproduce Bjerrum's equations exactly.
- **Free-ion size parameter.** The oracle uses the same a in the DH denominator as in b, which is the simplest self-consistent Bjerrum-DH model. Some treatments (following Bjerrum's own argument) use q for free ions when a < q. That is a documented modelling choice, not a hidden constant.
- **Validity.** Extended DH plus Bjerrum is a model, not experimental truth. At I up to 0.4 (2:2, 0.1 M) it is outside DH's quantitative comfort zone, but it stays deterministic and literature-defined.
- **Molar vs molal.** Concentrations are mol L⁻¹ throughout. Density corrections are omitted.
- **3:3 excluded.** For z = 3 the mass-action residual is non-monotone in α in parts of the domain (multiple roots possible), so value_d ∈ {1, 2}. Uniqueness for z ∈ {1, 2} was verified on 1500 random points.

## 9. Edge Cases

| id | value_a | value_b | value_c | value_d | output | Rationale |
|---|---|---|---|---|---|---|
| edge_control_1 | 5.0 | 0.0001 | 25 | 1 | −0.0050 | a > q; ~DH limiting law |
| edge_control_2 | 6.0 | 0.01 | 25 | 1 | −0.0427 | a > q; pure extended DH |
| edge_control_3 | 4.5 | 0.1 | 25 | 1 | −0.1102 | a > q; pairing exactly 0 (Fuoss: −0.1350) |
| edge_discrim_1 | 2.5 | 0.1 | 25 | 1 | −0.1446 | a < q: pairing lowers output by 0.016 (no-pair −0.1284) |
| edge_discrim_2 | 3.8 | 0.1 | 0 | 1 | −0.1119 | q = 3.49 Å < a: no pairs |
| edge_discrim_3 | 3.8 | 0.1 | 100 | 1 | −0.1397 | q = 4.02 Å > a: heating switches pairing on |
| edge_discrim_4 | 4.0 | 0.01 | 25 | 2 | −0.4507 | 2:2, K_A = 265, α = 0.67 (no-pair −0.3239) |
| edge_discrim_5 | 4.0 | 0.01 | 100 | 2 | −0.5420 | K_A rises to 509 on heating, α = 0.58 |
| edge_boundary_1 | 3.6 | 0.1 | 25 | 1 | −0.1177 | just above q = 3.579 Å: K_A = 0 |
| edge_boundary_2 | 3.5 | 0.1 | 25 | 1 | −0.1198 | just below q: pairing onset (kink) |
| edge_boundary_3 | 2.5 | 0.1 | 100 | 2 | −1.4009 | most associated corner, K_A ≈ 5700, α ≈ 0.11 |
| edge_boundary_4 | 8.0 | 0.0001 | 0 | 1 | −0.0048 | infinite-dilution limit → 0 |

The golden file adds 40 seeded random cases (seed 49), for 52 cases in total.

---

## GATE 1: Screening

| Test | Result | Notes |
|---|---|---|
| Algebraic collapse | PASS | a appears in both the DH denominator and b. α is the root of a transcendental equation coupled to γ_f. K_A has a non-analytic switch. No separable additive form exists. |
| Domain recall | PASS (CONCERN) | An expert may guess "activity coefficient with ion pairing". The specific Bjerrum cut-off (vs Fuoss), the shared a, the stoichiometric scale, and the ε_r(t) law are not named by the schema. |
| Twist survives | PASS | The cut-off is a hard zero of K_A for a ≥ q. It cannot be simplified away and is bit-exact detectable. |
| Genuine twist | PASS (1) | Bjerrum's critical-distance association departs from both the no-pairing textbook model and the widely used Fuoss form. α and ε_r(T) are rigour and are not counted. |
| Two-trap | PASS | value_a and value_c independently move the threshold. value_d switches pairing to always-on. |
| Tier fit | PASS | 4 steps ≥ 3; 1 twist ≥ 1 (Low). |
| PhD authenticity | PASS | Bjerrum vs Fuoss association, and pairing increasing with T in water, are graduate-level electrolyte topics (Robinson & Stokes ch. 14; Marcus & Hefter 2006 review). |

**Verdict: PASS** (one documented CONCERN on domain recall)

## GATE 2: Verification

| Check | Result |
|---|---|
| Twist integrity (freeze-one tests) | PASS: over 300 random points, pairing effect max 0.42 vs none and 0.16 vs Fuoss (`verify_task.py`) |
| I/O isolation | PASS: a ≥ q ⇒ bit-exact neutral (107 points); at a = 3.8 Å, 1:1, temperature alone moves pairing from 0 to −0.0041 |
| No invented constant | PASS: AST scan finds only whitelisted literals |
| Recall-reconstructability | PASS: the example solver recovers the function with 10 probes (≪ 50) once the family is hypothesised |
| Citations fetched & verified | PASS: five sources fetched (table above); primary Bjerrum text via two secondary reproductions (CONCERN noted) |
| Citation-code match | PASS: ε_r(25) = 78.30, A = 0.5116, B = 3.29 nm⁻¹, q = 3.58 Å, Q(b) closed form, all asserted |
| Output/schema neutralisation | PASS: key "output", float, value_a to value_d |
| Tier floor | PASS: 4 steps / 1 twist vs Low 3 / 1 |
| Edge cases ≥ 5 | PASS: 12 total (3 control, 5 discriminating, 4 boundary) |
| PhD-level QA | PASS: needs DH theory, the Bjerrum integral with its cut-off, and dielectric reasoning together |

**Verdict: PASS**

---

## Files

```
task-49/
├── oracle/implement.py          reference oracle (hidden)
├── solution/example_solver.py   probe + fit solver (10 probes, 52/52)
├── grader/grade_submission.py   tolerance 0.001 log10 units absolute, probe budget 200
├── golden/test_data.json        52 cases (12 edge + 40 random)
├── scripts/helpers.py           shared utilities, twist-freezing variant()
├── scripts/build_golden.py      regenerates golden data
├── scripts/verify_task.py       automated GATE 2 checks
├── scripts/run_tests.py         verify + oracle self-grade + solver grade
└── PROPOSAL.md
```

## References

- Bjerrum, N. (1926). Untersuchungen über Ionenassoziation. I. *Kgl. Danske Vidensk. Selsk. Math.-fys. Medd.* 7(9), 1-48.
- Debye, P. & Hückel, E. (1923). Zur Theorie der Elektrolyte. *Phys. Z.* 24, 185-206.
- Fuoss, R. M. (1958). Ionic association. III. The equilibrium between ion pairs and free ions. *J. Am. Chem. Soc.* 80, 5059-5061. (baseline only)
- Malmberg, C. G. & Maryott, A. A. (1956). Dielectric constant of water from 0° to 100 °C. *J. Res. Natl. Bur. Stand.* 56, 1-8. doi:10.6028/jres.056.001
- Robinson, R. A. & Stokes, R. H. (1959). *Electrolyte Solutions*, 2nd ed. Butterworths (ch. 9, 14).
- Marcus, Y. & Hefter, G. (2006). Ion pairing. *Chem. Rev.* 106, 4585-4621. doi:10.1021/cr040087x
- CODATA 2022 recommended values, NIST. https://physics.nist.gov/cuu/Constants/
- Open-access verification sources: LibreTexts "1.16.1 Ion Association" and "25.6 The Debye-Hückel Theory"; M. Dalal, *A Textbook of Physical Chemistry* Vol. 1, §4.12 (dalalinstitute.com preview PDF); archive.org jresv56n1p1.
