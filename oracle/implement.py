"""Reference oracle (hidden from solver).

Stoichiometric mean ionic activity coefficient, log10 gamma_pm, of a symmetric
z:z electrolyte in water, treated as free ions (extended Debye-Hueckel) in
equilibrium with Bjerrum ion pairs.

Pipeline (see PROPOSAL.md sections 2-4):
  S1  solvent + electrostatics  eps_r(t) = 87.740 - 0.40008 t + 9.398e-4 t^2 - 1.410e-6 t^3
                                                                  (Malmberg & Maryott 1956)
                                pair length  l = z^2 e^2 / (4 pi eps0 eps_r k T),  b = l / a
                                DH constants A(T), B(T) from CODATA constants and eps_r(t)
  S2  Bjerrum association       K_A = 4 pi N_A l^3 Q(b),  Q(b) = int_2^b x^-4 e^x dx   (b > 2)
                                K_A = 0 when b <= 2 (a >= q = l/2: no pairs can form)   (twist)
  S3  mass action + activity    K_A = (1 - alpha) / (alpha^2 c gamma_f^2),
                                log10 gamma_f = -A z^2 sqrt(I) / (1 + B a sqrt(I)),  I = alpha c z^2
  S4  stoichiometric scale      gamma_pm = alpha * gamma_f   (salt chemical potential = free-ion one)

Inputs (neutral names):
    value_a : distance of closest approach a, Angstrom     [2.5, 8.0]
    value_b : stoichiometric salt concentration, mol L^-1  [0.0001, 0.1]
    value_c : temperature, degC                            [0, 100]
    value_d : ionic charge number z of the z:z salt (int)  {1, 2}
Output:
    output  : log10 gamma_pm (stoichiometric), rounded to 4 decimals
"""
import math

# --- SI exact defining constants and CODATA 2022 --------------------------------
E_CHARGE = 1.602176634e-19      # C, elementary charge (exact, SI 2019)
K_B = 1.380649e-23              # J K^-1, Boltzmann constant (exact, SI 2019)
N_A = 6.02214076e23             # mol^-1, Avogadro constant (exact, SI 2019)
EPS_0 = 8.8541878188e-12        # F m^-1, vacuum permittivity (CODATA 2022)
T_ZERO = 273.15                 # K, 0 degC (exact)
L_PER_M3 = 1000.0               # L m^-3 (exact)
M_PER_ANGSTROM = 1e-10          # m per Angstrom (exact)

# --- Malmberg & Maryott (1956), J. Res. NBS 56, 1: eps_r(t), 0-100 degC ----------
MM_0, MM_1, MM_2, MM_3 = 87.740, 0.40008, 9.398e-4, 1.410e-6

# --- Bjerrum (1926): pairing cut-off at q = l/2, i.e. lower limit x = 2 of Q(b) --
BJERRUM_X_MIN = 2.0

N_SIMPSON = 4000                # fixed quadrature panels (even) -> bit-for-bit determinism
N_BISECT = 200                  # fixed bisection count on alpha in (0, 1]


def eps_r(t_c):
    return MM_0 - MM_1 * t_c + MM_2 * t_c ** 2 - MM_3 * t_c ** 3


def bjerrum_q(b):
    """Q(b) = integral_2^b x^-4 e^x dx (composite Simpson); zero when b <= 2."""
    if b <= BJERRUM_X_MIN:
        return 0.0
    h = (b - BJERRUM_X_MIN) / N_SIMPSON

    def f(x):
        return math.exp(x) / x ** 4

    s = f(BJERRUM_X_MIN) + f(b)
    for i in range(1, N_SIMPSON):
        s += (4 if i % 2 else 2) * f(BJERRUM_X_MIN + i * h)
    return s * h / 3


def electrostatics(t_c, z):
    """Return (l [m], A [log10, (mol L^-1)^-1/2], B [m^-1 (mol L^-1)^-1/2]) at t_c degC."""
    temp = t_c + T_ZERO
    ekt = EPS_0 * eps_r(t_c) * K_B * temp
    length = z * z * E_CHARGE ** 2 / (4 * math.pi * ekt)
    # kappa^2 = 2 N_A e^2 I / (eps kT) with I in mol m^-3 -> B = kappa / sqrt(I [mol L^-1])
    b_dh = math.sqrt(2 * N_A * L_PER_M3 * E_CHARGE ** 2 / ekt)
    # ln gamma = -z^2 e^2 kappa / (8 pi eps kT (1 + kappa a))
    a_dh = E_CHARGE ** 2 * b_dh / (8 * math.pi * ekt) / math.log(10)
    return length, a_dh, b_dh


def association_constant(a_m, length):
    """Bjerrum K_A in L mol^-1 for contact distance a_m [m] and pair length l [m]."""
    return 4 * math.pi * N_A * length ** 3 * bjerrum_q(length / a_m) * L_PER_M3


def compute(value_a, value_b, value_c, value_d):
    a_ang = float(value_a)
    conc = float(value_b)
    t_c = float(value_c)
    z = int(value_d)
    if not (2.5 <= a_ang <= 8.0 and 0.0001 <= conc <= 0.1 and 0 <= t_c <= 100 and z in (1, 2)
            and float(value_d) == z):
        raise ValueError("inputs out of domain")
    a_m = a_ang * M_PER_ANGSTROM

    # S1: solvent permittivity -> Bjerrum pair length and Debye-Hueckel constants
    length, a_dh, b_dh = electrostatics(t_c, z)

    # S2: Bjerrum association constant (zero below the critical distance)
    k_a = association_constant(a_m, length)

    # S3: free-ion fraction alpha from mass action with extended-DH free-ion activities
    def log_gamma_free(alpha):
        s = math.sqrt(alpha * conc * z * z)
        return -a_dh * z * z * s / (1 + b_dh * a_m * s)

    lo, hi = 0.0, 1.0
    for _ in range(N_BISECT):
        mid = 0.5 * (lo + hi)
        g = 10 ** log_gamma_free(mid)
        if k_a * mid * mid * conc * g * g - (1 - mid) > 0:
            hi = mid
        else:
            lo = mid
    alpha = 0.5 * (lo + hi)

    # S4: stoichiometric mean activity coefficient
    return round(math.log10(alpha) + log_gamma_free(alpha), 4)


def oracle(inputs):
    return {"output": compute(inputs["value_a"], inputs["value_b"], inputs["value_c"], inputs["value_d"])}


if __name__ == "__main__":
    import json
    import sys
    print(json.dumps(oracle(json.loads(sys.stdin.read()))))
