import math
E=1.602176634e-19;KB=1.380649e-23;NA=6.02214076e23;EPS0=8.8541878188e-12
def eps(t): return 87.740-0.40008*t+9.398e-4*t**2-1.410e-6*t**3
def Q(b,n=4000):
    if b<=2: return 0.0
    h=(b-2)/n; f=lambda x: math.exp(x)*x**-4
    s=f(2)+f(b)+sum((4 if i%2 else 2)*f(2+i*h) for i in range(1,n))
    return s*h/3
def run(a_A,c,t,z,pair=True,alpha_out=True,bjerrum=True):
    T=t+273.15; er=eps(t); ekt=4*math.pi*EPS0*er*KB*T
    l=z*z*E*E/ekt; a=a_A*1e-10; b=l/a
    kap_per_sqrtI=math.sqrt(2*NA*1000*E*E/(EPS0*er*KB*T))  # m^-1 per sqrt(mol/L)
    A=E*E*kap_per_sqrtI/(8*math.pi*EPS0*er*KB*T)/math.log(10)
    Bq=kap_per_sqrtI
    if bjerrum: KA=4*math.pi*NA*l**3*Q(b)*1000
    else: KA=4*math.pi*NA*a**3*math.exp(b)/3*1000
    if not pair: KA=0
    def lg(al):
        I=al*c*z*z; s=math.sqrt(I); return -A*z*z*s/(1+Bq*a*s)
    lo,hi=0.0,1.0
    for _ in range(200):
        al=0.5*(lo+hi); g=10**lg(al)
        # f = KA*al^2 c g^2 - (1-al), increasing in al
        if KA*al*al*c*g*g-(1-al)>0: hi=al
        else: lo=al
    al=0.5*(lo+hi)
    return (math.log10(al) if alpha_out else 0)+lg(al), al, KA, b, A, Bq*1e-10
for args in [(2.5,0.1,25,1),(3.0,0.1,25,1),(4,0.1,25,1),(4,0.1,100,1),(4,0.01,25,2),(4,0.1,25,2),(8,0.1,25,2),(4,0.01,25,3),(3,0.05,100,3),(8,1e-4,0,1),(2.5,0.1,100,3)]:
    r=run(*args); print(args, ["%.5g"%x for x in r], "noPair %.4f"%run(*args,pair=False)[0], "fuoss %.4f"%run(*args,bjerrum=False)[0], "free %.4f"%run(*args,alpha_out=False)[0])
