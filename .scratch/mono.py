import math,random
from proto import *
import proto
worst=0
rng=random.Random(1)
for trial in range(3000):
    a=rng.uniform(2.5,8); c=10**rng.uniform(-4,-1); t=rng.uniform(0,100); z=rng.choice([1,2,3])
    T=t+273.15; er=eps(t); ekt=4*math.pi*EPS0*er*KB*T; l=z*z*E*E/ekt; b=l/(a*1e-10)
    kap=math.sqrt(2*NA*1000*E*E/(EPS0*er*KB*T)); A=E*E*kap/(8*math.pi*EPS0*er*KB*T)/math.log(10)
    KA=4*math.pi*NA*l**3*Q(b,400)*1000
    if KA==0: continue
    prev=None; changes=0
    for i in range(1,20001):
        al=10**(-8+8*i/20000); I=al*c*z*z; s=math.sqrt(I); g=10**(-A*z*z*s/(1+kap*a*1e-10*s))
        F=KA*al*al*c*g*g-(1-al)
        if prev is not None and (F>0)!=(prev>0): changes+=1
        prev=F
    if changes!=1: print("nonmono",a,c,t,z,changes)
print("done")
