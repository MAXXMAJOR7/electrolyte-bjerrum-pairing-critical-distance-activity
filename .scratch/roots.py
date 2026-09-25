import math,random,sys
sys.path.insert(0,'../oracle'); import implement as o
rng=random.Random(2); bad=0; nonmono=0
for trial in range(1500):
    a=rng.uniform(2.5,8); c=10**rng.uniform(-4,-1); t=rng.uniform(0,100); z=rng.choice([1,2])
    l,A,B=o.electrostatics(t,z); KA=o.association_constant(a*1e-10,l)
    if KA==0: continue
    am=a*1e-10; prev=None; ch=0; dec=0
    for i in range(0,4001):
        al=10**(-6+6*i/4000); s=math.sqrt(al*c*z*z); g=10**(-A*z*z*s/(1+B*am*s))
        F=KA*al*al*c*g*g-(1-al)
        if prev is not None:
            ch+= (F>0)!=(prev>0); dec+= F<prev
        prev=F
    bad+= ch!=1; nonmono+= dec>0
print("multi-root",bad,"nonmono",nonmono)
