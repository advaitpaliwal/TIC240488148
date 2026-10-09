# per-sector mid-transit offset scan with fixed shape, to check the linear ephemeris
import numpy as np, sys
sys.argv=["x"]; sys.path.insert(0,"."); import fit as F_
import batman
th=np.array([F_.TREF0,F_.P0,0.0195,0.3,F_.RHO,0.43,0.36,1.0,-9.6,-10.6])
for s in F_.SECS:
    m=F_.S==s; t=F_.T[m]; f=F_.F[m]; e=F_.E[m]
    ex=1800 if s<20 else 200; ss=15 if s<20 else 3
    n=np.round((np.median(t)-F_.T00)/F_.P0); tc0=F_.T00+n*F_.P0
    p=F_.p; T0,P,k,b,rho,a,inc,u1,u2,f0,_,_=F_.unpack(th)
    p.per=P;p.rp=k;p.a=a;p.inc=inc;p.u=[u1,u2]
    mod=batman.TransitModel(p,t,supersample_factor=ss,exp_time=ex/86400)
    dts=np.linspace(-0.08,0.08,321); chi=[]
    sd=np.std(f)
    for dt in dts:
        p.t0=tc0+dt; chi.append(np.sum((f-mod.light_curve(p))**2)/sd**2)
    chi=np.array(chi); i=np.argmin(chi); ok=dts[chi<chi[i]+1]
    print(s,int(n),'best dt min %.1f'%(dts[i]*1440),'1sig range %.1f..%.1f'%(ok.min()*1440,ok.max()*1440),'dchi2 vs none %.1f'%(np.sum((f-1)**2)/sd**2-chi[i]))
