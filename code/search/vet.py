# adapted from ../vet.py: folded plot + per-sector depths for FFI light curves
import sys, json, numpy as np, pandas as pd
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
tic=int(sys.argv[1]); r=json.load(open(f'res/{tic}.json'))
star=pd.read_csv('pool_all.csv').set_index('ID').loc[tic]
d=np.load(f'lc/{tic}.npz'); t,f,sid=d['t'],d['f'],d['s']
P=float(sys.argv[2]) if len(sys.argv)>2 else r['P']; T0=float(sys.argv[3]) if len(sys.argv)>3 else r['T0']; D=float(sys.argv[4]) if len(sys.argv)>4 else None
G=6.674e-11; a=(G*star['mass']*1.989e30*(P*86400)**2/(4*np.pi**2))**(1/3)
Dexp=P/np.pi*np.arcsin(min(1,star.rad*6.957e8/a))
if D is None: D=max(0.85*Dexp,0.03)
ph=((t-T0+0.5*P)%P)-0.5*P; n=np.round((t-T0)/P).astype(int)
intr=np.abs(ph)<D/2*0.8; oot=(np.abs(ph)>D)&(np.abs(ph)<3*D)
def dep(m):
    a=f[m&intr]; b=f[m&oot]
    if len(a)<2 or len(b)<5: return [np.nan,np.nan,len(a)]
    return [float(np.median(b)-np.mean(a)), float(np.std(b)/np.sqrt(len(a))), int(len(a))]
out={'tic':tic,'P':P,'T0':T0,'dur_h':D*24,'used':r['used']}
out['depth_all']=dep(np.ones_like(f,bool)); out['odd']=dep(n%2==1); out['even']=dep(n%2==0)
o,e=out['odd'],out['even']; out['oe_sigma']=abs(o[0]-e[0])/np.hypot(o[1],e[1])
ph2=((t-T0)%P)-0.5*P; si=np.abs(ph2)<D/2; so=(np.abs(ph2)>D)&(np.abs(ph2)<3*D)
out['secondary']=[float(np.median(f[so])-np.mean(f[si])), float(np.std(f[so])/np.sqrt(si.sum()))]
out['per_sector']={int(s):dep(sid==s) for s in np.unique(sid)}
ev={}
for k in np.unique(n[intr]):
    m=n==k; a=f[m&intr]; b=f[m&oot]
    if len(a)>=2 and len(b)>=5: ev[int(k)]=[round(float((np.median(b)-np.mean(a))*1e6)),round(float(np.std(b)/np.sqrt(len(a))*1e6)),round(float(T0+k*P),3),int(sid[m][0])]
out['events']=ev
G=6.674e-11; a=(G*star['mass']*1.989e30*(P*86400)**2/(4*np.pi**2))**(1/3)
out['dur_expected_h']=P*24/np.pi*np.arcsin(min(1,star.rad*6.957e8/a))
out['Rp_earth']=np.sqrt(max(out['depth_all'][0],0))*star.rad*109.1
out['star']=dict(R=star.rad,M=star['mass'],Teff=star.Teff,Tmag=star.Tmag,d=star.d)
json.dump(out,open(f'vet_{tic}.json','w'),indent=1,default=float); print(json.dumps(out,default=float))
us=np.unique(sid); fig,ax=plt.subplots(2+len(us),1,figsize=(8,4+2.2*len(us)))
bins=np.linspace(-3*D,3*D,41); bc=(bins[1:]+bins[:-1])/2*24
def binned(x,y): idx=np.digitize(x,bins); return [np.mean(y[idx==i]) if (idx==i).sum()>0 else np.nan for i in range(1,len(bins))]
w=np.abs(ph)<3*D; ax[0].plot(ph[w]*24,f[w],'.',c='0.7',ms=2); ax[0].plot(bc,binned(ph,f),'o-k',ms=4)
dd=max(out['depth_all'][0],2e-4); yl=(1-4*dd,1+2.5*dd)
ax[0].set_ylim(yl); ax[0].set_title(f'TIC {tic}  P={P:.5f} d  depth {out["depth_all"][0]*1e6:.0f} ppm  all sectors (FFI)')
ax[1].plot(bc,binned(ph2,f),'o-k',ms=4); ax[1].set_ylim(yl); ax[1].set_title('phase 0.5 (secondary)')
for i,s in enumerate(us):
    m=sid==s; ax[2+i].plot(ph[m&w]*24,f[m&w],'.',c='0.7',ms=2); ax[2+i].plot(bc,binned(ph[m],f[m]),'o-k',ms=4); ax[2+i].set_ylim(yl)
    ps=out['per_sector'][int(s)]; ax[2+i].set_title(f'S{s}: {ps[0]*1e6:.0f} +- {ps[1]*1e6:.0f} ppm',fontsize=9)
ax[-1].set_xlabel('hours from mid-transit'); plt.tight_layout(); plt.savefig(f'vet_{tic}.png',dpi=90)
