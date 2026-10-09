import sys, numpy as np
tic=int(sys.argv[1]); P=float(sys.argv[2]); T0=float(sys.argv[3]); prod=sys.argv[4] if len(sys.argv)>4 else 'TESS-SPOC'
z=np.load(f'lc/{tic}_{prod}_masked.npz'); t,f=z['t'],z['f']
# refine T0/P/duration by grid on box SNR
best=None
for dP in np.linspace(-0.0004,0.0004,17):
  for dT in np.linspace(-0.03,0.03,13):
    for D in np.arange(0.5,4.01,0.25)/24:
      ph=((t-T0-dT+0.5*(P+dP))%(P+dP))-0.5*(P+dP); a=np.abs(ph)<D/2; b=(np.abs(ph)>D)&(np.abs(ph)<D+0.15)
      d=np.median(f[b])-np.mean(f[a]); e=np.std(f[b])/np.sqrt(a.sum()); s=d/e
      if best is None or s>best[0]: best=(s,P+dP,T0+dT,D,d,e)
s,P,T0,D,d,e=best; print(f'best SNR {s:.2f} P {P:.6f} T0 {T0:.4f} dur_h {D*24:.2f} depth {d*1e6:.0f}+-{e*1e6:.0f}')
ph=((t-T0+0.5*P)%P)-0.5*P; n=np.round((t-T0)/P).astype(int); a=np.abs(ph)<D/2; b=(np.abs(ph)>D)&(np.abs(ph)<D+0.15)
def dep(m): return (np.median(f[b&m])-np.mean(f[a&m]))*1e6, np.std(f[b&m])/np.sqrt((a&m).sum())*1e6
o=dep(n%2==1); ev=dep(n%2==0); print('odd %.0f+-%.0f even %.0f+-%.0f  diff %.2f sigma'%(o+ev+(abs(o[0]-ev[0])/np.hypot(o[1],ev[1]),)))
ph2=((t-T0)%P)-0.5*P; a2=np.abs(ph2)<D/2; b2=(np.abs(ph2)>D)&(np.abs(ph2)<D+0.15)
print('secondary %.0f+-%.0f'%((np.median(f[b2])-np.mean(f[a2]))*1e6,np.std(f[b2])/np.sqrt(a2.sum())*1e6))
for half in (0.25,0.5):  # ingress/egress shape: depth in inner vs outer half of box
  pass
inner=np.abs(ph)<D/4; outer=(np.abs(ph)>D/4)&(np.abs(ph)<D/2)
print('inner %.0f outer %.0f ppm (U-shape: similar; V: outer much smaller)'%((np.median(f[b])-np.mean(f[inner]))*1e6,(np.median(f[b])-np.mean(f[outer]))*1e6))
for k in (0.5,2,3,1/3,1.5):
  Pk=P*k; ph=((t-T0+0.5*Pk)%Pk)-0.5*Pk; a=np.abs(ph)<D/2; b=(np.abs(ph)>D)&(np.abs(ph)<D+0.15)
  print(f'P x {k:.2f}: depth {(np.median(f[b])-np.mean(f[a]))*1e6:.0f} +- {np.std(f[b])/np.sqrt(a.sum())*1e6:.0f}')
