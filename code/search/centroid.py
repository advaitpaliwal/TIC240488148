# per-sector flux-weighted centroid shift in transit vs random phases (TESS-SPOC MOM_CENTR)
import sys, glob, numpy as np, warnings; warnings.filterwarnings('ignore')
from astropy.io import fits
tic=int(sys.argv[1]); P=float(sys.argv[2]); T0=float(sys.argv[3]); D=float(sys.argv[4])/24
rng=np.random.default_rng(0)
for fn in sorted(glob.glob(f'lc/{tic}/spoc_s*.fits')):
    h=fits.open(fn); d=h[1].data; q=d['QUALITY']==0
    t=d['TIME'][q]; cx=d['MOM_CENTR1'][q]; cy=d['MOM_CENTR2'][q]; ok=np.isfinite(t)&np.isfinite(cx)&np.isfinite(cy); t,cx,cy=t[ok],cx[ok],cy[ok]
    def det(c):  # running median 0.75 d detrend
        out=np.empty_like(c)
        for i in range(len(t)): w=np.abs(t-t[i])<0.375; out[i]=c[i]-np.median(c[w])
        return out
    cx,cy=det(cx),det(cy)
    def shift(off):
        ph=((t-T0-off*P+0.5*P)%P)-0.5*P; a=np.abs(ph)<D/2; b=(np.abs(ph)>D)&(np.abs(ph)<D+0.15)
        return np.mean(cx[a])-np.mean(cx[b]), np.mean(cy[a])-np.mean(cy[b])
    s=np.array(shift(0)); r=np.array([shift(x) for x in rng.uniform(0.2,0.8,200)])
    z=(s-r.mean(0))/r.std(0)
    print(f'{fn.split("/")[-1]}: dx {s[0]*21e3:.1f} dy {s[1]*21e3:.1f} mas; z {z[0]:.2f} {z[1]:.2f}')
