import pandas as pd, numpy as np, astropy.units as u
from astropy.coordinates import SkyCoord
from astroquery.vizier import Vizier
t=pd.read_csv('pool_all.csv'); t=t[t.ns>=3].reset_index(drop=True)
v=Vizier(columns=['Source','RA_ICRS','DE_ICRS','Gmag','RUWE','_r'],row_limit=-1)
res=[]
for i in range(0,len(t),200):
    b=t.iloc[i:i+200]
    r=v.query_region(SkyCoord(b.ra.values*u.deg,b.dec.values*u.deg),radius=30*u.arcsec,catalog='I/355/gaiadr3')[0].to_pandas()
    r['tic']=b.ID.values[r['_q'].values-1]; res.append(r); print(i,len(r),flush=True)
r=pd.concat(res); r.to_csv('_sources/gaia_30as_vizier.csv',index=False)
fl=[]
for tic,g in r.groupby('tic'):
    g=g.sort_values('_r'); near=g[g._r<5]
    if len(near)==0: fl.append(dict(tic=tic)); continue
    tg=near.sort_values('Gmag').iloc[0]
    oth=g[g.Source!=tg.Source]
    fl.append(dict(tic=tic,gaia=tg.Source,G=tg.Gmag,ruwe=tg.RUWE,sep_t=tg._r,n_bright_30=int(((oth.Gmag-tg.Gmag)<5).sum()),nearest_bright=float(oth[(oth.Gmag-tg.Gmag)<5]._r.min()) if ((oth.Gmag-tg.Gmag)<5).any() else np.nan))
f=pd.DataFrame(fl); f.to_csv('gaia_flags.csv',index=False)
print(len(f),'ruwe<1.2',(f.ruwe<1.2).sum(),'clean',((f.ruwe<1.2)&(f.n_bright_30==0)).sum())
