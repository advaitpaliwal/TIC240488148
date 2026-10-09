import warnings; warnings.filterwarnings('ignore')
import numpy as np, glob, lightkurve as lk, astropy.units as u
from astropy.coordinates import SkyCoord
from astroquery.vizier import Vizier
P,T0,DUR=1.512529,1712.9250,1.5/24
RA0,DE0=323.628288,42.644034
def tpf_for(s):
    p=glob.glob(f'cut/tesscut/tess-s{s:04d}-*_15x15_astrocut.fits')[0]
    return lk.read(p)
def gaia_field(tpf, rad=150):
    import os, pandas as pd
    fn='gaia_field.csv'
    if not os.path.exists(fn):
        v=Vizier(columns=['Source','RA_ICRS','DE_ICRS','pmRA','pmDE','Gmag','Plx'],row_limit=-1)
        t=v.query_region(SkyCoord(RA0*u.deg,DE0*u.deg),radius=rad*u.arcsec,catalog='I/355/gaiadr3')[0].to_pandas()
        t.to_csv(fn,index=False)
    t=pd.read_csv(fn)
    dy=(tpf.time[0].jd-2457389.0)/365.25  # Gaia DR3 epoch 2016.0 = JD 2457389.0
    ra=(t.RA_ICRS+np.nan_to_num(t.pmRA)*dy/3.6e6/np.cos(np.radians(t.DE_ICRS))).values; de=(t.DE_ICRS+np.nan_to_num(t.pmDE)*dy/3.6e6).values
    x,y=tpf.wcs.all_world2pix(np.vstack([ra,de]).T,0).T
    c=SkyCoord(ra*u.deg,de*u.deg); i0=np.argmin(np.where(np.abs(t.Gmag.values-10.88)<0.2,c.separation(SkyCoord(RA0*u.deg,DE0*u.deg)).arcsec,1e9))
    t=t.assign(ra_ep=ra,de_ep=de,x=x,y=y,sep=c.separation(c[i0]).arcsec,pa=c[i0].position_angle(c).deg)
    return t.sort_values('sep').reset_index(drop=True)
