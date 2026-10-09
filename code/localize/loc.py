# TESS-Localize on a TESScut FFI stamp, using harmonics n/P of the transit period (Higgins & Bell 2023 treat EBs/transits as harmonic series).
from common import *
def main():
    import sys, json, numpy as np, astropy.units as u, matplotlib; matplotlib.use('Agg')
    import tess_localize as tl
    s=int(sys.argv[1]); N=int(sys.argv[2]) if len(sys.argv)>2 else 25
    tpf=tpf_for(s)
    ap=np.zeros(tpf.shape[1:],bool); ap[6:9,6:9]=True
    fr=np.arange(1,N+1)/P
    L=tl.Localize(targetpixelfile=tpf,frequencies=list(fr),frequnit=u.d**-1,principal_components=3,aperture=ap,gaia=True,magnitude_limit=21)
    sf=L.starfit
    g=gaia_field(tpf)
    tgt=sf.iloc[(np.hypot(sf.ra-RA0,(sf.dec-DE0))).argmin()]['source'] if False else None
    from astropy.coordinates import SkyCoord
    c=SkyCoord(sf.ra.values*u.deg,sf.dec.values*u.deg); sep=c.separation(SkyCoord(RA0*u.deg,DE0*u.deg)).arcsec
    sf=sf.assign(sep_arcsec=sep)
    sf.to_csv(f'out/loc_s{s}_starfit.csv',index=False)
    loc=SkyCoord(L.location_skycoord[0]*u.deg,L.location_skycoord[1]*u.deg)
    hh=np.array([L.result.params[f'height{j}'].value for j in range(N)]); he=np.array([L.result.params[f'height{j}'].stderr or np.nan for j in range(N)])
    out=dict(sector=s,N=N,loc_x=float(L.location[0]),loc_y=float(L.location[1]),loc_ra=float(L.location_skycoord[0]),loc_dec=float(L.location_skycoord[1]),
     loc_sep_from_target_arcsec=float(loc.separation(SkyCoord(RA0*u.deg,DE0*u.deg)).arcsec),
     x_err=float(L.result.params['column'].stderr or np.nan),y_err=float(L.result.params['row'].stderr or np.nan),
     height_snr=float(np.nansum(hh)/np.sqrt(np.nansum(he**2))),
     top=sf.head(8)[['source','Gmag','sep_arcsec','distance','pvalue','relative likelihood']].to_dict('records'))
    json.dump(out,open(f'out/loc_s{s}.json','w'),indent=1,default=float)
    L.plot(save=f'out/loc_s{s}.png') if False else None
    import matplotlib.pyplot as plt
    fig,ax=plt.subplots(figsize=(6,6)); ax.imshow(L.heatmap,origin='lower',cmap='viridis')
    w=sf.sep_arcsec<130; ax.scatter(sf.x[w],sf.y[w],s=np.clip(200/2**((sf.Gmag[w]-10)/2),2,200),facecolor='none',edgecolor='w')
    ax.plot(L.location[0],L.location[1],'rx',ms=12,mew=2); ax.set_title(f'S{s} TESS-Localize SNR heatmap, X=fit'); fig.savefig(f'out/loc_s{s}.png',dpi=90)
    print(json.dumps(out,default=float)[:600])
    
if __name__=='__main__':
    main()