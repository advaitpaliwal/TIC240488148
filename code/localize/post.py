# For each sector: fit OOT image to get target flux; for each TIC neighbour (TRICERATOPS list, tdepth>0) compute required flux deficit
# A_req = tdepth_k * F_k and test it against the in/out difference image (forced PRF fit at the star's Gaia position).
import warnings; warnings.filterwarnings('ignore')
import sys, json, numpy as np, pandas as pd
from common import *
def main():
    from PRF import TESS_PRF
    from astropy.coordinates import SkyCoord
    trc=pd.read_csv('out/trc_depths.csv'); trc=trc[trc.tdepth>0].reset_index(drop=True)
    TT=trc.Tmag[0]; rows=[]
    for s in map(int,sys.argv[1:]):
        tpf=tpf_for(s); z=np.load(f'out/diff_s{s}.npz'); D,sig,oot=z['D'],z['sig'],z['oot']; ny,nx=D.shape
        prf=TESS_PRF(tpf.camera,tpf.ccd,s,tpf.column+nx/2,tpf.row+ny/2)
        g=gaia_field(tpf); gc=SkyCoord(g.RA_ICRS.values*u.deg,g.DE_ICRS.values*u.deg)
        w=1/sig.ravel()**2
        def lin(models,img,W):
            A=np.vstack([m.ravel() for m in models]+[np.ones(img.size)]).T
            C=np.linalg.inv(A.T@(A*W[:,None])); b=C@(A.T@(img.ravel()*W)); r=img.ravel()-A@b
            return b,np.sqrt(np.diag(C)),float((r**2*W).sum())
        # target total flux from OOT image (PRF + const), photon-ish weights
        t0=g.iloc[0]; mT=prf.locate(t0.x,t0.y,(ny,nx))
        bT,_,_=lin([mT],oot,1/(np.abs(oot.ravel())+100)); FT=bT[0]
        bD,eD,c2T=lin([mT],D,w)   # target-only fit of difference image
        from astropy.io import fits
        with fits.open(f'work/hunt-ffi/lc/240488148/spoc_s{s}.fits') as h:
            a=h['APERTURE'].data.astype(int); c0=h['APERTURE'].header['CRVAL1P']; r0=h['APERTURE'].header['CRVAL2P']; crowd=h[1].header.get('CROWDSAP',1.0)
        rr,cc=np.nonzero(a&2); ap=np.zeros((ny,nx),bool)
        for R_,C_ in zip(rr+r0-tpf.row,cc+c0-tpf.column):
            if 0<=R_<ny and 0<=C_<nx: ap[R_,C_]=True
        dF_ap=366e-6*crowd*oot[ap].sum(); fT=mT[ap].sum()
        print('S',s,'ap npix',ap.sum(),'of',len(rr),'crowd',crowd,'dF_ap',round(dF_ap,2),'A_req target',round(dF_ap/fT,2),'target A fit',round(bD[0],2),'+-',round(eD[0],2),flush=True)
        for i,r in trc.iterrows():
            ic=SkyCoord(r.ra*u.deg,r.dec*u.deg) if 'ra' in r else None
            # match TIC star to Gaia by position (TIC coords are epoch 2000) and G mag
            sep=gc.separation(SkyCoord(*(lambda q:(q.ra,q.dec))(SkyCoord(0*u.deg,0*u.deg))) ) if False else None
            rows.append(None) if False else None
        trcpos=pd.read_csv('out/trc_depths_pos.csv')
        for i,r in trcpos.iterrows():
            d=gc.separation(SkyCoord(r.ra*u.deg,r.dec*u.deg)).arcsec; d=np.where(np.abs(g.Gmag-r.GAIAmag)<0.5,d,999); j=int(np.argmin(d))
            if d[j]>3: continue
            x,y=g.x[j],g.y[j]
            if not (-0.5<x<nx-0.5 and -0.5<y<ny-0.5): continue
            m=prf.locate(x,y,(ny,nx)); Fk=FT*10**(-0.4*(r.Tmag-TT)); fk=m[ap].sum()
            Areq=dF_ap/fk if fk>0 else np.inf
            b,e,c2=lin([m],D,w)
            c2req=float((((D.ravel()-Areq*m.ravel())-np.median(D.ravel()-Areq*m.ravel()))**2*w).sum())
            rows.append(dict(sector=s,ID=str(r.ID),Tmag=r.Tmag,sep=r['sep (arcsec)'],tdepth_trc=r.tdepth,Fk=Fk,ap_frac=fk,depth_req=Areq/Fk,FT=FT,A_req=Areq,A_fit=b[0],A_err=e[0],
                             z_req=(Areq-b[0])/e[0],chi2_free=c2,chi2_target=c2T,dchi2_free=c2-c2T,chi2_forced=c2req,dchi2_forced=c2req-c2T,targetA=bD[0],targetA_err=eD[0]))
    df=pd.DataFrame(rows); df.to_csv('out/neighbour_tests.csv',index=False)
    pd.set_option('display.width',250); print(df.round(2).to_string())
if __name__=='__main__': main()
