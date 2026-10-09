import warnings; warnings.filterwarnings('ignore')
import pickle, json, numpy as np, pandas as pd, sys
sys.path.insert(0,'work/hunt-ffi/trc')
H='work/hunt-ffi/trc/'
def main():
    tg=pickle.load(open(H+'out/target_240488148.pkl','rb'))
    r=json.load(open(H+'out/res_240488148_0.json'))
    tg.calc_depths(tdepth=r['depth_ppm']*1e-6, all_ap_pixels=tg._aps)
    s=tg.stars
    s[['ID','Tmag','GAIAmag','sep (arcsec)','PA (E of N)','fluxratio','tdepth','plx']].to_csv('out/trc_depths.csv',index=False)
    w=(s.tdepth>0)
    print('stars with tdepth>0:',w.sum(),'of',len(s))
    print(s[w][['ID','Tmag','sep (arcsec)','PA (E of N)','fluxratio','tdepth']].sort_values('sep (arcsec)').to_string())
if __name__=='__main__': main()
