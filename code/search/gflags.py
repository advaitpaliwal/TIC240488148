import pandas as pd, glob, numpy as np, io, os
fl=[]
for f in glob.glob('_sources/gaia/*.tsv'):
    tic=int(os.path.basename(f)[:-4])
    L=[l for l in open(f) if l.strip() and not l.startswith('#')]
    if len(L)<3: fl.append(dict(tic=tic)); continue
    g=pd.read_csv(io.StringIO(''.join(L)),sep='\t',skiprows=[1,2])
    g.columns=[c.strip() for c in g.columns]
    g=g[pd.to_numeric(g['_r'],errors='coerce').notna()]
    g=g[g['_r'].astype(str).str.strip().str.match(r'^[0-9.]+$')]
    for c in ['_r','Gmag','RUWE']: g[c]=pd.to_numeric(g[c],errors='coerce')
    g=g[g.Gmag.notna()]
    near=g[g._r<5]
    if len(near)==0: fl.append(dict(tic=tic)); continue
    tg=near.sort_values('Gmag').iloc[0]; oth=g[g.Source!=tg.Source]; b=oth[(oth.Gmag-tg.Gmag)<5]
    fl.append(dict(tic=tic,gaia=tg.Source,G=tg.Gmag,ruwe=tg.RUWE,sep_t=tg._r,n_bright_30=len(b),nearest_bright=b._r.min() if len(b) else np.nan))
f=pd.DataFrame(fl); f.to_csv('gaia_flags.csv',index=False)
print(len(f),'ruwe<1.2',(f.ruwe<1.2).sum(),'clean',((f.ruwe<1.2)&(f.n_bright_30==0)).sum())
