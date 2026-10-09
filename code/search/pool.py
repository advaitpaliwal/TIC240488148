# pool: first-hunt TIC dwarfs (tic_sel.csv) never on any 2-min target list S1-S100, not in TOI/CTOI/NEA
import pandas as pd, glob, numpy as np
H='../'
tic=pd.read_csv(H+'tic_sel.csv')
excl=set(pd.read_csv(H+'_sources/toi.csv')['TIC ID'].astype(int))
c=pd.read_csv(H+'_sources/ctoi.csv'); excl|=set(c[c.columns[0]].astype(int))
p=pd.read_csv(H+'_sources/ps.csv').dropna(subset=['tic_id']); excl|=set(p.tic_id.str.replace('TIC ','').astype(int))
two=set()
for f in sorted(glob.glob(H+'_sources/tl/S*.txt')):
    d=pd.read_csv(f,comment='#',sep=r'\s+',header=None,usecols=[0]); two|=set(d[0].astype(int))
print('tl files',len(glob.glob(H+'_sources/tl/S*.txt')),'2min TICs',len(two))
tic['two']=tic.ID.isin(two); tic['known']=tic.ID.isin(excl)
print('pool',len(tic),'ever 2-min',tic.two.sum(),'known',tic.known.sum())
t=tic[~tic.two & ~tic.known].copy()
print('no-2min & unknown',len(t))
# sectors from tess-point (S1-S134 pointing table); keep observed sectors <= 100 later via product availability
from tess_stars2px import tess_stars2px_function_entry as s2p
out=s2p(t.ID.values,t.ra.values,t.dec.values)
ids,secs=out[0],out[3]
m={}
for i,s in zip(ids,secs):
    if s<=110 or s>1000: m.setdefault(int(i),set()).add(int(s))
t['sectors']=t.ID.map(lambda i: sorted(m.get(int(i),[])))
t['ns']=t.sectors.map(len)
print(t.ns.value_counts().sort_index().to_string())
t.to_csv('pool_all.csv',index=False)
