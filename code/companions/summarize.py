import json, glob, numpy as np, pandas as pd
rows=[json.load(open(f)) for f in sorted(glob.glob('trc/out/res_*.json'))]
d=pd.DataFrame(rows)
for m,g in d.groupby('mode'):
    sc=[]
    for k in g['run']:
        p=pd.read_csv(f'trc/out/probs_{m}_{k}.csv'); sc.append(p.groupby('scenario')['prob'].sum())
    s=pd.concat(sc,axis=1).mean(axis=1)
    print(f"{m}: n={len(g)} FPP {g.FPP.mean():.4f} +/- {g.FPP.std(ddof=1) if len(g)>1 else float('nan'):.4f}  NFPP {g.NFPP.mean():.5f} +/- {g.NFPP.std(ddof=1) if len(g)>1 else float('nan'):.5f}  "
          f"TP {s.get('TP',0):.3f} PTP {s.get('PTP',0):.3f} DTP {s.get('DTP',0):.3f} STP {s.get('STP',0):.3f} SEB {s.get('SEB',0):.4f} PEB {s.get('PEB',0):.4f}  FPPs {np.round(g.FPP.values,4).tolist()}")
