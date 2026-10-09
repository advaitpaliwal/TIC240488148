import pandas as pd, numpy as np, glob
A=pd.concat([pd.read_csv(f) for f in sorted(glob.glob('molusc_out/m50k_s1?_all.csv'))],ignore_index=True)
A=A[A['mass ratio']*0.69>=0.08]
sb=[0,1.1,5,20,60,130,250,1e5]; gb=[-0.2,1,2,3,4,5,7.1,99]
A['s']=pd.cut(A['Projected Separation(AU)'],sb); A['g']=pd.cut(A['DeltaG'],gb)
t=A.groupby(['s','g'],observed=True)['Full Rejected'].mean().unstack().round(2)
print('Fraction of simulated stellar companions REJECTED by Gaia RUWE + Gaia resolved-source test (rows: projected sep AU; 1 arcsec = 58.8 AU)')
print(t.to_string())
n=A.groupby(['s'],observed=True)['Full Rejected'].agg(['mean','size']).round(3); print(n.to_string())
