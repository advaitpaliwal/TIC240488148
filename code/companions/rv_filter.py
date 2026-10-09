# Inferred extra test: Gaia DR3 RVS constant (rv_chisq_pvalue=0.828, rv_amplitude_robust=2.58 km/s, 32 transits, e_RV=0.18 km/s)
import pandas as pd, numpy as np, glob
A=pd.concat([pd.read_csv(f) for f in sorted(glob.glob('molusc_out/m50k_s1?_all.csv'))],ignore_index=True)
M=0.69; G=2.959122e-4 # AU^3/(Msun d^2)
q=A['mass ratio']; P=A['period(days)']; e=A['eccentricity']; ci=A['cos_i']; m2=q*M
K=(2*np.pi*G/P)**(1/3)*m2*np.sqrt(1-ci**2)/(M+m2)**(2/3)/np.sqrt(1-e**2)*1731.46  # AU/d -> km/s
sig_ep=np.sqrt((0.18**2-0.11**2)*2*32/np.pi); print('per-epoch sigma est %.2f km/s'%sig_ep)
rv_rej=(P<1038)&(K>3*sig_ep)
kept=A['Full Rejected']==0
print('kept before RV %.4f, after RV %.4f'%(kept.mean(),(kept&~rv_rej).mean()))
st=m2>=0.08
print('stellar kept after RV, frac of stellar generated %.4f'%((kept&~rv_rej&st).sum()/st.sum()))
w=A['semi-major axis(AU)']*(1-e)>10
print('RV-rejected among a(1-e)>10AU kept: %d'%((kept&rv_rej&w).sum()))
for f in (0.41,0.46):
    for s,lab in [((kept&~rv_rej).mean(),'all'),((kept&~rv_rej&st).sum()/st.sum()*st.mean(),'stellar')]:
        print('f=%.2f %s: s=%.3f posterior P(undetected companion)=%.3f'%(f,lab,s,f*s/(f*s+1-f)))
# separation/ contrast distribution of surviving stellar companions
k=A[kept&~rv_rej&st]
print('surviving stellar: median a %.1f AU; frac a<1 AU %.3f; frac proj sep<19mas(1.1AU) %.3f'%(k['semi-major axis(AU)'].median(),(k['semi-major axis(AU)']<1).mean(),(k['Projected Separation(AU)']<1.1).mean()))
print('surviving stellar: frac with DeltaG<3 %.3f, <5 %.3f'%((k['DeltaG']<3).mean(),(k['DeltaG']<5).mean()))
print('generated stellar: frac with DeltaG<3 %.3f'%((A[st]['DeltaG']<3).mean()))
