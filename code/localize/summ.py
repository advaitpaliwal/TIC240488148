import json, sys, numpy as np
for s in map(int,sys.argv[1:]):
    d=json.load(open(f'out/diff_s{s}.json')); st=d['stars']; t=[r for r in st if r['sep']<0.5][0]
    print(f"S{s} events {d['n_events']} cad {d['cadence_min']:.1f} ap3x3 {d['ap3x3_depth_ppm']:.0f}+-{d['ap3x3_err_ppm']:.0f} best {d['best']} target chi2 {t['chi2']:.1f} amp {t['amp']:.2f}+-{t['amp_err']:.2f}")
    for r in sorted(st,key=lambda r:r['chi2'])[:6]+[r for r in st if r['source'] in (1967569190261629440,1967568812304507776,1967568777945007744)]:
        print(f"   {r['source']} G{r['Gmag']:.2f} sep {r['sep']:.1f} pa {r['pa']:.0f} xy ({r['x']:.2f},{r['y']:.2f}) chi2 {r['chi2']:.1f} dchi2 {r['chi2']-t['chi2']:+.1f} amp {r['amp']:.2f}+-{r['amp_err']:.2f}")
