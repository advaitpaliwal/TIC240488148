# Archive sweep for TIC 240488148 / Gaia DR3 1967568881023980032; saves raw tables in _sources/
import warnings; warnings.filterwarnings('ignore')
import astropy.units as u
from astropy.coordinates import SkyCoord
from astroquery.vizier import Vizier
c=SkyCoord(323.62848767584302*u.deg, 42.643554860977403*u.deg)
cats={'gaiadr3':'I/355/gaiadr3','gaiadr3_param':'I/355/paramp','nss_twobody':'I/357/tbooc','nss_acc':'I/357/tbosb1','nss_sb1':'I/357/tbosb1','nss_sb2':'I/357/tbosb2','nss_accel':'I/357/accel7','nss_accel9':'I/357/accel9','nss_vim':'I/357/vimf','lamost_dr5_lrs':'V/164','lamost_dr8_mrs':'V/156','lamost_dr7':'V/153','rave_dr6':'III/283','apogee_dr17':'III/286','gaia_eso':'J/A+A/666/A120','galah_dr3':'J/MNRAS/506/150','elbadry_wide':'J/MNRAS/506/2269','wds':'B/wds/wds','tic82':'IV/39/tic82','gaiadr2':'I/345/gaia2','hipgaia_pma':'J/A+A/657/A7','kervella_pma':'J/A+A/623/A72','msc':'J/ApJS/235/6','2mass':'II/246/out'}
v=Vizier(columns=['**'],row_limit=50,timeout=120)
for k,cat in cats.items():
    r= 30*u.arcsec if k in('gaiadr3','elbadry_wide','wds','2mass','tic82') else 5*u.arcsec
    try:
        res=v.query_region(c,radius=r,catalog=cat)
    except Exception as e:
        print(k,cat,'ERR',e); continue
    if len(res)==0: print(k,cat,'no rows'); continue
    for i,t in enumerate(res):
        fn=f'_sources/vizier_{k}_{i}.csv'; t.write(fn,overwrite=True,format='ascii.csv')
        print(k,cat,t.meta.get('name',''),len(t),'rows ->',fn)
