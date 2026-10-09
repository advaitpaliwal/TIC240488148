import pandas as pd, numpy as np
from astropy.table import Table
from astroquery.gaia import Gaia
t=pd.read_csv('pool_all.csv'); t=t[t.ns>=3]
up=Table.from_pandas(t[['ID','ra','dec']].rename(columns={'ID':'tic'}))
q="""SELECT u.tic, g.source_id, g.phot_g_mean_mag, g.ruwe, g.ra, g.dec,
 DISTANCE(POINT(u.ra,u.dec),POINT(g.ra,g.dec))*3600 AS sep
 FROM tap_upload.u AS u JOIN gaiadr3.gaia_source AS g
 ON 1=CONTAINS(POINT(g.ra,g.dec),CIRCLE(u.ra,u.dec,30./3600))"""
j=Gaia.launch_job_async(q,upload_resource=up,upload_table_name='u')
r=j.get_results().to_pandas(); r.to_csv('_sources/gaia_30as.csv',index=False); print(len(r))
