"""
Retrieves aperture fluxes from legacy survey via noirlab.
"""
import requests_cache
requests_cache.install_cache('demo_cache', allowable_methods=('GET', 'POST')) #, expire_after=3600*24*7)

import numpy as np
from astropy.table import Table
import pandas as pd
import tqdm
import sys
from dl import queryClient as qc

t = Table.read(sys.argv[1])

dtypes = dict(
  ls_id=int,
  ref_id=int,
  brickid=int,
  objid=int,
  maskbits=int,
  ref_epoch=int,
  gaia_phot_g_n_obs=int,
  release=str,
  fitbits=int,
  gaia_phot_bp_n_obs=int,
  gaia_phot_rp_n_obs=int,
  gaia_phot_variable_flag=int,
  gaia_duplicated_source=int,
  nobs_g=int,
  nobs_r=int,
  nobs_i=int,
  nobs_z=int,
  nobs_w1=int,
  nobs_w2=int,
  nobs_w3=int,
  nobs_w4=int,
  ngood_g=int,
  ngood_r=int,
  ngood_i=int,
  ngood_z=int,
  ref_cat=str,
  wise_coadd_id=str,
)

elements = []
for row in tqdm.tqdm(t):
    ra0, dec0, radius = row['RA'], row['DEC'], 0.00028 # ra0,dec0,radius all in decimal degrees
    query = """SELECT TOP 1 * FROM ls_dr10.tractor 
    INNER JOIN ls_dr10.apflux ON ls_dr10.tractor.ls_id = ls_dr10.apflux.ls_id 
    WHERE  q3c_radial_query(ra,dec,{:f},{:f},{:f})""".format(ra0,dec0,radius)
    rowdf = qc.query(query, fmt='pandas')
    if len(rowdf) == 0:
        continue
    for k, newtype in dtypes.items():
        if rowdf[k].dtype == 'object':
            # convert if we got a mixed type "object"
            rowdf[k] = rowdf[k].astype(str).astype(newtype)
    rowdf['id'] = pd.Series([row['id'] * len(rowdf)])
    elements.append(rowdf)

df = pd.concat(elements)
data = Table.from_pandas(df)

for c in data.colnames:
    if data.columns[c].dtype == np.dtype('O'):
        print("  stripping complex column", c)
        del data[c]

data.write(sys.argv[2], overwrite=True)
