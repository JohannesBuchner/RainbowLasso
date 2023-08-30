"""
Retrieves aperture fluxes from legacy survey via noirlab.
"""
import requests_cache
requests_cache.install_cache('demo_cache', allowable_methods=('GET', 'POST'), expire_after=3600*24*7)

import numpy as np
from astropy.table import Table
import pandas as pd
import tqdm
import sys
from dl import queryClient as qc

t = Table.read(sys.argv[1])

elements = []
for row in tqdm.tqdm(t):
    ra0, dec0, radius = row['RA'], row['DEC'], 1./60/60 # ra0,dec0,radius all in decimal degrees
    query = """SELECT TOP 1 * FROM ukidss_dr11plus.lassource 
    WHERE  q3c_radial_query(ra,dec,{:f},{:f},{:f})""".format(ra0,dec0,radius)
    rowdf = qc.query(query, fmt='pandas')
    rowdf['id'] = pd.Series([row['id']])
    elements.append(rowdf)

data = Table.from_pandas(pd.concat(elements))

for c in data.colnames:
    if data.columns[c].dtype == np.dtype('O'):
        print("  stripping complex column", c)
        del data[c]

data.write(sys.argv[2], overwrite=True)
