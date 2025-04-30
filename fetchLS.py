"""
Retrieves aperture fluxes from legacy survey via noirlab.
"""
import os
import requests_cache
requests_cache.install_cache('demo_cache', allowable_methods=('GET', 'POST')) #, expire_after=3600*24*7)

import numpy as np
from astropy.table import Table
import pandas as pd
import tqdm
from dl import queryClient as qc

from joblib import Memory


# set data type to avoid conversion issues where no data are returned
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

# delete some colmns which we are surely not going to use
cols_to_delete = "htm9 nest4096 ring256 ref_cat wise_coadd_id g_r r_i r_z i_z z_w1 w1_w2 w2_w3 w3_w4 wise_x wise_y glon glat elon elat bx by".split()
col_stems_to_delete = "apflux_blobresid_ apflux_resid_ resid psfdepth_ blob_ nea_ galdepth_ rchisq_ fiberflux".split()

mem = Memory('.', verbose=False)

@mem.cache
def fetch_chunk(ids, ras, decs, radius):
    elements = []
    for id0, ra0, dec0 in zip(tqdm.tqdm(ids), ras, decs):
        if radius * 3600 < 1.01:  # find object
            query = """SELECT TOP 1 * FROM ls_dr10.tractor 
    INNER JOIN ls_dr10.apflux ON ls_dr10.tractor.ls_id = ls_dr10.apflux.ls_id 
    WHERE  q3c_radial_query(ra,dec,{:f},{:f},{:f})""".format(ra0,dec0,radius)
        else: # find all neighbours
            query = """SELECT FROM ls_dr10.tractor 
    INNER JOIN ls_dr10.apflux ON ls_dr10.tractor.ls_id = ls_dr10.apflux.ls_id 
    WHERE  q3c_radial_query(ra,dec,{:f},{:f},{:f})""".format(ra0,dec0,radius)

        rowdf = qc.query(query, fmt='pandas')

        if len(rowdf) == 0:
            continue

        for k, newtype in dtypes.items():
            if rowdf[k].dtype == 'object':
                # convert if we got a mixed type "object"
                rowdf[k] = rowdf[k].astype(str).astype(newtype)

        rowdf['id'] = pd.Series([id0] * len(rowdf))

        for col_to_delete in cols_to_delete:
            del rowdf[col_to_delete]
        for col_stem_to_delete in col_stems_to_delete:
            for col in rowdf.columns:
                if col.startswith(col_stem_to_delete):
                    del rowdf[col]
        elements.append(rowdf)
    if len(elements) > 0:
        print('fetch_chunk:', len(elements))
        return [pd.concat(elements)]
    else:
        return []


def main(query_radius, input_table, output_table):
    t = Table.read(sys.argv[1]).filled()

    chunklengths = 10000
    elements = []
    # ra0,dec0,radius all in decimal degrees
    radius = 0.00028 * query_radius
    for i in range(0, len(t), chunklengths):
        print(f"chunk {i}../{len(t)}...")
        elements += fetch_chunk(
            ids=t['id'][i : i + chunklengths],
            ras=t['RA'][i : i + chunklengths],
            decs=t['DEC'][i : i + chunklengths],
            radius=radius,
        )
    print("concatenating ...", len(elements))
    if len(elements) == 1:
        df = elements[0]
    else:
        df = pd.concat(elements)
    print("converting ...")
    data = Table.from_pandas(df)
    del df

    for c in data.colnames:
        if data.columns[c].dtype == np.dtype('O'):
            print("  stripping complex column", c)
            del data[c]

    print("writing ...")
    data.write(output_table, overwrite=True)

if __name__ == '__main__':
    import sys
    main(
        query_radius = float(os.environ.get('QUERY_RADIUS', '1')),
        input_table=sys.argv[1], output_table=sys.argv[2])
