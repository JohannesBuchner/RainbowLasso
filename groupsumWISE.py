"""
Retrieves aperture fluxes from legacy survey via noirlab.
"""
import numpy as np
from astropy.table import Table
import tqdm
import sys

t = Table.read(sys.argv[1])
t_ALLWISE = Table.read(sys.argv[2])
bands = ['WISE1', 'WISE2', 'WISE3', 'WISE4']

rows = []
for i in tqdm.tqdm(t['id']):
    row = [i]
    mask = t_ALLWISE['id_in'] == i
    if mask.any():
        # take sum over all nearby detections
        t_ALLWISE_sel = t_ALLWISE[mask]
        for band in bands:
            row.append(t_ALLWISE_sel[band].sum())
    else:
        # no entry means non-detection
        for band in bands:
            row.append(np.nan)
    rows.append(tuple(row))

data = Table(rows=rows, names=['id_in'] + bands)
data.write(sys.argv[3], overwrite=True)
