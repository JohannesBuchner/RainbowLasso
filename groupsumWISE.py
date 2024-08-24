"""

"""
import numpy as np
from astropy.table import Table
import tqdm
from uncertainties import unumpy
import sys

t = Table.read(sys.argv[1])
t_ALLWISE = Table.read(sys.argv[2])
bands = ['WISE1', 'WISE2', 'WISE3', 'WISE4']
# 5 sigma sensitivity from external check https://wise2.ipac.caltech.edu/docs/release/allwise/expsup/sec2_3a.html
allsky_upper_limits_mJy = {'WISE1':68, 'WISE2':111, 'WISE3':860, 'WISE4':5700}

rows = []
for i in tqdm.tqdm(t['id']):
    row = [i]
    for band in bands:
        mask = np.logical_and(t_ALLWISE['id_in'] == i, t_ALLWISE[band + '_err'] > 0)
        if mask.any():
            # take sum over all nearby detections
            t_ALLWISE_sel = t_ALLWISE[mask]
            total = unumpy.uarray(t_ALLWISE_sel[band], t_ALLWISE_sel[band + '_err']).sum()
            row += [total.nominal_value, total.std_dev]
        else:
            # no entry means non-detection
            row += [allsky_upper_limits_mJy[band], np.nan]
    rows.append(tuple(row))

column_names = ['id']
for b in bands:
    column_names += [b, b + '_err']

data = Table(rows=rows, names=column_names)
data.write(sys.argv[3], overwrite=True)
