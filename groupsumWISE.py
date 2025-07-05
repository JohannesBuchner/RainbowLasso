"""

"""
import numpy as np
from astropy.table import Table
import sys

t = Table.read(sys.argv[1])
t_ALLWISE = Table.read(sys.argv[2])
bands = ['WISE1', 'WISE2', 'WISE3', 'WISE4']
# 5 sigma sensitivity from external check https://wise2.ipac.caltech.edu/docs/release/allwise/expsup/sec2_3a.html
allsky_upper_limits_mJy = {'WISE1':68, 'WISE2':111, 'WISE3':860, 'WISE4':5700}

df = t.to_pandas()
df_allwise = t_ALLWISE.to_pandas()
df_allwise['id'] = df_allwise['id_in'].str.decode('ascii')
del df_allwise['id_in']
valid_mask = np.logical_or.reduce([
    df_allwise[f"{band}_err"] > 0 for band in bands
])
valid = df_allwise[valid_mask].copy()
print("aggregating...")
df_summed = valid.groupby('id').agg(
    **{band: (band, 'sum') for band in bands},
    **{f'{band}_err': (f'{band}_err', np.linalg.norm) for band in bands}
).reset_index()
print("merging and filling in NaN...")
df_result = df[['id']].copy().merge(df_summed, on='id', how='left')
for band in bands:
    bad = ~np.logical_and(df_result[band] > 0, df_result[f'{band}_err'] > 0)
    df_result.loc[bad, band] = allsky_upper_limits_mJy[band]
    df_result.loc[bad, f'{band}_err'] = np.nan

final_table = Table.from_pandas(df_result)
final_table.write(sys.argv[3], overwrite=True)
