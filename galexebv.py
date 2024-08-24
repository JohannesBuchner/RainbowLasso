import numpy as np
from astropy.units import deg
from astropy.coordinates import SkyCoord
from astropy.table import Table
from dustmaps.sfd import SFDQuery
import tqdm

sfd = SFDQuery()

table = Table.read('galex_ais_ctrs.fits')
coords = SkyCoord(ra=table['RAfdeg'], dec=table['DEfdeg'], unit='deg', frame='icrs')

ebv = sfd(coords)

"""
N = 1
for i in tqdm.trange(1):
    offsetx, offsety = np.random.uniform(-1, 1, size=2)
    r = np.sqrt(offsetx**2 + offsety**2)
    if r > 1: continue
    theta = np.arctan2(offsety, offsetx)
    offset_coords = coords.directional_offset_by(r * 1.1 / 2.0 * deg, theta * deg)
    newebv = sfd(offset_coords)
    ebv[newebv < ebv] = newebv[newebv < ebv]
    N += 1

ebv /= N
"""

table['EBV_min'] = ebv
table.write('galex_ais_ctrs_ebv.fits', overwrite=True)
