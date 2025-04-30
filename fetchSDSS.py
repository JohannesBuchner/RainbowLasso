"""
Retrieves aperture fluxes from legacy survey via noirlab.
"""
import os
import requests_cache
requests_cache.install_cache('sdss_cache', allowable_methods=('GET', 'POST')) #, expire_after=3600*24*7)

import numpy as np
from astropy.table import Table
from astropy.coordinates import SkyCoord
import astropy.units as u
import pandas as pd
import tqdm

from joblib import Memory
from SciServer import Authentication, CasJobs

flagnames = [
    'CANONICAL_CENTER', 'BRIGHT', 'EDGE', 'BLENDED', 'CHILD', 'PEAKCENTER', 'NODEBLEND',
    'NOPROFILE', 'NOPETRO', 'MANYPETRO', 'NOPETRO_BIG', 'DEBLEND_TOO_MANY_PEAKS',
    'COSMIC_RAY', 'MANYR50', 'MANYR90', 'BAD_RADIAL', 'INCOMPLETE_PROFILE', 'INTERP',
    'SATURATED', 'NOTCHECKED', 'SUBTRACTED', 'NOSTOKES', 'BADSKY', 'PETROFAINT',
    'TOO_LARGE', 'DEBLENDED_AS_PSF', 'DEBLEND_PRUNED', 'ELLIPFAINT',
    'BINNED1', 'BINNED2', 'BINNED4', 'MOVED', 'DEBLENDED_AS_MOVING', 'NODEBLEND_MOVING',
    'TOO_FEW_DETECTIONS', 'BAD_MOVING_FIT', 'STATIONARY', 'PEAKS_TOO_CLOSE',
    'MEDIAN_CENTER', 'LOCAL_EDGE', 'BAD_COUNTS_ERROR', 'BAD_MOVING_FIT_CHILD',
    'DEBLEND_UNASSIGNED_FLUX', 'SATUR_CENTER', 'INTERP_CENTER', 'DEBLENDED_AT_EDGE',
    'DEBLEND_NOPEAK', 'PSF_FLUX_INTERP', 'TOO_FEW_GOOD_DETECTIONS', 'CENTER_OFF_AIMAGE',
    'DEBLEND_DEGENERATE', 'BRIGHTEST_GALAXY_CHILD', 'CANONICAL_BAND', 'AMOMENT_FAINT',
    'AMOMENT_SHIFT', 'AMOMENT_MAXITER', 'MAYBE_CR', 'MAYBE_EGHOST', 'NOTCHECKED_CENTER',
    'OBJECT2_HAS_SATUR_DN', 'OBJECT2_DEBLEND_PEEPHOLE', 'GROWN_MERGED', 'HAS_CENTER', 'RESERVED'
]
#bad_flags = ['NOPROFILE', 'BAD_RADIAL', 'SATURATED', 'SATUR_CENTER']

def good_flags(flags, flagmask):
    """Check photometry flags

    Parameters
    ----------
    flags: int
        flags set, see https://cas.sdss.org/dr7/en/help/browser/enum.asp?n=PhotoFlags

    Returns
    -------
    status: bool
        whether the photometry is reliable.
    """
    return (flags & flagmask) == 0

# Softening parameters b for each band
b_values = {
    'u': 1.4e-10,
    'g': 0.9e-10,
    'r': 1.2e-10,
    'i': 1.8e-10,
    'z': 7.4e-10,
}
# Zero-point flux in Jy
f0_jy = 3631

def sdss_mag_to_flux_mjy(mag, mag_err, band):
    """Convert SDSS asinh magnitude to flux in mJy.

    Parameters
    ----------
    mag: float
        asinh magnitude.
    band: str
        Photometric band ('u', 'g', 'r', 'i', 'z').

    Returns
    -------
    flux_mjy: float
        Flux in mJy.
    flux_err_mjy: float
        Flux error in mJy.
    """
    if band not in b_values:
        raise ValueError(f"Unknown band '{band}'. Must be one of {list(b_values.keys())}.")
    b = b_values[band]
    flux_jy = 2 * b * f0_jy * np.sinh((mag * np.log(10) / -2.5) - np.log(b))
    fluxerr_jy = 2 * b * f0_jy * np.abs(flux_jy * mag_err * np.log(10) / -2.5 / np.tanh((mag * np.log(10) / -2.5) - np.log(b)))
    return flux_jy * 1000, fluxerr_jy * 1000

mem = Memory('.', verbose=False)

@mem.cache
def executeQuery(**kwargs):
    return CasJobs.executeQuery(**kwargs)


@mem.cache
def fetch_one(ra, dec, bad_flags):
    flagmask = 0
    for flagname in bad_flags:
        flagmask |= 1 << flagnames.index(flagname)

    results = executeQuery(sql=f"""
    SELECT p.objID, p.ra, p.dec,
    p.psfMag_u, p.psfMag_g, p.psfMag_r, p.psfMag_i, p.psfMag_z,
    p.psfMagErr_u, p.psfMagErr_g, p.psfMagErr_r, p.psfMagErr_i, p.psfMagErr_z,
    r.profMean * 28.27 * 1e9 * 3.631e-3 as profMean, r.profErr * 28.27 * 1e9 * 3.631e-3 as profErr,
    r.band, 
    p.extinction_u, p.extinction_g, p.extinction_r, p.extinction_i, p.extinction_z,
    p.flags_u, p.flags_g, p.flags_r, p.flags_i, p.flags_z,
    p.expRad_u, p.expRad_g, p.expRad_r, p.expRad_i, p.expRad_z
    FROM dbo.PhotoProfile as r
    CROSS APPLY dbo.fGetNearestObjEq({ra}, {dec}, 0.0167) AS n
    LEFT JOIN PhotoObj AS p ON n.objid=p.objid 
    WHERE p.objID = r.objID
    AND r.bin = 4
    """, context="DR7", format="pandas")
    for i, band in enumerate('ugriz'):
        if len(results) > 0:
            if good_flags(results['flags_' + band][0], flagmask):
                extinction = results['extinction_' + band]
                flux_and_err = sdss_mag_to_flux_mjy(results['psfMag_' + band][0], results['psfMagErr_' + band][0], band)
                results['psfMag_' + band] = flux_and_err[0] - extinction
                results['psfMagErr_' + band] = flux_and_err[1]
                assert flux_and_err[1] >= 0
                profMean = float(results[results['band'] == i]['profMean'].iloc[0])
                profErr = float(results[results['band'] == i]['profErr'].iloc[0])
                results['aper_' + band] = profMean * 10**(0.4*extinction)
                results['aper_' + band + '_err'] = profErr * 10**(0.4*extinction)
                assert profMean >= 0, profMean
                assert profErr >= 0, profErr
            else:
                results['psfMag_' + band] = np.nan
                results['psfMagErr_' + band] = np.nan
                results['aper_' + band] = np.nan
                results['aper_' + band + '_err'] = np.nan
        del results['extinction_' + band]
    del results['profMean'], results['profErr']
    return results[:1]


def main(query_radius, input_table, output_table):
    t = Table.read(sys.argv[1]).filled()
    try:
        Authentication.login(*open(os.path.expanduser('~/.config/sciserver/login.txt')).readline().strip().split(':'))
    except FileNotFoundError:
        print("Could not log into sciserver, not filling in SDSS information")
        tmock = t[['id']]
        for i, band in enumerate('ugriz'):
            tmock['psfMag_' + band] = np.nan
            tmock['psfMagErr_' + band] = np.nan
            tmock['aper_' + band] = np.nan
            tmock['aper_' + band + '_err'] = np.nan
        tmock.write(output_table, overwrite=True)
        print(tmock)
        return

    elements = []
    # ra0,dec0,radius all in decimal degrees
    radius = 3 / 3600. * query_radius
    for row in tqdm.tqdm(t):
        result = fetch_one(
            ra=row['RA'],
            dec=row['DEC'],
            bad_flags = ['NOPROFILE', 'SATURATED', 'SATUR_CENTER']
        )
        result['id'] = row['id']
        a = SkyCoord(row['RA'], row['DEC'], unit='deg')
        b = SkyCoord(result['ra'], result['dec'], unit='deg')
        del result['ra'], result['dec'], result['objID']
        sep = a.separation(b)
        mask = sep < radius * u.deg
        if mask.any():
            elements.append(result[mask])

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
