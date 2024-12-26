#!/bin/env python3
"""
Retrieves aperture fluxes from HSC survey via noirlab.
"""
import tqdm
import requests
import json
import os
import time
import sys
import csv
from astropy.table import Table, vstack
import astropy.io.fits as pyfits
import io
import joblib

mem = joblib.Memory('.', verbose=False)

api_url = 'https://hsc-release.mtk.nao.ac.jp/datasearch/api/catalog_jobs/'
release_version = 'pdr3-citus-columnar'
skip_syntax_check = False
nomail = True
version = 20190514.1

#cached_session = requests_cache.CachedSession('demo_cache', allowable_methods=('GET', 'POST')) #, expire_after=3600*24*7)

class QueryError(Exception):
    pass

def httpJsonPost(url, data, session=requests):
    data['clientVersion'] = version
    return session.post(url, data=json.dumps(data), headers={'Content-Type': 'application/json'})

def preview(credential, sql, out):
    url = api_url + 'preview'
    catalog_job = {
        'sql'             : sql,
        'release_version' : release_version,
    }
    postData = {'credential': credential, 'catalog_job': catalog_job}
    res = httpJsonPost(url, postData)
    result = res.json()

    writer = csv.writer(out)
    for row in result['result']['rows']:
        writer.writerow(row)

    if result['result']['count'] > len(result['result']['rows']):
        raise QueryError('only top %d records are displayed !' % len(result['result']['rows']))

def submitJob(credential, sql, out_format):
    url = api_url + 'submit'
    catalog_job = {
        'sql'                     : sql,
        'out_format'              : out_format,
        'include_metainfo_to_body': True,
        'release_version'         : release_version,
    }
    postData = {'credential': credential, 'catalog_job': catalog_job, 'nomail': nomail, 'skip_syntax_check': skip_syntax_check}
    res = httpJsonPost(url, postData)
    assert res.status_code == 200, ('ERROR in query, response status code:', res.status_code, 'message:', res.text)
    return res.json()


def jobStatus(credential, job_id):
    url = api_url + 'status'
    postData = {'credential': credential, 'id': job_id}
    res = httpJsonPost(url, postData)
    return res.json()


def jobCancel(credential, job_id):
    url = api_url + 'cancel'
    postData = {'credential': credential, 'id': job_id}
    httpJsonPost(url, postData)

verbose = False

def blockUntilJobFinishes(credential, job_id):
    max_interval = 5  # seconds
    interval = 1  # how long to wait at first
    last_status = 'submitted'
    while True:
        if verbose:
            sys.stderr.write(f"job {job_id}: {last_status}, checking again in {interval:.1f}s ...\r")
        time.sleep(interval)
        job = jobStatus(credential, job_id)
        if job['status'] == 'error':
            raise QueryError('query error: ' + job['error'])
        if job['status'] == 'done':
            break
        if job['status'] != last_status and verbose:
            last_status = job['status']
        interval = min(max_interval, 1.2 * interval)

def download(credential, job_id):
    url = api_url + 'download'
    postData = {'credential': credential, 'id': job_id}
    res = httpJsonPost(url, postData)
    return io.BytesIO(res.content)

def deleteJob(credential, job_id):
    url = api_url + 'delete'
    postData = {'credential': credential, 'id': job_id}
    httpJsonPost(url, postData)


try:
    with open(os.path.expanduser('~/.config/hsc-password')) as fpass:
        credential = {
            'account_name': fpass.readline().strip(),
            'password': fpass.readline().strip()
        }
except IOError:
    print("create '~/.config/hsc-password' with two lines: user name and password, from https://hsc-release.mtk.nao.ac.jp/datasearch/new_user/new")
    sys.exit(1)

@mem.cache
def fetchTable(query):
    job = submitJob(credential, query, out_format)
    blockUntilJobFinishes(credential, job['id'])
    ttmp = Table.read(download(credential, job['id']))
    deleteJob(credential, job['id'])
    return ttmp

columns = [
'_inputcount_value|Number of images contributing at center, not including anyclipping',
'_inputcount_flag|Set for any fatal failure',
'_localbackground_flag|General Failure Flag',
'_pixelflags|General failure flag, set if anything went wrong',
'_pixelflags_bad|Bad pixel in the Source footprint',
'_pixelflags_cr|Cosmic ray in the Source footprint',
'_pixelflags_saturatedcenter|Saturated pixel in the Source footprint',
'_pixelflags_saturated|Saturated pixel in the Source center',
'_pixelflags_edge|Source is outside usable exposure region (masked EDGE or NO_DATA)',
'_extendedness_value|Set to 1 for extended sources, 0 for point sources.',
'_extendedness_flag|Set to 1 for any fatal failure.',
'_kronflux_psf_radius|Radius of PSF',
'_kronflux_flag_bad_shape_no_psf|bad shape and no PSF',
'_sdssshape_psf_shape11|adaptive moments of the PSF model at the object position [arcsec^2]',
'_sdssshape_psf_shape12|adaptive moments of the PSF model at the object position [arcsec^2]',
'_sdssshape_psf_shape22|adaptive moments of the PSF model at the object position [arcsec^2]',
'_psfflux_flux|instFlux derived from linear least-squares fit of PSF model',
'_psfflux_fluxerr|1-sigma instFlux uncertainty',
'_psfflux_flag|General Failure Flag', 
'_psfflux_flag_nogoodpixels|not enough non-rejected pixels in data to attempt the fit',
'_psfflux_flag_edge|object was too close to the edge of the image to use the full PSF model',
'_psfflux_flag_apcorr|set if unable to aperture correct base_PsfFlux',
'_psfflux_flag_badcentroid|whether the reference centroid is marked as bad',
'_apertureflux_flag_badcentroid|whether the reference centroid is marked as bad',
'_sdsscentroid_flag|General Failure Flag',
]
r_columns = [
'_apertureflux_%d_flux|aperture flux',
'_apertureflux_%d_fluxerr|1-sigma instFlux uncertainty',
'_apertureflux_%d_flag|General Failure Flag',
'_apertureflux_%d_flag_aperturetruncated|General Failure Flag',
#'_apertureflux_%d_flag_sinccoeffstruncated|full sinc coefficient image did not fit within measurement image',
]
bands = 'grizy'
all_column_names = ['ra', 'dec','detect_ispatchinner']
for band in bands:
    all_column_names.append('a_' + band)
    all_column_names += [band + c.split("|")[0] for c in columns]
    for radius in 10, 15, 20, 30, 40, 57, 84, 118:
        all_column_names += [band + (c.split("|")[0] % radius) for c in r_columns]
out_format = 'fits'
t = Table.read(sys.argv[1])

id_dtype = 'int' if t['id'].dtype == int else 'text'
values = []
for row in t:
    i = row['id']
    ra0, dec0 = row['RA'], row['DEC'] # in decimal degrees
    if ra0 > 325 or ra0 < 45 and -8 < dec0 < 9:
        pass # HSC-WIDE overlapping with XMM-LSS
    elif 120 < ra0 < 230 and -3.5 < dec0 < 6.5:
        pass # HSC-WIDE overlapping with COSMOS
    elif 195 < ra0 < 225 and 41 < dec0 < 46:
        pass # Hectomap
    elif 210 < ra0 < 218 and 51 < dec0 < 54:
        pass # AEGIS
    else:
        continue # skip outside
    values.append(f"('{i}'::{id_dtype},'{row['RA']:.16e}'::double precision,'{row['DEC']:.16e}'::double precision)")

#new_column_names = [
#    ('forced2' if '_psfflux_' in colname else ('forced3' if '_apertureflux_' in colname else 'forced')) + '.' + colname + ' as ' + colname for colname in all_column_names]
id_query = 'WITH user_catalog("user.myid","user.RA","user.DEC") AS (VALUES' + ','.join(values) + """),
    match AS (
        SELECT
            object_id
        FROM
            user_catalog
        JOIN pdr3_wide.forced ON coneSearch(coord, "user.RA", "user.DEC", 1)
    )
SELECT *
FROM match
"""

if verbose:
    print(id_query)

object_ids = fetchTable(id_query)['object_id']
if verbose:
    print('object ids:', list(object_ids))

chunksize = 60
elements = []
for i in tqdm.trange(0, len(object_ids), chunksize):
    query = 'SELECT ' + ','.join(all_column_names) + """
FROM pdr3_wide.forced
LEFT JOIN pdr3_wide.forced2 USING (object_id)
LEFT JOIN pdr3_wide.forced3 USING (object_id)
WHERE object_id in (""" + ','.join(('%d' % i for i in object_ids[i:i+chunksize])) + """)
"""
    ttmp = fetchTable(query)
    if len(ttmp) > 0:
        for col in ttmp.colnames:
            if col.endswith('_isnull'):
                del ttmp[col]
    elements.append(ttmp)

print("\nstoring results...")
data = vstack(elements)
data.write(sys.argv[2], overwrite=True)

fout = pyfits.open(sys.argv[2])
for band in bands:
    for c in columns:
        colname, coldesc = c.split("|")
        i = fout[1].data.columns.names.index(band + colname)
        fout[1].header['TCOMM%d' % i] = coldesc
    for radius in 10, 15, 20, 30, 40, 57, 84, 118:
        for c in r_columns:
            colname, coldesc = c.split("|")
            i = fout[1].data.columns.names.index(band + colname % radius)
            fout[1].header['TCOMM%d' % i] = coldesc
fout.writeto(sys.argv[2], overwrite=True)

