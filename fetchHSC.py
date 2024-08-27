#!/bin/env python3
"""
Retrieves aperture fluxes from HSC survey via noirlab.
"""
import requests_cache
requests_cache.install_cache('demo_cache', allowable_methods=('GET', 'POST')) #, expire_after=3600*24*7)
import requests
import json
#import urllib.request, urllib.error, urllib.parse
#import astropy.io.fits
import time
import sys
import csv
from astropy.table import Table, vstack
import tqdm
import io

api_url = 'https://hsc-release.mtk.nao.ac.jp/datasearch/api/catalog_jobs/'
release_version = 'pdr3'
skip_syntax_check = False
nomail = False
version = 20190514.1

class QueryError(Exception):
    pass

"""
def httpPost(url, postData, headers):
    req = urllib.request.Request(url, postData.encode('utf-8'), headers)
    res = urllib.request.urlopen(req)
    return res

def httpJsonPost(url, data):
    data['clientVersion'] = version
    postData = json.dumps(data)
    return httpPost(url, postData, {'Content-type': 'application/json'})
"""

def httpPost(url, data, headers):
    response = requests.post(url, data=json.dumps(data), headers=headers)
    return response

def httpJsonPost(url, data):
    data['clientVersion'] = version
    response = httpPost(url, data, {'Content-Type': 'application/json'})
    return response

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


def blockUntilJobFinishes(credential, job_id):
    max_interval = 5 * 60 # sec.
    interval = 1
    while True:
        time.sleep(interval)
        job = jobStatus(credential, job_id)
        if job['status'] == 'error':
            raise QueryError('query error: ' + job['error'])
        if job['status'] == 'done':
            break
        interval *= 2
        if interval > max_interval:
            interval = max_interval

def download(credential, job_id):
    url = api_url + 'download'
    postData = {'credential': credential, 'id': job_id}
    res = httpJsonPost(url, postData)
    return io.BytesIO(res.content)
    #bufSize = 64 * 1<<10 # 64k
    #while True:
    #    buf = res.read(bufSize)
    #    out.write(buf)
    #    if len(buf) < bufSize:
    #        break

def deleteJob(credential, job_id):
    url = api_url + 'delete'
    postData = {'credential': credential, 'id': job_id}
    httpJsonPost(url, postData)



with open('/home/user/.config/hsc-password') as fpass:
    credential = {
        'account_name': fpass.readline().strip(),
        'password': fpass.readline().strip()
    }

columns = [
'_extendedness_value|Set to 1 for extended sources, 0 for point sources.',
'_psfflux_flux|instFlux derived from linear least-squares fit of PSF model',
'_psfflux_fluxerr|1-sigma instFlux uncertainty',
'_psfflux_flag|General Failure Flag', 
'_psfflux_flag_nogoodpixels|not enough non-rejected pixels in data to attempt the fit',
'_psfflux_flag_edge|object was too close to the edge of the image to use the full PSF model',
'_psfflux_flag_apcorr|set if unable to aperture correct base_PsfFlux',
'_psfflux_flag_badcentroid|whether the reference centroid is marked as bad',
'_apertureflux_flag_badcentroid|whether the reference centroid is marked as bad',
]
r_columns = [
'_apertureflux_%d_flux|aperture flux',
'_apertureflux_%d_fluxerr|1-sigma instFlux uncertainty',
'_apertureflux_%d_flag|General Failure Flag',
'_apertureflux_%d_flag_aperturetruncated|General Failure Flag',
#'_apertureflux_%d_flag_sinccoeffstruncated|full sinc coefficient image did not fit within measurement image',
]
bands = 'grizy'
#bands = 'gz'
all_column_names = ['ra', 'dec']
for band in bands:
    all_column_names.append('a_' + band)
    all_column_names += [band + c.split("|")[0] for c in columns]
    for radius in 10, 15, 20, 30, 40, 57, 84, 118:
        all_column_names += [band + (c.split("|")[0] % radius) for c in r_columns]
#print("\n".join(all_column_names))
out_format = 'fits'
t = Table.read(sys.argv[1])

elements = []
for row in tqdm.tqdm(t):
    ra0, dec0 = row['RA'], row['DEC'] # in decimal degrees
    radius = 1 # in arcsec
    query = 'SELECT ' + ','.join(all_column_names) + """
FROM pdr3_wide.forced
LEFT JOIN pdr3_wide.forced2 USING (object_id)
LEFT JOIN pdr3_wide.forced3 USING (object_id)
WHERE  coneSearch(coord,{:f},{:f},{:d})
LIMIT 1
  """.format(ra0,dec0,radius)
    print(query)

    job = submitJob(credential, query, out_format)
    blockUntilJobFinishes(credential, job['id'])
    #with open('tmp.fits', 'wb') as ftmp:
    #    download(credential, job['id'], ftmp)
    ttmp = Table.read(download(credential, job['id']))
    for col in ttmp.colnames:
        if col.endswith('_isnull'):
            del ttmp[col]
    ttmp['id'] = row['id']
    deleteJob(credential, job['id'])
    elements.append(ttmp)

print("storing...")
data = vstack(elements)
# add comments
data.write(sys.argv[2], overwrite=True)
