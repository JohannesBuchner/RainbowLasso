#!/bin/env python3
"""
Retrieves aperture fluxes from HSC survey via noirlab.
"""

import csv
import io
import json
import os
import sys
import time

import astropy.io.fits as pyfits
import requests
import tqdm
from astropy.table import Table, vstack

import joblib


class QueryError(Exception):
    pass


mem = joblib.Memory(".", verbose=False)

api_url = "https://hsc-release.mtk.nao.ac.jp/datasearch/api/catalog_jobs/"
release_version = "pdr3-citus-columnar"
skip_syntax_check = False
nomail = True
version = 20190514.1


def httpJsonPost(url, data, session=requests):
    """Submit to server.

    Parameters
    ----------
    url: str
        URL
    data: dict
        data to send
    session: object
        Request session to use

    Returns
    -------
    object
        whatever request.post returns
    """
    data["clientVersion"] = version
    return session.post(
        url, data=json.dumps(data), headers={"Content-Type": "application/json"}
    )


def submitJob(credential, sql, out_format):
    """Submit SQL query as a job to server.

    Parameters
    ----------
    credential: dict
        Login details
    sql: str
        SQL query
    out_format: str
        format to use for returned table, can by CSV, FITS, ...

    Returns
    -------
    info: dict
        information about submitted job
    """
    url = api_url + "submit"
    catalog_job = {
        "sql": sql,
        "out_format": out_format,
        "include_metainfo_to_body": True,
        "release_version": release_version,
    }
    postData = {
        "credential": credential,
        "catalog_job": catalog_job,
        "nomail": nomail,
        "skip_syntax_check": skip_syntax_check,
    }
    res = httpJsonPost(url, postData)
    assert res.status_code == 200, (
        "ERROR in query, response status code:",
        res.status_code,
        "message:",
        res.text,
    )
    return res.json()


def jobStatus(credential, job_id):
    """Check status of job running on the server.

    Parameters
    ----------
    credential: dict
        Login details
    job_id: str
        ID of job

    Returns
    -------
    info: dict
        information about the job
    """
    url = api_url + "status"
    postData = {"credential": credential, "id": job_id}
    res = httpJsonPost(url, postData)
    return res.json()


def jobCancel(credential, job_id):
    """Cancel job currently running on the server.

    Parameters
    ----------
    credential: dict
        Login details
    job_id: str
        ID of job
    """
    url = api_url + "cancel"
    postData = {"credential": credential, "id": job_id}
    httpJsonPost(url, postData)


verbose = os.environ.get("VERBOSE_DOWNLOAD", "1") == "1"
print_queries = os.environ.get("VERBOSE_DOWNLOAD", "1") == "2"


def blockUntilJobFinishes(credential, job_id):
    """Poll server until job is done.

    Parameters
    ----------
    credential: dict
        Login details
    job_id: str
        ID of job
    """
    max_interval = 5
    # seconds
    interval = 1
    # how long to wait at first
    last_status = "submitted"
    while True:
        if verbose:
            sys.stderr.write(
                f"job {job_id}: {last_status}, checking again in {interval:.1f}s ...\r"
            )
        time.sleep(interval)
        job = jobStatus(credential, job_id)
        if job["status"] == "error":
            raise QueryError("query error: " + job["error"])
        if job["status"] == "done":
            break
        if job["status"] != last_status and verbose:
            last_status = job["status"]
        interval = min(max_interval, 1.2 * interval)


def download(credential, job_id):
    """Download results of a completed job from server.

    Parameters
    ----------
    credential: dict
        Login details
    job_id: str
        ID of job

    Returns
    -------
    data: file
        binary data of the result
    """
    url = api_url + "download"
    postData = {"credential": credential, "id": job_id}
    res = httpJsonPost(url, postData)
    return io.BytesIO(res.content)


def deleteJob(credential, job_id):
    """Delete the job and its results from server.

    Parameters
    ----------
    credential: dict
        Login details
    job_id: str
        ID of job
    """
    url = api_url + "delete"
    postData = {"credential": credential, "id": job_id}
    httpJsonPost(url, postData)

try:
    with open(os.path.expanduser("~/.config/hsc-password")) as fpass:
        credential = {
            "account_name": fpass.readline().strip(),
            "password": fpass.readline().strip(),
        }
except IOError:
    print(
        "create '~/.config/hsc-password' with two lines: user name and password, from https://hsc-release.mtk.nao.ac.jp/datasearch/new_user/new"
    )
    sys.exit(1)


@mem.cache
def fetchTable(query):
    """Fetch results to a query from server.

    Parameters
    ----------
    sql: str
        SQL query

    Returns
    -------
    ttmp: astropy.Table
        Table of results
    """
    job = submitJob(credential, query, out_format)
    blockUntilJobFinishes(credential, job["id"])
    ttmp = Table.read(download(credential, job["id"]))
    deleteJob(credential, job["id"])
    return ttmp


chunksize = 20

columns = [
    "_inputcount_value|Number of images contributing at center, not including anyclipping",
    "_inputcount_flag|Set for any fatal failure",
    "_localbackground_flag|General Failure Flag",
    "_pixelflags|General failure flag, set if anything went wrong",
    "_pixelflags_bad|Bad pixel in the Source footprint",
    "_pixelflags_cr|Cosmic ray in the Source footprint",
    "_pixelflags_saturatedcenter|Saturated pixel in the Source footprint",
    "_pixelflags_saturated|Saturated pixel in the Source center",
    "_pixelflags_edge|Source is outside usable exposure region (masked EDGE or NO_DATA)",
    "_extendedness_value|Set to 1 for extended sources, 0 for point sources.",
    "_extendedness_flag|Set to 1 for any fatal failure.",
    "_kronflux_psf_radius|Radius of PSF",
    "_kronflux_flag_bad_shape_no_psf|bad shape and no PSF",
    "_sdssshape_psf_shape11|adaptive moments of the PSF model at the object position [arcsec^2]",
    "_sdssshape_psf_shape12|adaptive moments of the PSF model at the object position [arcsec^2]",
    "_sdssshape_psf_shape22|adaptive moments of the PSF model at the object position [arcsec^2]",
    "_psfflux_flux|instFlux derived from linear least-squares fit of PSF model",
    "_psfflux_fluxerr|1-sigma instFlux uncertainty",
    "_psfflux_flag|General Failure Flag",
    "_psfflux_flag_nogoodpixels|not enough non-rejected pixels in data to attempt the fit",
    "_psfflux_flag_edge|object was too close to the edge of the image to use the full PSF model",
    "_psfflux_flag_apcorr|set if unable to aperture correct base_PsfFlux",
    "_psfflux_flag_badcentroid|whether the reference centroid is marked as bad",
    "_apertureflux_flag_badcentroid|whether the reference centroid is marked as bad",
    "_sdsscentroid_flag|General Failure Flag",
]
r_columns = [
    "_apertureflux_%d_flux|aperture flux",
    "_apertureflux_%d_fluxerr|1-sigma instFlux uncertainty",
    "_apertureflux_%d_flag|General Failure Flag",
    "_apertureflux_%d_flag_aperturetruncated|General Failure Flag",
    #'_apertureflux_%d_flag_sinccoeffstruncated|full sinc coefficient image did not fit within measurement image',
]
bands = "grizy"
all_column_names = ["ra", "dec", "detect_ispatchinner"]
for band in bands:
    all_column_names.append("a_" + band)
    all_column_names += [band + c.split("|")[0] for c in columns]
    for radius in 10, 15, 20, 30, 40, 57, 84, 118:
        all_column_names += [band + (c.split("|")[0] % radius) for c in r_columns]
out_format = "fits"
t = Table.read(sys.argv[1])

id_dtype = "int" if t["id"].dtype == int else "text"
query_radius = os.environ.get("QUERY_RADIUS", "1")

elements = []
values = []
for k, row in enumerate(t):
    i = row["id"]
    ra0, dec0 = row["RA"], row["DEC"]  # in decimal degrees
    # HSC-WIDE overlapping with XMM-LSS
    in_XMMLSS = ra0 > 325 or ra0 < 45 and -8 < dec0 < 9
    # HSC-WIDE overlapping with COSMOS
    in_COSMOS = 120 < ra0 < 230 and -3.5 < dec0 < 6.5
    in_Hectomap = 195 < ra0 < 225 and 41 < dec0 < 46
    in_AEGIS = 210 < ra0 < 218 and 51 < dec0 < 54

    if in_XMMLSS or in_COSMOS or in_Hectomap or in_AEGIS:
        values.append(
            f"('{i}'::{id_dtype},'{row['RA']:.16e}'::double precision,'{row['DEC']:.16e}'::double precision)"
        )

    if len(values) > chunksize or k == len(t) - 1:
        id_query = (
            'WITH user_catalog("user.myid","user.RA","user.DEC") AS (VALUES'
            + ",".join(values)
            + """),
    match AS (
        SELECT
            object_id
        FROM
            user_catalog
        JOIN pdr3_wide.forced ON coneSearch(coord, "user.RA", "user.DEC", """
            + query_radius
            + """)
    )
SELECT *
FROM match
"""
        )

        if print_queries:
            print(id_query)

        object_ids = fetchTable(id_query)["object_id"]
        if verbose:
            print("object ids:", [int(i) for i in object_ids])

        query = (
            "SELECT "
            + ",".join(all_column_names)
            + """
FROM pdr3_wide.forced
LEFT JOIN pdr3_wide.forced2 USING (object_id)
LEFT JOIN pdr3_wide.forced3 USING (object_id)
WHERE object_id in ("""
            + ",".join(("%d" % i for i in object_ids))
            + """)
"""
        )
        ttmp = fetchTable(query)
        if len(ttmp) > 0:
            for col in ttmp.colnames:
                if col.endswith("_isnull"):
                    del ttmp[col]
        elements.append(ttmp)
        values = []

print("\nstoring results...")
data = vstack(elements)
data.write(sys.argv[2], overwrite=True)

fout = pyfits.open(sys.argv[2])
for band in bands:
    for c in columns:
        colname, coldesc = c.split("|")
        i = fout[1].data.columns.names.index(band + colname)
        fout[1].header["TCOMM%d" % i] = coldesc
    for radius in 10, 15, 20, 30, 40, 57, 84, 118:
        for c in r_columns:
            colname, coldesc = c.split("|")
            i = fout[1].data.columns.names.index(band + colname % radius)
            fout[1].header["TCOMM%d" % i] = coldesc
fout.writeto(sys.argv[2], overwrite=True)
