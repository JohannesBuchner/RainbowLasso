import numpy as np
import sys
from astropy.table import Table
from astropy.coordinates import SkyCoord
from astropy import units as u
from uncertainties import unumpy
import scipy.stats
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from PIL import Image
from requests_ratelimiter import LimiterMixin
from requests_cache import CacheMixin
from requests import Session
import astropy.io.fits as pyfits
from io import BytesIO
import os
from photutils.aperture import aperture_photometry, CircularAperture

print_urls = os.environ.get('PRINT_URLS', '0') == '1'

def image_viz_transform(img):
    with np.errstate(divide='ignore', invalid='ignore'):
        return np.log10(img)
    #return np.log10(np.clip(img, np.nanmax(img) / 1000., None))

def get_image_center(input_img):
    img = np.array(input_img)
    nx, ny = img.shape[0] - 1, img.shape[1] - 1
    x0, y0 = nx / 2, ny / 2
    #y0, x0 = np.unravel_index(img.argmax(), img.shape)
    # print(img.shape, x0, y0)
    return img, (x0, y0)

def get_aperture_curve_from_image(img, center, radii, ax=None, **circle_kwargs):
    x0, y0 = center
    #y0, x0 = np.unravel_index(img.argmax(), img.shape)
    # print(img.shape, x0, y0)
    positions = [center]
    values = []
    #last_apphot = 0.0
    for r in radii:
        phot_table = aperture_photometry(img, CircularAperture(positions, r=r))
        #values.append(float((phot_table['aperture_sum'] - last_apphot) / (2 * 3.14 * r**2)))
        values.append(float(phot_table['aperture_sum'][0]))
        #last_apphot = phot_table['aperture_sum']
        if ax is not None:
            ax.add_artist(Circle((x0, y0), r, **circle_kwargs))
    # print(values)
    return np.array(values)


browser_headers = {
    'User-Agent': 'Mozilla/5.0 (X11; Ubuntu; Linux x86_64; rv:109.0) Gecko/20100101 Firefox/112.0',
    'Accept': 'text/html,application/xhtml+xml,application/xml;q=0.9,image/avif,image/webp,*/*;q=0.8',
    'Accept-Language': 'en-US',
    'Accept-Encoding': 'gzip, deflate, br',
    'DNT': '1',
    'Connection': 'keep-alive',
}

class CachedLimiterSession(CacheMixin, LimiterMixin, Session):
    """Session class with caching and rate-limiting behavior.
    Accepts keyword arguments for both LimiterSession and CachedSession.
    """

my_session = CachedLimiterSession('demo_cache', allowable_methods=('GET', 'POST'), expire_after=3600 * 24 * 7, per_second=0.2)

def download_image(ra, dec, band, layer, verbose=print_urls):
    url = f'https://www.legacysurvey.org/viewer/jpeg-cutout?ra={ra}&dec={dec}&zoom=14&layer={layer}&bands={band}'
    if verbose:
        print(url)
    response = my_session.get(url)
    if verbose:
        print(url, response.status_code)
    image = Image.open(BytesIO(response.content))
    return image

from requests.auth import HTTPBasicAuth
try:
    with open(os.path.expanduser('~/.config/hsc-password')) as fpass:
        credential = {
            'account_name': fpass.readline().strip(),
            'password': fpass.readline().strip()
        }
        auth = HTTPBasicAuth(credential['account_name'], credential['password'])
except IOError:
    print("create '~/.config/hsc-password' with two lines: user name and password, from https://hsc-release.mtk.nao.ac.jp/datasearch/new_user/new")
    sys.exit(1)

def download_psf_image(ra, dec, band, layer, verbose=print_urls):
    url = f'https://hsc-release.mtk.nao.ac.jp/psf/pdr3/cgi/getpsf?ra={ra}&dec={dec}&filter={band}&rerun=pdr3_wide&tract=&patch=&centered=true&type=coadd'
    if verbose:
        print(url)
    response = my_session.get(url, auth=auth)
    if verbose:
        print(url, response.status_code)
    f = pyfits.open(BytesIO(response.content))[0]
    center_header = x0, y0 = (-(f.header['CRVAL1A']), -(f.header['CRVAL2A']))
    image = np.array(f.data[::-1,:])
    nx, ny = image.shape[0], image.shape[1]
    y0, x0 = np.unravel_index(image.argmax(), image.shape)
    center = x0, y0
    #print(nx, ny, x0, y0, (nx - 1)/2 - x0, (ny - 1)/2 - y0)
    if (nx - 1)/2 == x0 - 1:
        image = image[:-1,:]
    elif (nx - 1)/2 == x0 + 1:
        image = image[1:,:]
    if (ny - 1)/2 == y0 - 1:
        image = image[:,:-1]
    elif (ny - 1)/2 == y0 + 1:
        image = image[:,1:]
    #print('PSF:', center, image.shape, center_header)
    #assert center == (-20., -19.)
    #image = image[:-2, 2:]
    return image

def download_fits_image(ra, dec, band, layer, verbose=print_urls):
    url = f'https://hsc-release.mtk.nao.ac.jp/das_cutout/pdr3/cgi-bin/cutout?ra={ra}&dec={dec}&sw=0.001&sh=0.001&type=coadd&image=on&filter=HSC-{band.upper()}&tract=&rerun=pdr3_wide'
    if verbose:
        print(url)
    response = my_session.get(url, auth=auth)
    if verbose:
        print(url, response.status_code)
    f = pyfits.open(BytesIO(response.content))[1]
    image = f.data[::-1,:]
    return image

def plot_image(ax, pil_img):
    img = np.array(pil_img)
    nx, ny = img.shape[0] - 1, img.shape[1] - 1
    ax.imshow(img)
    ax.set_yticks([])
    ax.set_xticks([])
    ax.plot([nx / 2 + 10, nx / 2 + 4], [ny / 2, ny / 2], color='white', lw=0.4)
    ax.plot([nx / 2 - 4, nx / 2 - 10], [ny / 2, ny / 2], color='white', lw=0.4)
    ax.plot([nx / 2, nx / 2], [ny / 2 + 10, ny / 2 + 4], color='white', lw=0.4)
    ax.plot([nx / 2, nx / 2], [ny / 2 - 4, ny / 2 - 10], color='white', lw=0.4)

def monotonic_floodfill(array, start_pos):
    """
    Adjusts a 2D array so that the values monotonically decrease from a given start position.
    
    Arguments
    ---------
    array: np.ndarray
        The input 2D array of pixel values.
    start_pos: tuple
        The starting position (row, col) as a tuple.

    Returns
    -------
    np.ndarray:
        The modified 2D array with monotonically decreasing values.
    """
    rows, cols = array.shape
    start_row, start_col = int(start_pos[0]), int(start_pos[1])
    
    if not (0 <= start_row < rows and 0 <= start_col < cols):
        raise ValueError("Start position is out of bounds.")
    
    # Initialize BFS structures
    visited = np.zeros((rows, cols), dtype=bool)
    queue = [(start_row, start_col)]
    visited[start_row, start_col] = True

    # Directions: Up, Down, Left, Right
    directions = [(-1, 0), (1, 0), (0, -1), (0, 1)]

    while queue:
        r, c = queue.pop(0)
        current_value = array[r, c]
        
        for dr, dc in directions:
            nr, nc = r + dr, c + dc
            if not (0 <= nr < rows and 0 <= nc < cols):
                continue
            # mark neighbours as +1 further away from center
            #visit_depth[nr, nc] = min(visit_depth[nr, nc], visit_depth[r, c] + 1)
            #if visit_depth[nr, nc] > visit_depth[r, c] + 1:
            #    # only go from center outwards, not back inside
            #    continue
            if visited[nr, nc]:
                continue

            # print('filling pixel:', len(queue), nr, nc)
            # Enforce monotonic decrease on neighbours
            array[nr, nc] = min(array[nr, nc], current_value)
            visited[nr, nc] = True

            # propagate monotonicity from this position
            queue.append((nr, nc))

    return array

AP_UNITS = 'mJy'

def plot_aperture_flux_model(ax, radius, apfluxes, plot_areas, **kwargs):
    return ax.plot(
        radius,
        np.diff(np.array([0] + list(apfluxes))),
        **kwargs)

def plot_aperture_flux_data(ax, x, y, yerr, plot_areas, **kwargs):
    return ax.errorbar(
        x=x,
        y=np.diff(np.array([0] + list(y))),
        yerr=np.diff(np.array([0] + list(yerr))),
        **kwargs)

# read sources we should handle
table = Table.read(sys.argv[1]).filled()
# and their neighbours
neighbour_table = Table.read(sys.argv[2]).filled()

assert 'y_apertureflux_118_flux' in table.colnames
layer = 'hsc-dr3'
bands = 'grizy'
annuli = np.array([10, 15, 20, 30, 40, 57, 84, 118])
radius = annuli / 10.0 / 2
areas = np.pi * radius**2
plot_areas = areas
for band in bands:
    for i in annuli:
        table['apfluxext_%s_%d' % (band, i)] = np.nan
        table['apfluxext_err_%s_%d' % (band, i)] = np.nan
    table['neighbours_' + band + '_saturated'] = False
    table['neighbours_' + band + '_flux_30'] = np.nan
    table['neighbours_' + band + '_flux_50'] = np.nan

# handle each object one by one
for k, obj in enumerate(table):
    # skip objects saturated in z, because the centering is wrong, and then the aperture fluxes are wrong
    if not (obj['z_pixelflags_saturated'] == 0 and obj['z_pixelflags_saturatedcenter'] == 0):
        print("skipping: saturated in z, so centering will not be reliable", obj['id'])
        continue

    # get the neighbours within 5 arcsec:
    pos = SkyCoord(obj['ra'], obj['dec'], unit='deg')
    neighbour_sep = SkyCoord(neighbour_table['ra'], neighbour_table['dec'], unit='deg').separation(pos)
    verynearby_neighbours = neighbour_table[neighbour_sep < 3 * u.arcsec]
    nearby_neighbours = neighbour_table[neighbour_sep < 5 * u.arcsec]
    fig, axs = plt.subplots(
        len(bands), 5, figsize=(8, 1 + 2.5 * len(bands)),
        gridspec_kw=dict(hspace=0, wspace=0), sharex='col', sharey=False
    )
    for i, band in enumerate(bands):
        # plot each band
        print(f'{k * 100 / len(table):.1f}%% [{k}/{len(table)}] {obj["id"]} {band} -------')
        # if neighbours are saturated, do not use this object
        plot_image(axs[i,1], download_image(obj['ra'], obj['dec'], band, layer))
        fits_apertures = None
        psf_apertures = None
        try:
            if obj[band + '_inputcount_value'] > 0:
                origfitsimage, center = get_image_center(download_fits_image(obj['ra'], obj['dec'], band, layer))
                psfimage, psfcenter = get_image_center(download_psf_image(obj['ra'], obj['dec'], band, layer))
                fitsimage = monotonic_floodfill(origfitsimage.copy(), center)
                #fitsimage = origfitsimage.copy()
                plot_image(axs[i,2], image_viz_transform(origfitsimage))
                plot_image(axs[i,3], image_viz_transform(fitsimage))
                plot_image(axs[i,4], image_viz_transform(psfimage))
                axs[i,2].text(
                    0.02, 0.98, 'FITS image', size=6,
                    color='white', transform=axs[i,2].transAxes, va='top')
                axs[i,3].text(
                    0.02, 0.98, 'floodfill', size=6,
                    color='white', transform=axs[i,2].transAxes, va='top')
                axs[i,4].text(
                    0.02, 0.98, 'PSF', size=6,
                    color='white', transform=axs[i,3].transAxes, va='top')
                arcsec_per_px = 0.17
                # get apertures
                fits_apertures = get_aperture_curve_from_image(
                    fitsimage, center, radius / arcsec_per_px, axs[i,3], edgecolor='white', fill=False, linewidth=0.1) / arcsec_per_px**2 * 1.7
                psf_apertures = get_aperture_curve_from_image(
                    psfimage, psfcenter, radius / arcsec_per_px, axs[i,4], edgecolor='white', fill=False, linewidth=0.1) / arcsec_per_px**2 * 1.7
        except OSError as e:
            print('ERROR:', e)
            pass

        # take note of potential issues
        diagnostic_str = f"{obj['ra']:6f} {obj['dec']:6f}\n"
        diagnostic_str += f"#:{obj[band + '_inputcount_value']*1}[{obj[band + '_inputcount_flag']*1}] bkg:{obj[band + '_localbackground_flag']*1}\n"
        diagnostic_str += f"pix:{obj[band + '_pixelflags']*1}: bad:{obj[band + '_pixelflags_bad']*1}  edge:{obj[band + '_pixelflags_edge']*1}\n"
        diagnostic_str += f"sat:{obj[band + '_pixelflags_saturated']*1}: ctr:{obj[band + '_pixelflags_saturatedcenter']*1}\n"
        diagnostic_str += f"flag:{obj[band + '_apertureflux_10_flag']*1}:{obj[band + '_apertureflux_40_flag']*1} trunc:{obj[band + '_apertureflux_10_flag_aperturetruncated']*1}:{obj[band + '_apertureflux_40_flag_aperturetruncated']*1}\n"
        diagnostic_str += f"PSF:{obj[band + '_kronflux_psf_radius']:.2f} {obj[band + '_sdssshape_psf_shape11']**0.5:.2f}/{obj[band + '_sdssshape_psf_shape12']/obj[band + '_sdssshape_psf_shape11']**0.5/obj[band + '_sdssshape_psf_shape22']**0.5:.2f}/{obj[band + '_sdssshape_psf_shape22']**0.5:.2f} flag:{obj[band + '_kronflux_flag_bad_shape_no_psf']*1}\n"

        # store flags, which will be considered in the Makefile.
        # first, we should get the neighbours to know how contaminated the image is, so we can know how far from the center we can safely integrate
        #usable = True
        obj['neighbours_' + band + '_saturated'] = not (np.all(nearby_neighbours[band + '_pixelflags_saturated'] == 0) and np.all(nearby_neighbours[band + '_pixelflags_saturatedcenter'] == 0))
        #if not np.all(nearby_neighbours[band + '_pixelflags_saturated'] == 0 and nearby_neighbours[band + '_pixelflags_saturatedcenter'] == 0):
        #    usable = False
        obj['neighbours_' + band + '_flux_30'] = verynearby_neighbours[band + '_psfflux_flux'].sum()
        #if not verynearby_neighbours[band + '_cmodel_flux'].sum() > obj[band + '_cmodel_flux']:
        #    print("skipping: neighbouring objects within 3'' are brighter", obj['id'])
        #    diagnostic_str += "closeby objects outshine"
        #    usable = False
        obj['neighbours_' + band + '_flux_50'] = nearby_neighbours[band + '_psfflux_flux'].sum()
        #if not nearby_neighbours[band + '_cmodel_flux'].sum() > 10 * obj[band + '_cmodel_flux']:
        #    print("skipping: neighbouring objects within 5'' extremely bright", obj['id'])
        #    diagnostic_str += "neighbours extremely bright"
        #    usable = False
        #aperture_radius = 40
        #if not nearby_neighbours[band + '_cmodel_flux'].sum() > obj[band + '_cmodel_flux']:
        #    print("neighbouring objects are bright, only measuring until 2''", obj['id'])
        #    diagnostic_str += "neighbouring objects are bright, using 2arcsec"
        #    aperture_radius = 20
        
        axs[i,1].text(
            0.02, 0.98, diagnostic_str, size=6,
            color='white', transform=axs[i,1].transAxes, va='top')
        ax = axs[i,0]
        ax.set_ylabel(f'{band} flux [{AP_UNITS}]')
        ax.set_xlabel('Radius [arcsec]')
        # get aperture fluxes from catalog
        apflux = np.array([obj['%s_apertureflux_%d_flux' % (band, i)] for i in annuli])
        apflux_err = np.array([obj['%s_apertureflux_%d_fluxerr' % (band, i)] for i in annuli])

        if not np.isfinite(apflux_err).any():
            continue

        # compare aperture fluxes to image FITS apertures and PSF aperture fluxes (below)
        plot_aperture_flux_data(ax, x=radius, y=apflux, yerr=apflux_err, plot_areas=plot_areas, label='Total', ls='--', color='gray')

        ax.set_xlim(0, None)
        ymin, ymax = ax.get_ylim()
        ax.set_ylim(max(ymin, -ymax*0.05), ymax)

        if psf_apertures is None:
            continue

        plot_aperture_flux_model(ax, radius, psf_apertures / psf_apertures[0] * fits_apertures[0], plot_areas=plot_areas, label='imgPSF', color='red')
        plot_aperture_flux_model(ax, radius, fits_apertures, plot_areas=plot_areas, label='img', color='k')

        # compute surely extended flux:
        # flux in the inner-most aperture, for normalisation
        nctr = unumpy.uarray(apflux[0], apflux_err[0])
        # normalise psf aperture fluxes to first bin
        psf_flux_renormalised = psf_apertures * nctr / psf_apertures[0]
        # get galaxy fluxes
        fluxgal = unumpy.uarray(apflux, apflux_err)

        for j in range(1, len(annuli)):
            # residual flux, after subtracting normalised gaussian
            c = fluxgal[j] - psf_flux_renormalised[j]
            # store residual flux columns
            finite = np.isfinite(unumpy.std_devs(c))
            obj['apfluxext_%s_%d' % (band, annuli[i])] = np.where(finite, unumpy.nominal_values(c), np.nan)
            obj['apfluxext_err_%s_%d' % (band, annuli[i])] = np.where(finite, unumpy.std_devs(c), np.nan)

    # print(f'writing "extflux_{obj["id"]}.pdf"')
    axs[0,0].legend(fontsize=8)
    plt.savefig(f'extflux_{obj["id"]}.pdf')
    plt.close()

table.write(sys.argv[3], overwrite=True)
