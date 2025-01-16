import numpy as np
import sys
from astropy.table import Table
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
import scipy.stats
from PIL import Image
from requests_ratelimiter import LimiterMixin
from requests_cache import CacheMixin
from requests import Session
import astropy.io.fits as pyfits
from io import BytesIO
import os
from photutils.aperture import aperture_photometry, CircularAperture

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

def download_image(ra, dec, band, layer):
    url = f'https://www.legacysurvey.org/viewer/jpeg-cutout?ra={ra}&dec={dec}&zoom=14&layer={layer}&bands={band}'
    print(url)
    response = my_session.get(url)
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

def download_psf_image(ra, dec, band, layer):
    url = f'https://hsc-release.mtk.nao.ac.jp/psf/pdr3/cgi/getpsf?ra={ra}&dec={dec}&filter={band}&rerun=pdr3_wide&tract=&patch=&centered=true&type=coadd'
    print(url)
    response = my_session.get(url, auth=auth)
    print(url, response.status_code)
    f = pyfits.open(BytesIO(response.content))[0]
    center_header = x0, y0 = (-(f.header['CRVAL1A']), -(f.header['CRVAL2A']))
    image = np.array(f.data[::-1,:])
    nx, ny = image.shape[0], image.shape[1]
    y0, x0 = np.unravel_index(image.argmax(), image.shape)
    center = x0, y0
    print(nx, ny, x0, y0, (nx - 1)/2 - x0, (ny - 1)/2 - y0)
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

def download_fits_image(ra, dec, band, layer):
    url = f'https://hsc-release.mtk.nao.ac.jp/das_cutout/pdr3/cgi-bin/cutout?ra={ra}&dec={dec}&sw=0.001&sh=0.001&type=coadd&image=on&filter=HSC-{band.upper()}&tract=&rerun=pdr3_wide'
    print(url)
    response = my_session.get(url, auth=auth)
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
    #ax.plot([nx / 2 - 4, nx / 2 - 10], [ny / 2, ny / 2])
    #ax.plot([nx / 2 + 10, nx / 2 + 4], [ny / 2, ny / 2])
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

# cumulative area-normalised:
#AP_UNITS = 'mJy / arcsec$^2$'
# differential:
AP_UNITS = 'mJy'

def plot_aperture_flux_model(ax, radius, apfluxes, plot_areas, **kwargs):
    # cumulative style:
    #return ax.plot(
    #    radius, np.array(apfluxes), **kwargs)
    # area-normalised cumulative style:
    #return ax.plot(
    #    radius, np.array(apfluxes) / plot_areas, **kwargs)
    # differential style
    #aper_areas = np.diff(np.array([0] + list(plot_areas)))
    return ax.plot(
        radius,
        np.diff(np.array([0] + list(apfluxes))),
        **kwargs)

def plot_aperture_flux_data(ax, x, y, yerr, plot_areas, **kwargs):
    # cumulative style:
    #return ax.errorbar(
    #    x=x, y=np.array(y), yerr=np.array(yerr),
    #    **kwargs)
    # area-normalised cumulative style:
    #return ax.errorbar(
    #    x=x, y=np.array(y) / plot_areas, yerr=np.array(yerr) / plot_areas,
    #    **kwargs)
    # differential style
    #aper_areas = np.diff(np.array([0] + list(plot_areas)))
    return ax.errorbar(
        x=x,
        y=np.diff(np.array([0] + list(y))),
        yerr=np.diff(np.array([0] + list(yerr))),
        **kwargs)

table = Table.read(sys.argv[1]).filled()

if 'y_apertureflux_118_flux' in table.colnames:
    layer = 'hsc-dr3'
    bands = 'grizy'
    annuli = np.array([10, 15, 20, 30, 40, 57, 84, 118])
    radius = annuli / 10.0 / 2
    areas = np.pi * radius**2
elif 'apflux_g_1' in table.colnames:
    layer = 'ls-dr10-grz'
    bands = 'grz'
    annuli = np.arange(1, 9)
    radius = np.array([0.5, 0.75, 1.0, 1.5, 2.0, 3.5, 5.0, 7.0])
    areas = np.pi * radius**2
else:
    assert False, table.colnames

if len(sys.argv) == 3:
    if table['id'].dtype.kind == 'i':
        objs = table[table['id'] == int(sys.argv[2])]
    else:
        objs = table[table['id'] == sys.argv[2]]
else:
    if layer == 'ls-dr10-grz':
        objs = table[table['maskbits'] == 0]
    else:
        objs = table
        for band in bands:
            objs = objs[objs[band + '_pixelflags_saturated'] == 0]
            objs = objs[objs[band + '_pixelflags_saturatedcenter'] == 0]

plot_areas = areas
for obj in objs[:20]:
    fig, axs = plt.subplots(
        len(bands), 5, figsize=(8, 1 + 2.5 * len(bands)),
        gridspec_kw=dict(hspace=0, wspace=0), sharex='col', sharey=False
    )
    for i, band in enumerate(bands):
        #if obj[band + '_inputcount_value'] == 0:
        #    continue
        plot_image(axs[i,1], download_image(obj['ra'], obj['dec'], band, layer))
        fits_apertures = None
        psf_apertures = None
        try:
            if obj[band + '_inputcount_value'] > 0:
                origfitsimage, center = get_image_center(download_fits_image(obj['ra'], obj['dec'], band, layer))
                psfimage, psfcenter = get_image_center(download_psf_image(obj['ra'], obj['dec'], band, layer))
                fitsimage = monotonic_floodfill(origfitsimage.copy(), center)
                #fitsimage = origfitsimage.copy()
                plot_image(axs[i,2], np.log10(origfitsimage))
                plot_image(axs[i,3], np.log10(fitsimage))
                plot_image(axs[i,4], np.log10(psfimage))
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
                fits_apertures = get_aperture_curve_from_image(
                    fitsimage, center, radius / arcsec_per_px, axs[i,3], edgecolor='white', fill=False, linewidth=0.1) / arcsec_per_px**2 * 1.7
                psf_apertures = get_aperture_curve_from_image(
                    psfimage, psfcenter, radius / arcsec_per_px, axs[i,4], edgecolor='white', fill=False, linewidth=0.1) / arcsec_per_px**2 * 1.7
        except OSError as e:
            print('ERROR:', e)
            pass

        diagnostic_str = f"{obj['ra']:6f} {obj['dec']:6f}\n"
        if layer == 'hsc-dr3':
            diagnostic_str += f"#:{obj[band + '_inputcount_value']*1}[{obj[band + '_inputcount_flag']*1}] bkg:{obj[band + '_localbackground_flag']*1}\n"
            diagnostic_str += f"pix:{obj[band + '_pixelflags']*1}: bad:{obj[band + '_pixelflags_bad']*1}  edge:{obj[band + '_pixelflags_edge']*1}\n"
            diagnostic_str += f"sat:{obj[band + '_pixelflags_saturated']*1}: ctr:{obj[band + '_pixelflags_saturatedcenter']*1}\n"
            diagnostic_str += f"flag:{obj[band + '_apertureflux_10_flag']*1}:{obj[band + '_apertureflux_40_flag']*1} trunc:{obj[band + '_apertureflux_10_flag_aperturetruncated']*1}:{obj[band + '_apertureflux_40_flag_aperturetruncated']*1}\n"
            diagnostic_str += f"PSF:{obj[band + '_kronflux_psf_radius']:.2f} {obj[band + '_sdssshape_psf_shape11']**0.5:.2f}/{obj[band + '_sdssshape_psf_shape12']/obj[band + '_sdssshape_psf_shape11']**0.5/obj[band + '_sdssshape_psf_shape22']**0.5:.2f}/{obj[band + '_sdssshape_psf_shape22']**0.5:.2f} flag:{obj[band + '_kronflux_flag_bad_shape_no_psf']*1}\n"
        else:
            diagnostic_str += f"mask:{obj['fracmasked_' + band]*100:.1f}% in:{obj['fracin_' + band]*100:.1f}% flux:{obj['fracflux_' + band]*100:.1f}%\n"
            diagnostic_str += f"mask:{obj['maskbits']} fit:{obj['fitbits']}"
        axs[i,1].text(
            0.02, 0.98, diagnostic_str, size=6,
            color='white', transform=axs[i,1].transAxes, va='top')
        ax = axs[i,0]
        ax.set_ylabel(f'{band} flux [{AP_UNITS}]')
        ax.set_xlabel('Radius [arcsec]')
        if layer == 'hsc-dr3':
            apflux = np.array([obj['%s_apertureflux_%d_flux' % (band, i)] for i in annuli])
            apflux_err = np.array([obj['%s_apertureflux_%d_fluxerr' % (band, i)] for i in annuli])
        else:
            apflux = np.array([obj['apflux_%s_%d' % (band, i)] for i in annuli])
            apflux_err = np.array([obj['apflux_ivar_%s_%d' % (band, i)]**-0.5 for i in annuli])
        if not np.isfinite(apflux_err).any():
            continue
        plot_aperture_flux_data(ax, x=radius, y=apflux, yerr=apflux_err, plot_areas=plot_areas, label='Total', ls='--', color='gray')

        # Gaussian PSF, convert PSF FWHM to sigma
        if layer == 'hsc-dr3':
            # from Ihara+18
            #sigma = dict(g=0.72, r=0.67, i=0.56, z=0.63, y=0.64)[band] / 2.355
            #print(f'sigma[{band}]: {sigma:.3f} {obj[band + "_convolvedflux_seeing"]:.3f} {obj[band + "_kronflux_psf_radius"]:.3f}')
            #sigma = obj[band + '_convolvedflux_seeing'] # / 2.355
            sigma = obj[band + '_kronflux_psf_radius']
            #sigma = np.sqrt(sigma**2 + obj[band + '_convolvedflux_seeing']**2 + obj[band + '_kronflux_psf_radius']**2)
            #sigma = np.sqrt(sigma**2 + obj[band + '_kronflux_psf_radius']**2)
            # sigma = np.sqrt(sigma**2 + obj[band + '_convolvedflux_seeing']**2) / 2.355
            #sigma = 1.0
            # print(f'sigma[{band}]: {sigma:.3f}')
        else:
            sigma = obj['psfsize_' + band] / 2.355
        # flux in the inner aperture for a Gaussian PSF
        if not sigma > 0:
            continue
        # match gaussian PSF in the center
        gaussflux = scipy.stats.rayleigh.cdf(radius, loc=0, scale=sigma)
        gaussflux2 = gaussflux # * 0.85 + 0.15 * scipy.stats.rayleigh.cdf(radius, loc=0, scale=sigma * 3.0)
        plot_aperture_flux_model(ax, radius, apflux[0] * gaussflux / gaussflux[0], plot_areas=plot_areas, label='PSF', color='pink', ls=':')
        # bkgflux = (apflux[-2] - apflux[-1]) / (areas[-2] - areas[-1]) * areas
        # match flux in outer-most annulus as a flat background
        # ax.plot(radius, (apflux[0] / gaussflux[0] * gaussflux + bkgflux) / plot_areas, label='PSF+bkg')

        if layer == 'ls-dr10-grz':
            extflux = [obj['apfluxext_%s_%d' % (band, i)] for i in annuli]
            extflux_err = [obj['apfluxext_err_%s_%d' % (band, i)] for i in annuli]
            # print(gaussflux, apflux[0], extflux)

            plot_aperture_flux_data(ax, x=radius, y=extflux, plot_areas=plot_areas, yerr=np.array(extflux_err) / plot_areas, label='Extended', color='green')
        ax.set_xlim(0, None)
        ymin, ymax = ax.get_ylim()
        ax.set_ylim(max(ymin, -ymax*0.05), ymax)

        if psf_apertures is not None:
            plot_aperture_flux_model(ax, radius, psf_apertures / psf_apertures[0] * fits_apertures[0], plot_areas=plot_areas, label='imgPSF', color='red')
            np.savetxt(f'psfapertures_{obj["id"]}_{band}.txt', np.transpose([radius, psf_apertures / psf_apertures[0], plot_areas]))
        if fits_apertures is not None:
            plot_aperture_flux_model(ax, radius, fits_apertures, plot_areas=plot_areas, label='img', color='k')

    print(f'writing "extflux_{obj["id"]}.pdf"')
    axs[0,0].legend(fontsize=8)
    plt.savefig(f'extflux_{obj["id"]}.pdf')
    plt.close()
