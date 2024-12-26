import numpy as np
import sys
from astropy.table import Table
from uncertainties import unumpy
import scipy.stats

annuli = np.array([10, 15, 20, 30, 40, 57, 84, 118])  # diameter in 0.1 arcsec units
radius = annuli / 10.0 / 2
areas = np.pi * radius**2

table = Table.read(sys.argv[1])

for band in 'grizy':
    # Galaxy aperture fluxes
    apflux = [table['%s_apertureflux_%d_flux' % (band, i)] for i in annuli]
    apflux_err = [table['%s_apertureflux_%d_fluxerr' % (band, i)] for i in annuli]
    # Gaussian PSF, convert PSF FWHM to sigma
    sigma = table[band + '_kronflux_psf_radius']
    # normalisation of the galaxy to match, in the innermost aperture
    ng = unumpy.uarray(apflux[0], apflux_err[0])

    # flux in the inner aperture for a Gaussian PSF
    gaussflux0 = scipy.stats.rayleigh.cdf(radius[0], 0, sigma)

    for i in range(1, len(annuli)):
        fluxg = unumpy.uarray(np.array(apflux[i-1]), np.array(apflux_err[i-1]))
        # flux in this aperture for a Gaussian PSF
        gaussflux = scipy.stats.rayleigh.cdf(radius[i-1], 0, sigma)
        # normalisation of gaussflux to first bin
        normg = ng / (gaussflux0)
        # residual flux, after subtracting normalised gaussian
        c = fluxg - gaussflux * normg #/ areas[i-1]
        # store residual flux columns
        finite = np.isfinite(unumpy.std_devs(c))
        table['apfluxext_%s_%d' % (band, annuli[i])] = np.where(finite, unumpy.nominal_values(c), np.nan)
        table['apfluxext_err_%s_%d' % (band, annuli[i])] = np.where(finite, unumpy.std_devs(c), np.nan)
        
table.write(sys.argv[2], overwrite=True)
