.. _apertures:

=================
Matched apertures
=================

To model the same object at all wavelengths, the same aperture needs to
be extracted ("matched aperture"). Depending on whether the source is
extended or point-like, RainbowLasso uses either PSF-like photometry
(weighted by the survey point-spread-function, best for stars and
distant quasars) or a total flux based on an extraction aperture (for
galaxies). The choice is made with the ``pointlike`` flag, taken from
the Legacy Survey morphology.

This table lists the imaging resolutions and the aperture radii used:

.. list-table::
   :header-rows: 1
   :widths: 14 16 24 28 18

   * - survey / instrument
     - imaging resolution (PSF FWHM)
     - point-like measurement
     - extended-source measurement
     - final columns
   * - GALEX AIS
     - ~4.2'' (FUV), ~5.3'' (NUV)
     - PSF-weighted total flux (all objects are treated as point
       sources); additionally the 6'' radius aperture (``FUV_4``,
       4 pixel x 1.5'') is stored
     - (same as point-like)
     - ``FUV``, ``NUV``
   * - Legacy Survey (DECam)
     - seeing ~1.0-1.5''
     - model flux (PSF-like fit, ``flux_*``)
     - 3.5'' radius aperture (``apflux_*_6``, 7'' diameter)
     - ``decam_g/r/i/z``
   * - Legacy Survey (90prime, MzLS)
     - seeing ~1.5''
     - model flux
     - 3.5'' radius aperture
     - ``90prime_g/r``, ``zd_mosaic``
   * - Legacy Survey WISE (unWISE deblended)
     - ~6.1'' (W1) .. 12'' (W4)
     - profile-fit flux (model, ``flux_w*``)
     - 3.0'' radius aperture (``apflux_w1_1``, 6'' diameter)
     - ``WISE1..4``
   * - VHS (VISTA)
     - seeing ~0.7-1.0''
     - 2.8'' diameter aperture, aperture-corrected (``apermag4``, 1.4''
       radius)
     - 5.7'' diameter aperture, no aperture correction
       (``apermagnoapercorr6``)
     - ``UV_Y/J/H/K``
   * - UKIDSS (WFCAM)
     - seeing ~0.7-1.0''
     - 2.8'' diameter aperture, aperture-corrected (``apermag4``)
     - 5.7'' diameter aperture, aperture-corrected (``apermag6``)
     - ``WFCAM_Y/J/H/Ks``
   * - AllWISE
     - ~6.1'' (W1), ~6.4'' (W2), ~7.4'' (W3), ~12'' (W4)
     - flux summed over all matches within 20'' (origin ``AllWISE``),
       used as fallback and for the upper limits
     - (same as point-like)
     - ``WISE1..4`` (fallback)
   * - HSC PDR3 (optional)
     - seeing ~0.6-0.7''
     - PSF model flux (``*_psfflux_flux``)
     - 5.75'' radius aperture (``*_apertureflux_57_flux``, 11.5''
       diameter)
     - ``HSC_g..y``
   * - SDSS (optional)
     - seeing ~1.4''
     - PSF magnitudes (``psfMag_*``)
     - 3.0'' radius aperture (``profMean`` bin 4, 6'' diameter)
     - ``u..z_sdss``

Notes
-----

* The VISTA/WFCAM aperture magnitudes follow the CASU aperture
  sequence (radii 0.5, 0.71, 1.0, 1.4, 2.0, 2.8, 4.0, 5.66 arcsec,
  steps of sqrt(2)): ``apermag4`` is a 1.4'' radius aperture (2.8''
  diameter) and ``apermag6`` a 2.8'' radius aperture (5.7'' diameter).
  The VHS aperture sizes are documented in the VizieR `II/367 column
  descriptions <https://vizier.cds.unistra.fr/viz-bin/VizieR?-source=II/367>`_,
  and AperMag3 (1.0'' radius) is the UKIDSS default aperture
  (`II/319/las9 <https://vizier.cds.unistra.fr/viz-bin/VizieR?-source=II/319/las9>`_).
* The GALEX aperture magnitudes are stored with a 4 pixel radius
  (1.5'' pixel, i.e. a 6'' radius aperture); with its very large PSF,
  all objects are treated as point sources (see :doc:`galex`).
* Euclid Q1 (optional) provides no aperture fluxes with suitable
  diameters; it only contributes an extended-fraction prior to the
  optical flux columns (see :doc:`euclid`).
* The extended-source apertures are chosen so that the same physical
  aperture is measured at all wavelengths as far as possible, given the
  resolutions above.
