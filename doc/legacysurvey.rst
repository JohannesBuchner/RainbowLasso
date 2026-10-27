.. _legacysurvey:

=========================
Legacy Survey (optical)
=========================

*Catalog*: `Legacy Survey DR10 <https://www.legacysurvey.org/dr10/>`_,
queried through the NOAO Data Lab (see :doc:`noirlab`):
``ls_dr10.tractor`` joined with the aperture flux table
``ls_dr10.apflux``.

*Wavelengths*: g, r, i, z from DECam (decam), r, g from 90prime
(BASS), z from MzLS (zd_mosaic), plus W1-W4 from WISE (deblended by the
Legacy Survey).

Only sources within the Legacy Survey area are handled.

Adaptive fluxes
---------------

.. literalinclude:: segments/legacysurvey-fluxes.mk
   :language: make

The flux extraction is adapted to the source type:

* for **point sources**: the survey model/PSF flux (``flux_*``),
* for **extended sources**: the flux within a 3.5 arcsec radius aperture
  (7 arcsec diameter; ``apflux_*_6`` for the optical bands,
  ``apflux_*_1`` for the WISE bands, whose aperture radii are
  band-dependent),

corrected for Milky Way extinction with the ``MW_TRANSMISSION_*``
columns. The flux unit is nanomaggy; multiplying by ``10^-28.44``
(= 3.63e-29) converts to erg/s/cm^2/Hz (the internal LU unit).

The resulting columns are ``adflux_*`` and ``LU_flux_*`` (with errors).

Cuts applied
------------

The quality cuts are applied in the final merge
(:ref:`the final rule <final>`), where the per-band columns
``decam_*``, ``90prime_*``, ``zd_mosaic`` are computed:

.. literalinclude:: segments/final.mk
   :language: make

* **Image fitting problems** (``fitbits_LS``): for point sources, the
  bits 1 | 4 | 8 of the Legacy Survey fitbits must be zero (otherwise
  the band is discarded).
* **Mask bits** (``maskbits_LS``): bits 1024 | 2048 (bright star
  haloes/diffraction spikes regions) exclude the source for all bands;
  bit 4 (g), 8 (r), 16 (z), 16384 (i) exclude the respective band.
* **Fraction in aperture** (``fracin_*``): only sources with
  ``fracin_* > 0.5`` are used; the error is conservatively inflated by
  ``1/fracin_*`` because up to that fraction of the flux may be missed.
* **Fractional flux** (``fracflux_*``): the fraction of the flux of
  neighbouring sources; sources with fracflux > 0.1 in g, r, z, W1, W2
  (or > 10 in W3, W4) are flagged as not isolated (``isolatedLS``).
* **Sky region** (``inMzLSBASS``): sources in the BASS/MzLS area
  (DEC>32, 90<RA<300) get their r, g, z fluxes in the ``90prime_r``,
  ``90prime_g``, ``zd_mosaic`` columns instead of ``decam_r``,
  ``decam_g``, ``decam_z``.

Extended-source flux priors
---------------------------

For the ``%_all_extflux.fits`` products, additional galaxy total-flux
priors (``prior_GALflux_decam_*``) are computed from the LS aperture
fluxes (aperture radius 4, i.e. 11.4 arcsec diameter), for extended
sources only:

.. literalinclude:: segments/legacysurvey-extflux.mk
   :language: make

.. literalinclude:: segments/legacysurvey-extflux-merge.mk
   :language: make

These require ``TYPE != PSF``, ``fracflux < 0.1`` and positive errors,
and are only used if they exceed 1% of the decam flux.

.. _final:

Final rule
----------

The ``%_all.fits`` rule merges all catalogs together (GALEX, LS, UKIDSS,
VHS, AllWISE, GALEX coverage) and computes the final flux columns; see
:doc:`outputs` for the column convention.

References
----------

* `Dey et al. 2019 <https://arxiv.org/abs/1804.08657>`_ (AJ 157, 168):
  the Legacy Surveys overview paper, describing the Tractor
  photometry, the source types and the mask bits.
* `Legacy Survey DR10 documentation <https://www.legacysurvey.org/dr10/catalogs/>`_:
  the catalog column descriptions (fluxes in nanomaggy, aperture radii,
  MW_TRANSMISSION), and the `bitmask descriptions
  <https://www.legacysurvey.org/dr10/bitmasks/>`_ for ``maskbits`` and
  ``fitbits``.
* `Lang, Hogg & Schlegel 2014 <https://arxiv.org/abs/1405.7368>`_
  (AJ 148, 26): the unWISE coadds underlying the WISE deblended
  photometry in the Legacy Survey.
* `Schlafly & Finkbeiner 2011 <https://arxiv.org/abs/1009.2514>`_
  (ApJ 737, 103): the MW transmission coefficients used for the
  extinction correction.
