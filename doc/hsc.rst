.. _hsc:

================================
HSC-WIDE (optional)
================================

*Catalog*: `Hyper Suprime-Cam PDR3 <https://hsc-release.mtk.nao.ac.jp/doc/>`_
(WIDE), queried through the hsc-release catalog_jobs API (see
``fetchHSC.py``).

*Wavelengths*: g, r, i, z, y.

Matching
--------

.. literalinclude:: _segments/hsc-fetch.mk
   :language: make

The nearest source within 1 arcsec is retrieved (and separately the
neighbours within 5 arcsec, ``%_HSC_neighbours_aper.fits``, to judge the
contamination of the images).

Adaptive fluxes
---------------

.. literalinclude:: _segments/hsc-fluxes.mk
   :language: make

The flux extraction is adapted to the source type (``extended_*``):
point sources use the 5.75 arcsec radius aperture flux
(``*_apertureflux_57_flux``), extended sources the PSF flux
(``*_psfflux_flux``).

.. warning::

   This assignment (``extended ? psfflux : apertureflux_57``) disagrees
   with the comment in the Makefile (which says point sources use the
   model flux and extended sources the 5.7 arcsec aperture, as for the
   Legacy Survey). Treat the HSC fluxes with care until this is
   resolved.

The Milky Way extinction correction uses the extinction coefficients
``a_*`` and the conversion ``10^-29`` to flux units.

Cuts applied
------------

In the merged ``%_HSCall.fits`` file, a band is only used when:

* the pixel flags (saturated, saturated-center, interpolation,
  local-background, input-count) are clean (``psf_*_usable``),
* the photometry flags (``*_psfflux_flag``,
  ``*_apertureflux_57_flag``) are zero.

Extended-source flux priors
---------------------------

For the ``%_HSCall_extflux.fits`` products, galaxy total-flux priors
(``prior_GALflux_HSC_*``) are computed from the large 11.4 arcsec
aperture fluxes (``apfluxext_*_40``), for extended, unblended,
unsaturated sources:

.. literalinclude:: _segments/hsc-extflux.mk
   :language: make

Output columns
--------------

``HSC_g``, ``HSC_r``, ``HSC_i``, ``HSC_z``, ``HSC_y`` (and ``_err``):
fluxes in mJy in the ``%_HSCall.fits`` file.
