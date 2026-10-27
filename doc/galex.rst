.. _galex:

=====
GALEX
=====

*Catalog*: GUVcat_AIS, the revised catalog of GALEX ultraviolet sources
(`Bianchi, Shiao & Thilker 2017 <https://arxiv.org/abs/1704.05903>`_),
queried through VizieR as table ``II/335/galex_ais``
(see :doc:`vizier`).

*Wavelengths*: far-UV (FUV, ~1540 A) and near-UV (NUV, ~2320 A).

*Depth*: "The typical depth of the GUVcat_AIS catalog is
FUV=19.9, NUV=20.8 ABmag" (from the paper abstract). Measurements
fainter than these depths are unreliable and are discarded.

Matching
--------

.. literalinclude:: _segments/galex.mk
   :language: make

Each input position is matched against GUVcat_AIS with the CDS X-Match
service (``stilts cdsskymatch``) within a 4 arcsec radius, and the
closest match is used (``find=each``).

Cuts applied
------------

* **Depth cuts**: only detections brighter than the catalog depth are
  used, FUVmag < 19.9 and NUVmag < 20.8 (AB). Fainter measurements are
  set to -99 (missing).
* **Flux cuts**: the total fluxes must be brighter than the depth
  equivalent, Fflux > 39.7 uJy (FUV) and Nflux > 17.4 uJy (NUV).
* **Extraction quality**: the artifact flag (``Fafl``, ``Nafl``) and the
  extraction flag (``Fexf``, ``Nexf``) must be zero
  (``Fgood_extraction``, ``Ngood_extraction``).
* **Milky Way extinction**: the magnitudes are dereddened with
  ``EBmV`` from the catalog (8.286*EBmV in FUV, 8.621*EBmV in NUV), and
  the energy fluxes are corrected by ``10^(0.4*A_lam)``.

Aperture fluxes
---------------

The GALEX PSF is very large (FWHM ~4-5 arcsec), so all objects at
z >> 0 are treated as point sources. The pipeline stores the aperture
magnitudes within a 4 pixel radius (pixel size 1.5 arcsec, i.e. a 6
arcsec radius aperture, ``FUV_4``/``NUV_4``), converted to AB
magnitudes by adding the AB zeropoints (+18.82 in FUV, +20.08 in NUV),
alongside the total magnitudes.

No upper limits
---------------

GALEX upper limits are deliberately *not* provided, because:

* only no-entry cases (no match within the search radius) would be
  suitable;
* in the Milky Way, the upper limit would need to be corrected up, and
  thus be drastically weaker (higher);
* GALEX is very shallow already, so those upper limits are probably not
  so informative.

Instead, for sources without a GALEX detection, a 3 sigma NUV upper
limit is estimated from the survey coverage, see :doc:`galex-coverage`.

Output columns
--------------

The Makefile rule produces, per band (prefix F = FUV, N = NUV):

* ``FUV``, ``NUV``: dereddened total magnitudes (AB).
* ``FUV_aper``, ``NUV_aper``: dereddened aperture magnitudes (AB).
* ``Fflux``, ``Nflux``, ``Fflux_aper``, ``Nflux_aper``: energy fluxes
  and errors (uJy, AB-based).
* ``Fflux_real_LU``, ``Nflux_real_LU`` (and ``_aper`` variants):
  extinction-corrected energy fluxes in the internal LU unit
  (erg/s/cm^2/Hz), which feed the final ``FUV``/``NUV`` mJy columns of
  the merged catalog (see :doc:`outputs`).

References
----------

* `Bianchi, Shiao & Thilker 2017 <https://arxiv.org/abs/1704.05903>`_
  (GUVcat_AIS, ApJS 230, 24): the catalog, its depth and columns
  (documented in the `VizieR II/335 ReadMe
  <http://vizier.cds.unistra.fr/viz-bin/ReadMe/II/335>`_).
* `Morrissey et al. 2007 <https://arxiv.org/abs/astro-ph/0703531>`_
  (ApJS 173, 682): the GALEX instrument, calibration and AB zeropoints
  (see also the `GALEX FAQ on counts and background
  <https://asd.gsfc.nasa.gov/archive/galex/FAQ/counts_background.html>`_).
* `Schlafly & Finkbeiner 2011 <https://arxiv.org/abs/1009.2514>`_
  (ApJ 737, 103): the per-band Milky Way extinction coefficients used
  for the dereddening (8.286 in FUV, 8.621 in NUV).
