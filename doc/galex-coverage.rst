.. _galex-coverage:

======================
GALEX upper limits
======================

For sources without a GALEX detection, the pipeline estimates a 3 sigma
NUV upper limit from the GALEX AIS coverage and the Milky Way dust
column (no FUV upper limit is provided, see :doc:`galex`).

Coverage table
--------------

``galex_ais_ctrs.fits`` contains the center positions of the GALEX AIS
tiles (columns ``RAfdeg``, ``DEfdeg``).

``galexebv.py`` adds the SFD Milky Way E(B-V) at each tile center
(via the ``dustmaps`` package) and writes ``galex_ais_ctrs_ebv.fits``:

.. literalinclude:: _segments/galexebv.mk
   :language: make

Upper limit per position
------------------------

Each input position is matched to the nearest GALEX AIS tile center
within 0.56 degrees (1980 arcsec, the tile diagonal). The NUV flux limit
is scaled up with the Galactic extinction towards the tile center:

.. literalinclude:: _segments/galex-coverage.mk
   :language: make

The rule computes ``NUV_fluxlim = 0.02 * 10^(EBV_min * 1.35)`` (mJy),
i.e. the 3 sigma sensitivity (0.02 mJy) degraded by the dust
transmission. In the final merge (see :doc:`legacysurvey` and
:doc:`outputs`), sources without a GALEX detection get
``NUV = 3 * NUV_fluxlim`` (a 3 sigma upper limit, negative error
convention) and ``NUV_err = -NUV``.

References
----------

* `Schlegel, Finkbeiner & Davis 1998 <https://arxiv.org/abs/astro-ph/9803161>`_
  (ApJ 500, 525): the SFD Milky Way dust maps that ``galexebv.py``
  queries (via the `dustmaps <https://github.com/gregreen/dustmaps>`_
  package, `Green 2018 <https://doi.org/10.21105/joss.00695>`_).
* `Bianchi, Shiao & Thilker 2017 <https://arxiv.org/abs/1704.05903>`_
  (GUVcat_AIS): the AIS tile grid that the coverage table
  ``galex_ais_ctrs.fits`` is built from.
