.. _allwise:

=======
AllWISE
=======

*Catalog*: `AllWISE <https://wise2.ipac.caltech.edu/docs/release/allwise/>`_,
queried through VizieR as table ``II/328/allwise``
(see :doc:`vizier`).

*Wavelengths*: W1 (3.4 um), W2 (4.6 um), W3 (12 um), W4 (22 um).

AllWISE provides detections based on WISE data alone, and is helpful for
providing upper limits by finding the brightest neighbour.

Matching and aggregation
------------------------

.. literalinclude:: _segments/allwise-fetch.mk
   :language: make

For the upper-limit handling, all AllWISE sources within 20 arcsec are
matched (``find=all``); ``groupsumWISE.py`` then sums the fluxes of all
matches per input source (flux sum, error as norm):

.. literalinclude:: _segments/allwise-sum.mk
   :language: make

* The Vega magnitudes are converted to fluxes (AB offsets 2.699, 3.339,
  5.174, 6.620 for W1-W4, and 48.60 to convert to flux units).
* If none of the matches has a detection in a band, the 5 sigma all-sky
  sensitivity is used as an upper limit (68, 111, 860, 5700 mJy for
  W1-W4, `from the AllWISE exposure support
  <https://wise2.ipac.caltech.edu/docs/release/allwise/expsup/sec2_3a.html>`_).

Cuts applied
------------

In the final merge (:ref:`the final rule <final>`), the WISE fluxes of
the Legacy Survey deblending (``LU_flux_w*_LS``) are preferred, and the
AllWISE measurements are only used as fallback:

* **Isolated, good extraction** (``isolatedLS && goodfitsW1/W2``): the
  LS-deblended flux is used (origin ``LS10``).
* **Blended** (``fracflux_w* > 1``) or AllWISE brighter than the LS
  flux plus its error budget: the summed AllWISE flux of all neighbours
  is used (origin ``AllWISE``; negative error convention).
* **Otherwise**: the LS flux inflated by its error and the fracflux
  fraction is used as a conservative 3 sigma upper limit (origin
  ``LS10UL``).
* **W3/W4** are only used when W1/W2 are unblended and isolated
  (``!W34_blended``), because the LS W3/W4 deblending is unreliable.

The origin of each WISE1/WISE2 value is recorded in the ``WISE1_origin``
and ``WISE2_origin`` columns (``LS10``, ``AllWISE`` or ``LS10UL``); see
:doc:`outputs` for the column convention.
