.. _ukidss:

================================
UKIDSS (near-infrared)
================================

*Catalog*: `UKIDSS DR11plus Large Area Survey
<http://wsa.roe.ac.uk/dr11plus.html>`_, queried through the NOAO Data
Lab (see :doc:`noirlab`): ``ukidss_dr11plus.lassource``.

*Wavelengths*: Y, J, H, K.

Matching
--------

.. literalinclude:: segments/ukidss-fetch.mk
   :language: make

The nearest source within 1 arcsec is retrieved (SQL ``q3c_join`` via
Data Lab).

Aperture fluxes
---------------

.. literalinclude:: segments/ukidss-fluxes.mk
   :language: make

The flux extraction is adapted to the source type:

* for **point sources**: aperture fluxes with a 2.8 arcsec diameter
  aperture (``apermag4``, 1.4 arcsec radius),
* for **extended sources**: aperture fluxes with a 5.7 arcsec diameter
  (``apermag6``, 2.8 arcsec radius).

These apertures follow the CASU aperture sequence (radii 0.5, 0.71,
1.0, 1.4, 2.0, 2.8, 4.0, 5.66 arcsec, steps of sqrt(2); AperMag3 is
the 1.0'' radius default aperture). See :doc:`apertures` for the
matched-aperture table.

.. note::

   There are ONLY aperture-corrected fluxes provided in the UKIDSS
   catalog, also for the large aperture.

The per-source Milky Way extinction (``ay``, ``aj``, ``ah``, ``ak``
columns) is subtracted from the magnitudes, the constants (Y: 0.634,
J: 0.938, H: 1.379, K: 1.9) are the conversion from Vega to AB
magnitudes (from
`Hewett et al. <http://wsa.roe.ac.uk/pubs/hewett-ukidss.pdf>`_), and
48.6 converts to flux units.

Cuts applied
------------

* **Error bits** (``Yerrbits`` .. ``Kerrbits``): the band is only used
  when the error bits are zero.
* The choice between point-like and extended aperture is made with the
  ``pointlike`` flag (from the Legacy Survey morphology), see
  :ref:`the final rule <final>`.

.. note::

   Blended or problematic photometry in UKIDSS is not discarded yet
   (see the TODO list in the README), and upper limits are not provided
   (they would need coverage information).

Output columns
--------------

``WFCAM_Yapc4flux``, ``WFCAM_Japc4flux`` .. (and ``_err``): aperture
fluxes in the internal LU unit, which feed the final ``WFCAM_Y``,
``WFCAM_J``, ``WFCAM_H``, ``WFCAM_Ks`` mJy columns of the merged
catalog (see :doc:`outputs`).

References
----------

* `Lawrence et al. 2007 <https://arxiv.org/abs/astro-ph/0604426>`_
  (MNRAS 379, 1599): the UKIDSS survey.
* `Hewett et al. 2006 <https://arxiv.org/abs/astro-ph/0604427>`_
  (MNRAS 367, 454): the UKIDSS photometric system and the Vega to AB
  conversion (`WSA copy
  <http://wsa.roe.ac.uk/pubs/hewett-ukidss.pdf>`_).
* `WSA lassource column documentation
  <http://wsa.roe.ac.uk/las_lassource.html>`_: the catalog columns
  (aperture magnitudes, error bits, per-source extinction ``ay`` ..
  ``ak``).
