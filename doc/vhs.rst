.. _vhs:

================================
VISTA Hemisphere Survey (VHS)
================================

*Catalog*: `VISTA Hemisphere Survey DR5
<http://casu.roe.ac.uk/vhsDR5/>`_ (VHS), queried through the NOAO Data
Lab (see :doc:`noirlab`): ``vhs_dr5.vhs_cat_v3``.

*Wavelengths*: Y, J, H, Ks.

Matching
--------

.. literalinclude:: segments/vhs-fetch.mk
   :language: make

The nearest source within 1 arcsec is retrieved (SQL ``q3c_join`` via
Data Lab).

Aperture fluxes
---------------

.. literalinclude:: segments/vhs-fluxes.mk
   :language: make

The flux extraction is adapted to the source type:

* for **point sources**: aperture-corrected fluxes within a 2.8 arcsec
  diameter aperture (``apermag4``, 1.4 arcsec radius),
* for **extended sources**: not-aperture-corrected fluxes within a
  5.7 arcsec diameter aperture (``apermagnoapercorr6``, 2.8 arcsec
  radius).

These apertures follow the CASU aperture sequence (radii 0.5, 0.71,
1.0, 1.4, 2.0, 2.8, 4.0, 5.66 arcsec, steps of sqrt(2)); the sizes are
documented in the VizieR `II/367 column descriptions
<https://vizier.cds.unistra.fr/viz-bin/VizieR?-source=II/367>`_
("2.8\" aperture diam", "5.7\" aperture diam"). See
:doc:`apertures` for the matched-aperture table.

The per-source Milky Way extinction (``ay``, ``aj``, ``ah``, ``aks``
columns) is subtracted from the magnitudes, the constants (0.634, 0.938,
1.379, 1.9) are the conversion from Vega to AB, and 48.6 converts to
flux units.

Cuts applied
------------

* **Error bits** (``Yerrbits`` .. ``Kserrbits``): the band is only used
  when the error bits are zero.
* The choice between point-like and extended aperture is made with the
  ``pointlike`` flag (from the Legacy Survey morphology), see
  :ref:`the final rule <final>`.

.. note::

   Blended or problematic photometry in VHS is not discarded yet (see
   the TODO list in the README), and upper limits are not provided
   (they would need coverage information).

Output columns
--------------

``UV_Yapc4flux``, ``UV_Japc4flux`` .. (and ``_err``): aperture fluxes
in the internal LU unit, which feed the final ``UV_Y``, ``UV_J``,
``UV_H``, ``UV_K`` mJy columns of the merged catalog (see
:doc:`outputs`).

References
----------

* `McMahon et al. 2013 <https://ui.adsabs.harvard.edu/abs/2013Msngr.154...35M/abstract>`_
  (The Messenger 154, 35): the VHS survey description.
* `VSA/WSA documentation <http://wsa.roe.ac.uk/>`_ and the `VHS DR5
  pages <http://wsa.roe.ac.uk/vhsDR5.html>`_: the catalog columns
  (aperture magnitudes, error bits, per-source extinction ``ay``
  .. ``aks``).
