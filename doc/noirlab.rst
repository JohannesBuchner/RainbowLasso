.. _noirlab:

================
NOAO Data Lab
================

`Astro Data Lab <https://datalab.noirlab.edu/>`_ hosts the Legacy Survey,
VISTA Hemisphere Survey, UKIDSS and Euclid catalogs that RainbowLasso
queries.

Setup
-----

Data Lab requires a free account. Log in once::

    datalab login

The credentials are stored in your home directory and shared between
environments.

How the pipeline queries Data Lab
---------------------------------

For each survey, the Makefile follows the same pattern (shown here for
VHS; the segments are regenerated with ``python3 makefilesegments.py``,
see :ref:`makefilesegments`):

.. literalinclude:: _segments/vhs-fetch.mk
   :language: make

1. The input catalog is reduced to id, RA, DEC columns
   (``%_coords.csv``, with :doc:`stilts <vizier>`).
2. The coordinates are uploaded into a temporary ``mydb`` table
   (``mydb_import``).
3. A SQL query cross-matches the uploaded coordinates against the survey
   table with ``q3c_join`` within 1 arcsec (the survey source closest to
   each input position is selected). For neighbour queries
   (``%_LS_neighbours_noirlab.fits``), a 20 arcsec radius is used.
4. The result is fetched as a FITS file (``--fmt=fits``).
5. The temporary ``mydb`` tables are dropped again.

Tables used
-----------

.. list-table::
   :header-rows: 1

   * - survey
     - table
     - bands
   * - Legacy Survey
     - ``ls_dr10.tractor`` + ``ls_dr10.apflux``
     - g, r, i, z, W1, W2, W3, W4 (see :doc:`legacysurvey`)
   * - VHS
     - ``vhs_dr5.vhs_cat_v3``
     - Y, J, H, Ks (see :doc:`vhs`)
   * - UKIDSS
     - ``ukidss_dr11plus.lassource``
     - Y, J, H, K (see :doc:`ukidss`)
   * - Euclid Q1
     - ``euclid_q1.object``
     - VIS, NISP (see :doc:`euclid`)

For the Legacy Survey, the tractor photometry is joined with the
aperture flux table ``ls_dr10.apflux`` (via ``ls_id``) to obtain the
aperture fluxes and inverse variances (``apflux_*_1`` .. ``apflux_*_8``,
``apflux_ivar_*_1`` .. ``apflux_ivar_*_8``) that the adaptive flux
calculation needs.
