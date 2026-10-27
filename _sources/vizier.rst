.. _vizier:

.. highlight:: shell

=============
CDS / VizieR
=============

The `CDS X-Match service <https://cdsxmatch.u-strasbg.fr/xmatch>`_ and
`VizieR <https://vizier.cds.unistra.fr/>`_ are queried through the
``stilts cdsskymatch`` task for the all-sky catalogs (GALEX AIS, AllWISE)
which do not require a login.

How the pipeline queries VizieR
-------------------------------

For example, the GALEX match (see :doc:`galex`):

.. literalinclude:: _segments/galex.mk
   :language: make

and the AllWISE match (see :doc:`allwise`):

.. literalinclude:: _segments/allwise-fetch.mk
   :language: make

Tables used
-----------

.. list-table::
   :header-rows: 1

   * - catalog
     - table
     - matched within
     - find
   * - GALEX AIS
     - ``II/335/galex_ais``
     - 4 arcsec
     - each
   * - AllWISE
     - ``II/328/allwise``
     - 4 arcsec (``%_ALLWISE.fits``)
     - each
   * - AllWISE
     - ``II/328/allwise``
     - 20 arcsec (``%_ALLWISE_all.fits``)
     - all

Notes
-----

* ``cdstable`` needs the full table identifier including the table name
  (e.g. ``II/335/galex_ais``, not just ``II/335``).
* The ``radius`` of the X-Match service is capped at 180 arcsec
  server-side.
* ``find=all`` output carries an ``angDist`` column (arcsec) with the
  distance to each match, used by the AllWISE aggregation
  (:doc:`allwise`).

References
----------

* `Ochsenbein, Bauer & Marcout 2000 <https://ui.adsabs.harvard.edu/abs/2000A%26AS..143...23O/abstract>`_
  (A&AS 143, 23): the VizieR catalog service.
* `CDS X-Match service <https://cdsxmatch.u-strasbg.fr/xmatch>`_:
  the position cross-matching API behind ``cdsskymatch``
  (`Pineau et al. 2017 <https://ui.adsabs.harvard.edu/abs/2017A%26C....1...51P>`_).
* `Taylor 2006 <https://arxiv.org/abs/astro-ph/0607353>`_ (ADASS XV):
  the STILTS tool that drives the queries.
