.. _outputs:

=============
Output files
=============

The pipeline produces a series of FITS files next to the input catalog.
For an input ``myfilename.fits``, the main products are:

.. list-table::
   :header-rows: 1

   * - file
     - contents
   * - ``myfilename_all_lite.fits``
     - reduced merged catalog, ready for use with GRAHSP/CIGALE
   * - ``myfilename_all.fits``
     - merged catalog with all columns from the input catalogs
   * - ``myfilename_all_extflux.fits``
     - merged catalog with galaxy total-flux priors (``prior_GALflux_decam_*``) from LS apertures
   * - ``myfilename_all_extflux_withEuclid.fits``
     - with the priors overridden by Euclid Q1 constraints
   * - ``myfilename_HSCall.fits``
     - merged catalog using HSC-WIDE (requires both LS and HSC detections)
   * - ``myfilename_HSCall_extflux.fits``
     - with galaxy total-flux priors from HSC apertures
   * - ``myfilename_withSDSS.fits``
     - merged catalog with SDSS magnitudes

Additionally, intermediate files are kept (``.SECONDARY`` in the
Makefile): ``myfilename_coords.csv``, ``myfilename_GALEX.fits``,
``myfilename_GALEX_UL.fits``, ``myfilename_LS_noirlab.fits``,
``myfilename_LS.fits``, ``myfilename_UKIDSS_noirlab.fits``,
``myfilename_UKIDSS.fits``, ``myfilename_VHS_noirlab.fits``,
``myfilename_VHS.fits``, ``myfilename_ALLWISE_all.fits``,
``myfilename_ALLWISE_sum.fits``. They document the per-catalog cuts and
can be deleted to save space.

Format
------

All merged outputs are FITS binary tables with **one row per input
source**. The input columns (``id``, ``RA``, ``DEC``, and any additional
input catalog columns) are copied to the output. Only the id, RA and DEC
columns are uploaded to the external services; additional input columns
are not uploaded anywhere.

Merging is done by exact id match (``stilts tmatchn``), so missing
surveys simply leave the respective columns at the missing value.

Column convention
-----------------

Flux columns
~~~~~~~~~~~~

* All final flux columns are in **mJy** (1e-26 erg/s/cm^2/Hz), e.g.
  ``decam_g``, ``UV_J``, ``WFCAM_Ks``, ``WISE1``, ``FUV``, ``HSC_g``.
* Every flux column ``X`` has an error column ``X_err`` in the same
  unit.
* Intermediate per-catalog columns use the internal **LU** unit
  (erg/s/cm^2/Hz; the name goes back to lephare units), e.g.
  ``LU_flux_g_LS``, ``UV_Yapc4flux_VHS``, ``Fflux_real_LU_GALEX``.
* **Missing** information is indicated with the value **-99**.
* **Upper limits**: a **negative error value** indicates that the value
  is a (3 sigma) upper limit, not a detection. For example, a source
  without GALEX detection gets ``NUV = 3 * NUV_fluxlim`` and
  ``NUV_err = -NUV``; the AllWISE fallback and the WISE3/4 upper limits
  use the same convention. Missing values (-99) can carry either sign
  and should be discarded first.

Band naming
~~~~~~~~~~~

Each telescope/instrument contributes columns named
``<prefix>_<band>``:

.. list-table::
   :header-rows: 1

   * - prefix
     - telescope / instrument
     - bands
   * - ``decam_``
     - DECam (Legacy Survey)
     - g, r, i, z
   * - ``90prime_``
     - 90prime (BASS)
     - g, r
   * - ``zd_mosaic``
     - MzLS
     - z (no ``_<band>`` suffix)
   * - ``UV_``
     - VISTA (VHS)
     - Y, J, H, K
   * - ``WFCAM_``
     - WFCAM (UKIDSS)
     - Y, J, H, Ks
   * - ``WISE1..4``
     - WISE (LS deblending or AllWISE)
     - W1, W2, W3, W4
   * - ``FUV``, ``NUV``
     - GALEX AIS
     - FUV, NUV
   * - ``HSC_``
     - Hyper Suprime-Cam (PDR3 WIDE)
     - g, r, i, z, y

The sky region decides between DECam and BASS/MzLS: sources in the
BASS/MzLS area (DEC > 32, 90 < RA < 300, column ``inMzLSBASS``) get
``90prime_g``, ``90prime_r`` and ``zd_mosaic`` instead of ``decam_g``,
``decam_r``, ``decam_z``.

Diagnostic columns
~~~~~~~~~~~~~~~~~~

* ``pointlike``: whether the source is treated as a point source
  (PSF-like photometry) or extended (aperture photometry).
* ``goodfitsLS``, ``goodfitsg``, ``goodfitsr``, ``goodfitsz``,
  ``goodfitsi``, ``goodfitsW1``, ``goodfitsW2``: whether the band
  extraction is trustworthy (blending and image fitting problems).
* ``isolatedLS``: whether the source is isolated (no bright neighbours,
  from the LS fracflux columns). For a clean sample, sub-select sources
  with ``isolatedLS == True``.
* ``W34_blended``, ``WISE1_origin``, ``WISE2_origin``, ``WISE34_origin``:
  the origin of the WISE values (``LS10`` = LS-deblended measurement,
  ``AllWISE`` = summed AllWISE flux of the neighbours, ``LS10UL`` =
  conservative 3 sigma upper limit).
* ``prior_GALflux_*``: galaxy total-flux priors (mJy, with ``_errlo``
  and ``_errhi`` bounds) for the extended-flux products.

Reading large files
-------------------

If you get the following warning::

    WARNING: Using non-standard extended column convention
    WARNING: Other FITS software may not see columns 999-1188

then ``_all.fits`` has over 1000 columns, and astropy will not be able
to read the last columns. Use ``_all_lite.fits`` instead, which only
contains the final flux columns.
