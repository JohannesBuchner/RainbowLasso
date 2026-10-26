.. _euclid:

================================
Euclid Q1 (optional)
================================

*Catalog*: `Euclid Q1 <https://www.euclid-ec.org/public/earlyrelease/>`_,
queried through the NOAO Data Lab (see :doc:`noirlab`):
``euclid_q1.object``.

*Wavelengths*: VIS (5400-9300 A), NISP (Y, J, H).

Matching
--------

.. literalinclude:: _segments/euclid-fetch.mk
   :language: make

The nearest source within 1 arcsec is retrieved (SQL ``q3c_join`` via
Data Lab).

Extended-source flux priors
---------------------------

Currently, Euclid catalogs do not provide aperture fluxes with suitable
diameters. We follow other studies and assume the morphology is
~wavelength-independent, so the Euclid extended flux fraction is
applied to the optical flux columns as a prior:

.. literalinclude:: _segments/euclid-prior.mk
   :language: make

* **Reliability** (``DET_QUALITY_FLAG < 4 && SPURIOUS_FLAG == 0``):
  only reliable detections are used.
* **Extended fraction** (``frac_extended_euclid``): estimated from the
  surface brightness, ``min(10^(0.4*(mumax_minus_mag)) - 0.06 - 0.03,
  1)``.
* Euclid VIS (5400-9300 A) overlaps with the DECam g, r, i bands (not
  z), so ``prior_GALflux_decam_g`` and ``prior_GALflux_decam_r`` are
  overridden if the Euclid constraint is better; for the error a 1%
  systematic uncertainty is assumed.

Output columns
--------------

``prior_GALflux_decam_g`` and ``prior_GALflux_decam_r`` (with
``_errlo``, ``_errhi``) in the ``%_all_extflux_withEuclid.fits`` file.
