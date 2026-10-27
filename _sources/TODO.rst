.. _todo:

====
TODO
====

Findings from the code review that should be looked at:

* **NUV extinction inconsistency**: the ``A_lam_nuv`` coefficient
  (``A_v*(0.674 + 1.248/3.1)`` ≈ 3.34*EBmV) disagrees with the NUV
  magnitude dereddening (8.621*EBmV, Schlafly & Finkbeiner 2011). The
  FUV coefficient (8.10*EBmV) is consistent with its 8.286*EBmV. This
  propagates to ``Nflux_real_LU`` (the final NUV flux, corrected by
  ``10^(1.34*EBmV)`` instead of ``~10^(3.45*EBmV)``: an NUV flux
  deficit of ~10% at EBmV=0.025, ~60% at EBmV=0.1) and to the
  ``NUV_fluxlim`` scaling factor 1.35 (= 0.4*3.34). Decide which
  extinction law is intended and make the magnitude and flux
  corrections consistent.
* **fitbits FORCED_POINTSOURCE**: for point-like sources, the
  ``goodfitsLS`` cut excludes fitbits bits 1|4|8 =
  FORCED_POINTSOURCE | HIT_RADIUS_LIMIT | HIT_SERSIC_LIMIT. The radius
  and sersic limit bits are clearly problematic, but FORCED_POINTSOURCE
  only means the source was fit as a point source from the reference
  catalog - normal for quasar candidates. Including it may discard
  valid point sources; check whether bit 1 should stay.
* **LS WISE aperture comment**: the comment says the extended-source
  aperture is "3.5'' radius = 7'' diameter ... IR:apflux_*_1", but the
  LS DR10 WISE apertures have radii 3, 5, 7, 9, 11 arcsec, so
  ``apflux_w1_1`` is 3.0'' radius (6.0'' diameter). The code is fine;
  the comment could be corrected.
* **UKIDSS aperture sizes**: the VHS aperture sizes are confirmed by
  the VizieR II/367 column descriptions (apermag4 = 2.8'' aperture
  diameter, apermag6 = 5.7'' diameter). UKIDSS uses the same CASU
  aperture sequence (radii 0.5, 0.71, 1.0, 1.4, 2.0, 2.8, 4.0, 5.66
  arcsec, steps of sqrt(2); AperMag3 = 1.0'' radius is the UKIDSS
  default aperture), so apermag4/ap6 should be 2.8''/5.7'' diameter
  there too, but the WSA schema pages could not be reached to confirm
  it - worth double-checking.
* **VHS error column pairing**: ``UV_Yap6flux_err`` uses
  ``Yapermag6err`` as the error of the ``apermagnoapercorr6``
  magnitude. Verify that ``apermag6err`` is the error of the
  no-aperture-correction aperture in VHS DR5 (and not only of the
  aperture-corrected one).
* **Dangling Makefile rule**: ``%_LS_extaper_in.fits`` depends on
  ``%_LS_staraper.fits`` and ``%_LS_aper.fits``, which are commented
  out (the rule is unused in the current flow and would fail if built).
* **Redundant columns**: ``UV_Yapc6flux`` duplicates ``UV_Yap6flux``
  (same formula, unused in the final merge); the ``z_usable`` check in
  the HSC extflux rule tests ``z_inputcount_value>0`` twice (harmless).
* Discarding blended or problematic photometry in VHS and UKIDSS, and
  upper limits for them (need coverage information), as listed in the
  README.

The bugs already fixed on this branch: the swapped GALEX FUV/NUV depth
cuts (FUV=19.9, NUV=20.8 ABmag), and the swapped HSC point/extended
flux assignment (point sources use the PSF model flux, extended sources
the 5.75'' radius aperture).
