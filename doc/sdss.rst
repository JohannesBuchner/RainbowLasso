.. _sdss:

================================
SDSS (optional)
================================

*Catalog*: `SDSS DR16 <https://www.sdss.org/dr16/>`_ (PhotoObj), queried
through SciServer CasJobs (requires a sciserver account).

*Wavelengths*: u, g, r, i, z.

Setup
-----

To use SDSS fluxes, put a sciserver username and password into a text
file ``~/.config/sciserver/login.txt``::

    username:password

Matching and fluxes
-------------------

.. literalinclude:: _segments/sdss-fetch.mk
   :language: make

The nearest source within 1 arcsec is retrieved
(``fetchSDSS.py``). The photometry flags are decoded with the SDSS
flag names, and the usable flags are computed (see ``fetchSDSS.py``).

In the final ``%_withSDSS.fits`` file:

.. literalinclude:: _segments/sdss-fluxes.mk
   :language: make

the flux extraction is adapted to the source type:

* for **point sources**: the PSF magnitudes (``psfMag_*``),
* for **extended sources**: the aperture magnitudes (``aper_*``).

Output columns
--------------

``u_sdss``, ``g_sdss`` .. (and ``_err``): magnitudes and errors in the
``%_withSDSS.fits`` file.
