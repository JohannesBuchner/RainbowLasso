RainbowLasso
============

RainbowLasso compiles reliable ultraviolet to infrared for 
modelling stars, galaxies and Active Galactic Nuclei
from public all-sky astronomical surveys.

.. image:: logo.png

How it works
------------

Understanding the stellar population in galaxies requires constraining 
a model with reliable multi-wavelength (ultraviolet to infrared) 
flux measurements. Since galaxies are extended, to model the same object 
at all wavelength, the same aperture needs to be extracted. This is called 
matched aperture. For point sources (stars or distant quasars),
measuring the flux weighted by the survey point-spread-function provides
the best constraints.

Depending on whether the source is extended or point-like, 
information from the Legacy Survey multi-wavelength source fitting from optical to mid-infrared,
either PSF-like photometry or a total flux based on an extraction aperture is used.
Only sources within the Legacy Survey area are handled.

RainbowLasso fetches and calculates fluxes from these surveys:

* GALEX (far-UV and near-UV), 
* optical (griz bands from Legacy Survey), 
* near-infrared (JKH from VISTA hemisphere survey and UKIDSS)
* infrared (WISE W1,W2,W3,W4, deblended by Legacy Survey)

RainbowLasso takes a input FITS catalog with coordinates (id, RA, DEC) and
automatically fetches the necessary photometry from publicly available sources (noirlab, vizier).
Additional input catalog columns are copied to the output and are not uploaded anywhere.
The output files should be usable with `GRAHSP <https://arxiv.org/abs/2405.19297>`_, Cigale, LePHARE (TODO).

Features
--------

* Milky Way attenuation by dust is corrected for each band up until the near-infrared. The E(B-V) is taken from `here <https://www.legacysurvey.org/dr10/catalogs/#galactic-extinction-coefficients>`_ and for other wavelengths as listed to the COSMOS2015 Laigle+ paper.
* Reliable flux errors are computed. If needed, they are conservatively estimated from magnitude errors.
* Unreliable flux estimates are discarded, indicated by flags on blending and image fitting problems
* Discards blended or problematic photometry in GALEX (artifact and extraction flag, magnitude error)
* Discards blended or problematic photometry in LS&WISE (fracin, fracflux)
* Provides upper limits for blended WISE cases
  * conservative 3 sigma upper limits are computed for WISE3&4 based on the total flux of all nearby sources (fracflux) and the error bars.
* Missing information (lacking information) is indicated (e.g. -99 entries).
* Automated
* For a clean sample, you can sub-select sources flagged "isolated"


Prerequisites
-------------

You need to have the programs installed:

* `stilts <http://www.star.bristol.ac.uk/~mbt/stilts/sun256/sun256.html>`_
* make

and the following python libraries:

* requests_cache, astropy, numpy, pandas, tqdm, uncertainties, dustmaps
* astro-datalab (`see noirlab <https://datalab.noirlab.edu/docs/manual/UsingAstroDataLab/InstallDatalab/InstallDatalab/InstallDatalab.html>`_)

which you could get with something like::

	pip install requests_cache astropy numpy pandas tqdm astro-datalab uncertainties dustmaps

* dustmaps need to be `set up and SFD downloaded <https://github.com/gregreen/dustmaps/?tab=readme-ov-file#getting-the-data>`_

Usage
-----

Prepare an input fits file with columns 'id', 'RA' and 'DEC'. Additional columns are fine.

If for example the input file is dr16QWX_selection.fits, then run::

	make dr16QWX_selection_all_lite.fits

Files can also be outside the RainbowLasso directory.

In that case you can run from the RainbowLasso directory::

	$ make path/to/data/dr16QWX_selection_all_lite.fits

Or you can run from the data directory::

	$ make -C path/to/RainbowLasso/ $PWD/dr16QWX_selection_all_lite.fits

The Makefile contains the steps performed to fetch the various multiwavelength surveys.

Important output files:

 * $PWD/dr16QWX_selection_all_lite.fits : File ready for use with GRAHSP/CIGALE.
 * $PWD/dr16QWX_selection_all.fits : this contains all columns from the input catalogs.
 
If you get the following warning::

	WARNING: Using non-standard extended column convention
	WARNING: Other FITS software may not see columns 999-1188

Then _all.fits has over 1000 columns, and astropy will not be able to read the last columns.

If you also want HSC-WIDE fluxes, you can run:

  $ make path/to/data/dr16QWX_selection_allHSC_lite.fits

Quality control
---------------

Create a visualisation of the errors and fluxes::

	make dr16QWX_selection_all_lite.fits_errors.pdf dr16QWX_selection_all_lite.fits_fluxes.pdf

Examples are uploaded to this repository.

How to read the diagnostic plot:

* each page is a filter. 
* fluxes.pdf compares the flux (mJy) from one band to the next. If the spectrum is relatively flat, and bands are close together, they should follow the 1:1 line.
* errors.pdf compares the flux error to the significance of the value (value / error ratio). 
* The previous filter is shown in gray. Within one telescope and filter system, the error properties should be comparable.

Known issues
------------

* If none of the sources in the input catalog are in VHS/UKIDSS, the query retrieval returns with no file produced (pandas data frames cannot be concatenated), and the pipeline fails

  * workaround: insert some source coordinates from the provided example, then delete them after. (thanks to Pietro Baldini)

TODO
----

* ☐ discard blended or problematic photometry in VHS
* ☐ discard blended or problematic photometry in UKIDSS
* ☐ provide upper limits for VHS, UKIDSS (need coverage information for that)
* ☐ improve efficiency of noirlab fetching of aperture fluxes. They were not able yet to suggest a better working solution.

To add more surveys, contributions are welcome.

**Scope**: The current scope are large-area surveys (ten to ten-thousands of square degrees),
which provide state-of-the-art photometry to an already-existing selection of sources.
Cross-matching (see `NWAY <https://github.com/JohannesBuchner/nway/>`_) or image analyses 
are outside the scope of this project.

Contributors
------------

* Suraj D Shankar
* Mara Salvato
* Johannes Buchner
* Isabelle Gauger

Citing
------

1. Cite the `accompaning paper <https://arxiv.org/abs/2405.19297>`_. 
2. You can also include the repository URL as a footnote.
3. Cite the data products of the surveys used. See the accompaning paper for a list of references.

Licence
-------

AGPL-3 (see LICENCE file).

Logo
-------

The logo is based on work by Ivan Abirawa, Those Icons and Freepik.
