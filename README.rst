RainbowLasso
============

RainbowLasso compiles matched aperture fluxes from ultraviolet to infrared
all-sky astronomical surveys of stars, galaxies and Active Galactic Nuclei.

.. image:: logo.png

How it works
------------

Matched aperture flux means that the same extraction region is used in 
each wavelength, which includes 

* GALEX (far-UV and near-UV), 
* optical (griz bands from Legacy Survey), 
* near-infrared (UltraVISTA hemisphere survey, UKIDSS)
* infrared (WISE, decomposed by Legacy Survey)

Depending on whether the source is extended or point-like 
(information from the Legacy Survey multi-wavelength source fitting from optical to mid-infrared),
either PSF-like photometry or a total flux based on an extraction aperture is used.
If the latter, the same aperture radius is used also for the UV and near-infrared (matched aperture).
Only sources within the Legacy Survey area are handled.

This allows modelling a physically consistent region of a extended source.

RainbowLasso takes a input FITS catalog with coordinates (id, RA, DEC) and
automatically fetches the necessary photometry from publicly available sources (noirlab, vizier).
Additional input catalog columns are copied to the output and are not uploaded anywhere.
The output files are directly usable with LePhare (TODO), Cigale (TODO) and GRAHSP.

Prerequisites
-------------

You need to have the programs installed:

* `stilts <http://www.star.bristol.ac.uk/~mbt/stilts/sun256/sun256.html>`_
* make

and the following python libraries:

* requests_cache, astropy, numpy, pandas, tqdm
* dl (from `noirlab <https://datalab.noirlab.edu/docs/manual/UsingAstroDataLab/InstallDatalab/InstallDatalab/InstallDatalab.html>`_)

which you could get with something like::

	pip install requests_cache astropy numpy pandas tqdm noirlab


Usage
-----

If for example the input file is dr16QWX_selection.fits, then run::

	make dr16QWX_selection_lite.fits

Some steps are manual (instructions are shown).

The Makefile contains the steps performed to fetch the various multiwavelength surveys.

Checking photometry
--------------------

Create a visualisation of the errors and fluxes::

	make dr16QWX_selection_lite.fits_errors.pdf dr16QWX_selection_lite.fits_fluxes.pdf

How to read the diagnostic plot:

* each page is a filter. 
* fluxes.pdf compares the flux (mJy) from one band to the next. If the spectrum is relatively flat, and bands are close together, they should follow the 1:1 line.
* errors.pdf compares the flux error to the significance of the value (value / error ratio). 
* The previous filter is shown in gray. Within one telescope and filter system, the error properties should be comparable.

Contributors
------------

* Suraj D Shankar
* Mara Salvato
* Johannes Buchner
* Isabelle Gauger

Citing
------

1. Cite the accompaning paper. 
2. You can also include the repository URL as a footnote.
3. Cite the relevant data products. See the accompaning paper for a list of references.

Licence
-------

AGPL-3 (see LICENCE file).

Logo
-------

The logo is based on work by Ivan Abirawa, Those Icons and Freepik.
