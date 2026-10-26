.. _install:
.. highlight:: shell

============
Installation
============

RainbowLasso is driven by a Makefile, and needs a few external programs
and python packages.

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

Example: creating a conda environment
-------------------------------------

With conda (or mamba), a working environment can be created with::

    mamba create -n rainbowlasso -c conda-forge python=3.12 \
        numpy matplotlib joblib astropy requests requests_cache scipy \
        tqdm uncertainties pillow dustmaps make pip
    pip install astro-datalab

Downloading RainbowLasso
------------------------

Get RainbowLasso onto your computer::

    git clone git@github.com:JohannesBuchner/RainbowLasso.git
    cd RainbowLasso

Update RainbowLasso to the latest version::

    git pull

This directory contains several scripts, the Makefile and this
documentation (in the doc/ folder).

Logging in to the services
--------------------------

The pipeline queries two public catalog services (see :doc:`noirlab` and
:doc:`vizier`):

* Data Lab at noirlab requires a free account. Log in to datalab::

    datalab login

* The VizieR/CDS X-Match service (used for GALEX and AllWISE) needs no
  login.

Verifying the installation
--------------------------

Run the test pipeline on the small example catalog (50 sources)::

    make test

This runs ``make dr16QWX_selection_all_lite.fits`` end-to-end, including
the Data Lab and VizieR queries. The first run takes a few minutes;
results are cached by ``requests_cache`` and joblib.

.. _makefilesegments:

Makefile segments for the documentation
---------------------------------------

Each documentation page for a survey/catalog includes the relevant
Makefile rule (target + recipe) via ``literalinclude``. The script
``makefilesegments.py`` extracts these rule segments from the Makefile
into the ``doc/_segments/`` folder:

.. code-block:: console

    $ python3 makefilesegments.py            # extract all segments
    $ python3 makefilesegments.py galex      # extract one segment
    $ python3 makefilesegments.py --list     # list available segments

The segments are regenerated whenever the Makefile changes; the ``docs``
Makefile target runs the extraction and builds this documentation::

    make docs
