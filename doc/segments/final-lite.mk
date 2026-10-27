%_lite.fits: %.fits
	stilts tpipe in=$^ out=$@ cmd='delcols "adflux*_LS *_LS LU_flux*_LS *_GALEX *_VHS *_UKIDSS *_ALLWISE *_GALEXUL *_HSC_extaper *_HSC psf*_SDSS aper*_SDSS id*_SDSS ra*_SDSS dec*_SDSS"'
	stilts tpipe in=$@ omode=stats
