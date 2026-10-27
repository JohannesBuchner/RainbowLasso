%_LS_extaper.fits: addextflux.py %_LS_noirlab.fits
	python3 $^ $@
