%_SDSS.fits: fetchSDSS.py %.fits
	python3 $^ $@ || stilts tpipe in=$^ out=$@ cmd='keepcols "id"'
