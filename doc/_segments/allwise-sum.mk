%_ALLWISE_sum.fits: groupsumWISE.py %_coords.csv %_ALLWISE_all.fits
	python3 $^ $@
