%_coords.csv: %.fits
	stilts tpipe in=$^ cmd='keepcols "ID RA DEC"' out=$@
