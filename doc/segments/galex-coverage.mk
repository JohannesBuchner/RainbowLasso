%_GALEX_UL.fits: %_coords.csv galex_ais_ctrs_ebv.fits
	# search within 0.56 degrees
	stilts tmatch2 out=$@ matcher=sky find=best params=1980 \
		in1=$< values1="RA Dec" \
		in2=galex_ais_ctrs_ebv.fits values2="RAfdeg DEfdeg" \
		ocmd='addcol NUV_fluxlim "0.02 * pow(10, EBV_min * 1.35)"'
		#ocmd='addcol NUV_maglim "20.8 - 8.621 * EBV_min"' \
		#ocmd='addcol NUV_fluxlim "3.631 * pow(10, (NUV_maglim - 22.5) / 2.5)"'
