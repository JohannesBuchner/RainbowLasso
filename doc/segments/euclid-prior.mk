%_all_extflux_withEuclid.fits: %_all_extflux.fits %_EuclidQ1_noirlab.fits
	# currently, Euclid catalogs do not provide aperture fluxes with suitable diameters
	# We follow other studies and assume the morphology ~wavelength-independent
	# so we first compute the Euclid extended flux, then apply it to optical flux columns.
	# Euclid VIS is 5400A-9300A, which overlaps with decam g,r,i (not z)
	# we add a gal_decam_r+z prior based on the fraction of Euclid extended flux.
	# min(pow(10, 0.4*(mumax_minus_mag)) - 0.06 - 0.03, 1)
	# for the error, we assume 1% 1sigma error (for systematics).
	# we override the decam columns if the Euclid constraints are better (they should be).
	stilts tmatch2 out=$@ join=all1 \
		in1=$< suffix1= values1=id \
		in2=$*_EuclidQ1_noirlab.fits suffix2=_EuclidQ1 values2=id \
		fixcols=all matcher=exact \
		ocmd='addcol reliable_euclid "DET_QUALITY_FLAG_EuclidQ1 < 4 && SPURIOUS_FLAG_EuclidQ1==0"' \
		ocmd='addcol frac_extended_euclid "reliable_euclid ? min(pow(10, 0.4*(mumax_minus_mag_EuclidQ1)) - 0.06 - 0.03, 1):-99"' \
		ocmd='addcol prior_GALflux_decam_g_keep "!(frac_extended_euclid>0)||prior_GALflux_decam_g_errlo>0&&prior_GALflux_decam_g>decam_g*frac_extended_euclid"' \
		ocmd='addcol prior_GALflux_decam_r_keep "!(frac_extended_euclid>0)||prior_GALflux_decam_r_errlo>0&&prior_GALflux_decam_r>decam_r*frac_extended_euclid"' \
		ocmd='replacecol prior_GALflux_decam_g "prior_GALflux_decam_g_keep?prior_GALflux_decam_g:(decam_g>0&&frac_extended_euclid>0?decam_g*frac_extended_euclid:-99)"' \
		ocmd='replacecol prior_GALflux_decam_r "prior_GALflux_decam_r_keep?prior_GALflux_decam_r:(decam_r>0&&frac_extended_euclid>0?decam_r*frac_extended_euclid:-99)"' \
		ocmd='replacecol prior_GALflux_decam_g_errlo "prior_GALflux_decam_g_keep?prior_GALflux_decam_g_errlo:(prior_GALflux_decam_g>0?prior_GALflux_decam_g*0.01:-99)"' \
		ocmd='replacecol prior_GALflux_decam_r_errlo "prior_GALflux_decam_r_keep?prior_GALflux_decam_r_errlo:(prior_GALflux_decam_r>0?prior_GALflux_decam_r*0.01:-99)"' \
		ocmd='replacecol prior_GALflux_decam_g_errhi "prior_GALflux_decam_g_keep?prior_GALflux_decam_g_errhi:1e10"' \
		ocmd='replacecol prior_GALflux_decam_r_errhi "prior_GALflux_decam_r_keep?prior_GALflux_decam_r_errhi:1e10"' \
