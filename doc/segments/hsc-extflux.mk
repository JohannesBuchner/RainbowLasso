%_HSCall_extflux.fits: %_HSCall.fits %_HSC_extaper.fits
	# first, we should get the neighbours to know how contaminated the image is, so we can know how far from the center we can safely integrate
	#usable = True
	# obj['neighbours_' + band + '_saturated'] = not (np.all(nearby_neighbours[band + '_pixelflags_saturated'] == 0) and np.all(nearby_neighbours[band + '_pixelflags_saturatedcenter'] == 0))
	#if not np.all(nearby_neighbours[band + '_pixelflags_saturated'] == 0 and nearby_neighbours[band + '_pixelflags_saturatedcenter'] == 0):
	#    usable = False
	# obj['neighbours_' + band + '_flux_30'] = verynearby_neighbours[band + '_psfflux_flux'].sum()
	#if not verynearby_neighbours[band + '_cmodel_flux'].sum() > obj[band + '_cmodel_flux']:
	#    print("skipping: neighbouring objects within 3'' are brighter", obj['id'])
	#    diagnostic_str += "closeby objects outshine"
	#    usable = False
	# obj['neighbours_' + band + '_flux_50'] = nearby_neighbours[band + '_psfflux_flux'].sum()
	#if not nearby_neighbours[band + '_cmodel_flux'].sum() > 10 * obj[band + '_cmodel_flux']:
	#    print("skipping: neighbouring objects within 5'' extremely bright", obj['id'])
	#    diagnostic_str += "neighbours extremely bright"
	#    usable = False
	#aperture_radius = 40
	#if not nearby_neighbours[band + '_cmodel_flux'].sum() > obj[band + '_cmodel_flux']:
	#    print("neighbouring objects are bright, only measuring until 2''", obj['id'])
	#    diagnostic_str += "neighbouring objects are bright, using 2arcsec"
	#    aperture_radius = 20
	# also the object should not be saturated
	stilts tmatch2 out=$@ matcher=sky find=best params=1 fixcols=all \
		in1=$*_HSCall.fits values1="RA DEC" suffix1= \
		in2=$*_HSC_extaper.fits values2="RA DEC" suffix2=_HSC_extaper \
		icmd2='addcol g_blended "neighbours_g_saturated||neighbours_g_flux_50 > 10*g_psfflux_flux"' \
		icmd2='addcol r_blended "neighbours_r_saturated||neighbours_r_flux_50 > 10*r_psfflux_flux"' \
		icmd2='addcol i_blended "neighbours_i_saturated||neighbours_i_flux_50 > 10*i_psfflux_flux"' \
		icmd2='addcol z_blended "neighbours_z_saturated||neighbours_z_flux_50 > 10*z_psfflux_flux"' \
		icmd2='addcol y_blended "neighbours_y_saturated||neighbours_y_flux_50 > 10*y_psfflux_flux"' \
		icmd2='addcol g_usable "!z_pixelflags_saturated&&!g_pixelflags_saturated&&!g_pixelflags_saturatedcenter&&g_inputcount_value>0&&!g_inputcount_flag&&!g_localbackground_flag&&!g_pixelflags&&!g_apertureflux_40_flag&&apfluxext_err_g_40>0&&g_extendedness_value>0&&!g_extendedness_flag"' \
		icmd2='addcol r_usable "!z_pixelflags_saturated&&!r_pixelflags_saturated&&!r_pixelflags_saturatedcenter&&r_inputcount_value>0&&!r_inputcount_flag&&!r_localbackground_flag&&!r_pixelflags&&!r_apertureflux_40_flag&&apfluxext_err_r_40>0&&r_extendedness_value>0&&!r_extendedness_flag"' \
		icmd2='addcol i_usable "!z_pixelflags_saturated&&!i_pixelflags_saturated&&!i_pixelflags_saturatedcenter&&i_inputcount_value>0&&!i_inputcount_flag&&!i_localbackground_flag&&!i_pixelflags&&!i_apertureflux_40_flag&&apfluxext_err_i_40>0&&i_extendedness_value>0&&!i_extendedness_flag"' \
		icmd2='addcol z_usable "!z_pixelflags_saturated&&!z_pixelflags_saturatedcenter&&z_inputcount_value>0&&z_inputcount_value>0&&!z_inputcount_flag&&!z_localbackground_flag&&!z_pixelflags&&!z_apertureflux_40_flag&&apfluxext_err_z_40>0&&z_extendedness_value>0&&!z_extendedness_flag"' \
		icmd2='addcol y_usable "!z_pixelflags_saturated&&!y_pixelflags_saturated&&!y_pixelflags_saturatedcenter&&y_inputcount_value>0&&!y_inputcount_flag&&!y_localbackground_flag&&!y_pixelflags&&!y_apertureflux_40_flag&&apfluxext_err_y_40>0&&y_extendedness_value>0&&!y_extendedness_flag"' \
		ocmd='delcols "*_apertureflux_*HSC* *_pixelflags*HSC* *_kron*HSC* MW_TRANSMISSION_*"' \
		ocmd='addcol prior_GALflux_HSC_g "(!g_blended_HSC_extaper&&g_usable_HSC_extaper)?(apfluxext_g_40_HSC_extaper/pow(10, -0.4*a_g_HSC_extaper))/1000*pow(10,-29)*1e26:-99"' \
		ocmd='addcol prior_GALflux_HSC_r "(!r_blended_HSC_extaper&&r_usable_HSC_extaper)?(apfluxext_r_40_HSC_extaper/pow(10, -0.4*a_r_HSC_extaper))/1000*pow(10,-29)*1e26:-99"' \
		ocmd='addcol prior_GALflux_HSC_i "(!i_blended_HSC_extaper&&i_usable_HSC_extaper)?(apfluxext_i_40_HSC_extaper/pow(10, -0.4*a_i_HSC_extaper))/1000*pow(10,-29)*1e26:-99"' \
		ocmd='addcol prior_GALflux_HSC_z "(!z_blended_HSC_extaper&&z_usable_HSC_extaper)?(apfluxext_z_40_HSC_extaper/pow(10, -0.4*a_z_HSC_extaper))/1000*pow(10,-29)*1e26:-99"' \
		ocmd='addcol prior_GALflux_HSC_y "(!y_blended_HSC_extaper&&y_usable_HSC_extaper)?(apfluxext_y_40_HSC_extaper/pow(10, -0.4*a_y_HSC_extaper))/1000*pow(10,-29)*1e26:-99"' \
		ocmd='addcol prior_GALflux_HSC_g_errlo "(!g_blended_HSC_extaper&&g_usable_HSC_extaper&&fracflux_g_LS<0.1)?(apfluxext_err_g_40_HSC_extaper/pow(10, -0.4*a_g_HSC_extaper))/1000*pow(10,-29)*1e26:-99"' \
		ocmd='addcol prior_GALflux_HSC_r_errlo "(!r_blended_HSC_extaper&&r_usable_HSC_extaper&&fracflux_r_LS<0.1)?(apfluxext_err_r_40_HSC_extaper/pow(10, -0.4*a_r_HSC_extaper))/1000*pow(10,-29)*1e26:-99"' \
		ocmd='addcol prior_GALflux_HSC_i_errlo "(!i_blended_HSC_extaper&&i_usable_HSC_extaper&&fracflux_r_LS<0.1)?(apfluxext_err_i_40_HSC_extaper/pow(10, -0.4*a_i_HSC_extaper))/1000*pow(10,-29)*1e26:-99"' \
		ocmd='addcol prior_GALflux_HSC_z_errlo "(!z_blended_HSC_extaper&&z_usable_HSC_extaper&&fracflux_z_LS<0.1)?(apfluxext_err_z_40_HSC_extaper/pow(10, -0.4*a_z_HSC_extaper))/1000*pow(10,-29)*1e26:-99"' \
		ocmd='addcol prior_GALflux_HSC_y_errlo "(!y_blended_HSC_extaper&&y_usable_HSC_extaper&&fracflux_z_LS<0.1)?(apfluxext_err_y_40_HSC_extaper/pow(10, -0.4*a_y_HSC_extaper))/1000*pow(10,-29)*1e26:-99"' \
		ocmd='addcol prior_GALflux_HSC_g_errhi "1e10"' \
		ocmd='addcol prior_GALflux_HSC_r_errhi "1e10"' \
		ocmd='addcol prior_GALflux_HSC_i_errhi "1e10"' \
		ocmd='addcol prior_GALflux_HSC_z_errhi "1e10"' \
		ocmd='addcol prior_GALflux_HSC_y_errhi "1e10"' \
		ocmd='delcols "apfluxext_*HSC*"' \
