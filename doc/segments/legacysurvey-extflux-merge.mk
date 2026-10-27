%_all_extflux.fits: %_all.fits %_LS_extaper.fits
	# use aperture photometry until 1.5'', after that the PSF wings are not Gaussian anymore
	stilts tmatch2 out=$@ matcher=exact find=best join=all1 fixcols=dups \
		in1=$*_all.fits values1=id suffix1= \
		in2=$*_LS_extaper.fits values2=id  icmd2='keepcols "id apfluxext*_4 MW_TRANSMISSION_*"' suffix2=_LS_extaper \
		ocmd='addcol usable_extflux_g "(TYPE_LS!=\"PSF\"&&fracflux_g_LS<0.1&&decam_g_err>0&&apfluxext_err_g_4>0&&(apfluxext_g_4/MW_TRANSMISSION_G)*pow(10, -28.44)*1e26>decam_g*0.01)"' \
		ocmd='addcol usable_extflux_r "(TYPE_LS!=\"PSF\"&&fracflux_r_LS<0.1&&decam_r_err>0&&apfluxext_err_r_4>0&&(apfluxext_r_4/MW_TRANSMISSION_R)*pow(10, -28.44)*1e26>decam_r*0.01)"' \
		ocmd='addcol usable_extflux_z "(TYPE_LS!=\"PSF\"&&fracflux_z_LS<0.1&&decam_z_err>0&&apfluxext_err_z_4>0&&(apfluxext_z_4/MW_TRANSMISSION_Z)*pow(10, -28.44)*1e26>decam_z*0.01)"' \
		ocmd='addcol prior_GALflux_decam_g "usable_extflux_g?(apfluxext_g_4/MW_TRANSMISSION_G)*pow(10, -28.44)*1e26:-99"' \
		ocmd='addcol prior_GALflux_decam_r "usable_extflux_r?(apfluxext_r_4/MW_TRANSMISSION_R)*pow(10, -28.44)*1e26:-99"' \
		ocmd='addcol prior_GALflux_decam_z "usable_extflux_z?(apfluxext_z_4/MW_TRANSMISSION_Z)*pow(10, -28.44)*1e26:-99"' \
		ocmd='addcol prior_GALflux_decam_g_errlo "usable_extflux_g?(apfluxext_err_g_4/MW_TRANSMISSION_G)*pow(10, -28.44)*1e26:-99"' \
		ocmd='addcol prior_GALflux_decam_r_errlo "usable_extflux_r?(apfluxext_err_r_4/MW_TRANSMISSION_R)*pow(10, -28.44)*1e26:-99"' \
		ocmd='addcol prior_GALflux_decam_z_errlo "usable_extflux_z?(apfluxext_err_z_4/MW_TRANSMISSION_Z)*pow(10, -28.44)*1e26:-99"' \
		ocmd='addcol prior_GALflux_decam_g_errhi "1e10"' \
		ocmd='addcol prior_GALflux_decam_r_errhi "1e10"' \
		ocmd='addcol prior_GALflux_decam_z_errhi "1e10"' \
		ocmd='delcols "apfluxext* usable_extflux* MW_TRANSMISSION_* *_LS_extaper"'
