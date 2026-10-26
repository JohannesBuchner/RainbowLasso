%_withSDSS.fits: %.fits %_SDSS.fits
	stilts tmatch2 out=$@ join=all1 \
		in1=$*.fits suffix1= values1=id \
		in2=$*_SDSS.fits suffix2=_SDSS values2=id \
		fixcols=all matcher=exact \
		ocmd='addcol u_sdss "pointlike?psfMag_u_SDSS:aper_u_SDSS"' \
		ocmd='addcol u_sdss_err "pointlike?psfMagErr_u_SDSS:aper_u_err_SDSS"' \
		ocmd='addcol g_sdss "pointlike?psfMag_g_SDSS:aper_g_SDSS"' \
		ocmd='addcol g_sdss_err "pointlike?psfMagErr_g_SDSS:aper_g_err_SDSS"' \
		ocmd='addcol r_sdss "pointlike?psfMag_r_SDSS:aper_r_SDSS"' \
		ocmd='addcol r_sdss_err "pointlike?psfMagErr_r_SDSS:aper_r_err_SDSS"' \
		ocmd='addcol i_sdss "pointlike?psfMag_i_SDSS:aper_i_SDSS"' \
		ocmd='addcol i_sdss_err "pointlike?psfMagErr_i_SDSS:aper_i_err_SDSS"' \
		ocmd='addcol z_sdss "pointlike?psfMag_z_SDSS:aper_z_SDSS"' \
		ocmd='addcol z_sdss_err "pointlike?psfMagErr_z_SDSS:aper_z_err_SDSS"' \
