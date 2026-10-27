%_HSC.fits: %_HSC_aper.fits %.fits
	# adflux_* = adaptive flux column, depending on source type
	# for extended sources: 5.7'' aperture = g_apertureflux_57_flux
	# for point sources: model flux
	# Milky way extinction correction
	# 28.44 is to convert to erg/s/cm^2/Hz
	stilts tmatch2 out=$@ matcher=sky find=best params=1 fixcols=all \
		in1=$*.fits values1="RA DEC" icmd1='keepcols "id RA DEC"' suffix1=_in \
		in2=$< values2="RA DEC" suffix2= \
		ocmd='addcol extended "g_extendedness_value==1||r_extendedness_value==1||i_extendedness_value==1||z_extendedness_value==1||y_extendedness_value==1"' \
		ocmd='addcol adflux_g "extended ? g_apertureflux_57_flux : g_psfflux_flux"' \
		ocmd='addcol adflux_r "extended ? r_apertureflux_57_flux : r_psfflux_flux"' \
		ocmd='addcol adflux_i "extended ? i_apertureflux_57_flux : i_psfflux_flux"' \
		ocmd='addcol adflux_z "extended ? z_apertureflux_57_flux : z_psfflux_flux"' \
		ocmd='addcol adflux_y "extended ? y_apertureflux_57_flux : y_psfflux_flux"' \
		ocmd='addcol adflux_g_err "extended ? g_apertureflux_57_fluxerr : g_psfflux_fluxerr"' \
		ocmd='addcol adflux_r_err "extended ? r_apertureflux_57_fluxerr : r_psfflux_fluxerr"' \
		ocmd='addcol adflux_i_err "extended ? i_apertureflux_57_fluxerr : i_psfflux_fluxerr"' \
		ocmd='addcol adflux_z_err "extended ? z_apertureflux_57_fluxerr : z_psfflux_fluxerr"' \
		ocmd='addcol adflux_y_err "extended ? y_apertureflux_57_fluxerr : y_psfflux_fluxerr"' \
		ocmd='addcol LU_flux_g "ADFLUX_G>-10000000000L ? ((ADFLUX_G/1000*(pow(10, -(-a_g)/(2.5)))))*pow(10, -29):-99."' \
		ocmd='addcol LU_flux_r "ADFLUX_R>-10000000000L ? ((ADFLUX_R/1000*(pow(10, -(-a_r)/(2.5)))))*pow(10, -29):-99."' \
		ocmd='addcol LU_flux_i "ADFLUX_I>-10000000000L ? ((ADFLUX_I/1000*(pow(10, -(-a_i)/(2.5)))))*pow(10, -29):-99."' \
		ocmd='addcol LU_flux_z "ADFLUX_Z>-10000000000L ? ((ADFLUX_Z/1000*(pow(10, -(-a_z)/(2.5)))))*pow(10, -29):-99."' \
		ocmd='addcol LU_flux_y "ADFLUX_Y>-10000000000L ? ((ADFLUX_Y/1000*(pow(10, -(-a_y)/(2.5)))))*pow(10, -29):-99."' \
		ocmd='addcol LU_flux_g_err "ADFLUX_G_err>-10000000000L ? ((ADFLUX_G_err/1000*(pow(10, -(-a_g)/(2.5)))))*pow(10, -29):-99."' \
		ocmd='addcol LU_flux_r_err "ADFLUX_R_err>-10000000000L ? ((ADFLUX_R_err/1000*(pow(10, -(-a_r)/(2.5)))))*pow(10, -29):-99."' \
		ocmd='addcol LU_flux_i_err "ADFLUX_I_err>-10000000000L ? ((ADFLUX_I_err/1000*(pow(10, -(-a_i)/(2.5)))))*pow(10, -29):-99."' \
		ocmd='addcol LU_flux_z_err "ADFLUX_Z_err>-10000000000L ? ((ADFLUX_Z_err/1000*(pow(10, -(-a_z)/(2.5)))))*pow(10, -29):-99."' \
		ocmd='addcol LU_flux_y_err "ADFLUX_Y_err>-10000000000L ? ((ADFLUX_Y_err/1000*(pow(10, -(-a_y)/(2.5)))))*pow(10, -29):-99."' \
		ocmd='delcols "RA_in DEC_in"; colmeta -name id id_in'
