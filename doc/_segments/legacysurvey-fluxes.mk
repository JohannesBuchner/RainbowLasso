%_LS.fits: %_LS_noirlab.fits
	# adflux_* = adaptive flux column, depending on source type
	# for extended sources: 3.5'' radius aperture = 7'' diameter aperture = OPT:apflux_*_6 and IR:apflux_*_1
	# for point sources: model flux
	# Milky way extinction correction
	# 28.44 is to convert to erg/s/cm^2/Hz
	stilts tpipe in=$^ out=$@ \
		cmd='addcol adflux_g "type==\"PSF\" ? flux_g : apflux_g_6"' \
		cmd='addcol adflux_r "type==\"PSF\" ? flux_r : apflux_r_6"' \
		cmd='addcol adflux_i "type==\"PSF\" ? flux_i : apflux_i_6"' \
		cmd='addcol adflux_z "type==\"PSF\" ? flux_z : apflux_z_6"' \
		cmd='addcol adflux_W1 "type==\"PSF\" ? flux_W1 : apflux_W1_1"' \
		cmd='addcol adflux_W2 "type==\"PSF\" ? flux_W2 : apflux_W2_1"' \
		cmd='addcol adflux_W3 "type==\"PSF\" ? flux_W3 : apflux_W3_1"' \
		cmd='addcol adflux_W4 "type==\"PSF\" ? flux_W4 : apflux_W4_1"' \
		cmd='addcol adflux_ivar_g  "type==\"PSF\" ? flux_ivar_g : apflux_ivar_g_6"' \
		cmd='addcol adflux_ivar_r  "type==\"PSF\" ? flux_ivar_r : apflux_ivar_r_6"' \
		cmd='addcol adflux_ivar_i  "type==\"PSF\" ? flux_ivar_i : apflux_ivar_i_6"' \
		cmd='addcol adflux_ivar_z  "type==\"PSF\" ? flux_ivar_z : apflux_ivar_z_6"' \
		cmd='addcol adflux_ivar_W1 "type==\"PSF\" ? flux_ivar_W1 : apflux_ivar_W1_1"' \
		cmd='addcol adflux_ivar_W2 "type==\"PSF\" ? flux_ivar_W2 : apflux_ivar_W2_1"' \
		cmd='addcol adflux_ivar_W3 "type==\"PSF\" ? flux_ivar_W3 : apflux_ivar_W3_1"' \
		cmd='addcol adflux_ivar_W4 "type==\"PSF\" ? flux_ivar_W4 : apflux_ivar_W4_1"' \
		cmd='addcol LU_flux_g "FLUX_G>-10000000000L ? ((ADFLUX_G/MW_TRANSMISSION_G))*pow(10, -28.44):-99."' \
		cmd='addcol LU_flux_r "FLUX_R>-10000000000L ? ((ADFLUX_R/MW_TRANSMISSION_R))*pow(10, -28.44):-99."' \
		cmd='addcol LU_flux_i "FLUX_I>-10000000000L ? ((ADFLUX_I/MW_TRANSMISSION_I))*pow(10, -28.44):-99."' \
		cmd='addcol LU_flux_z "FLUX_Z>-10000000000L ? ((ADFLUX_Z/MW_TRANSMISSION_Z))*pow(10, -28.44):-99."' \
		cmd='addcol LU_flux_w1 "FLUX_W1>-10000000000L ? ((ADFLUX_W1/MW_TRANSMISSION_W1))*pow(10, -28.44):-99."' \
		cmd='addcol LU_flux_w2 "FLUX_W2>-10000000000L ? ((ADFLUX_W2/MW_TRANSMISSION_W2))*pow(10, -28.44):-99."' \
		cmd='addcol LU_flux_w3 "FLUX_W3>-10000000000L ? ((ADFLUX_W3/MW_TRANSMISSION_W3))*pow(10, -28.44):-99."' \
		cmd='addcol LU_flux_w4 "FLUX_W4>-10000000000L ? ((ADFLUX_W4/MW_TRANSMISSION_W4))*pow(10, -28.44):-99."' \
		cmd='addcol LU_flux_g_err "ADFLUX_IVAR_G>0 ? (((1/pow(ADFLUX_IVAR_G, 0.5))/MW_TRANSMISSION_G))*pow(10, -28.44):-99."' \
		cmd='addcol LU_flux_r_err "ADFLUX_IVAR_R>0 ? (((1/pow(ADFLUX_IVAR_R, 0.5))/MW_TRANSMISSION_R))*pow(10, -28.44):-99."' \
		cmd='addcol LU_flux_i_err "ADFLUX_IVAR_I>0 ? (((1/pow(ADFLUX_IVAR_I, 0.5))/MW_TRANSMISSION_I))*pow(10, -28.44):-99."' \
		cmd='addcol LU_flux_z_err "ADFLUX_IVAR_Z>0 ? (((1/pow(ADFLUX_IVAR_Z, 0.5))/MW_TRANSMISSION_Z))*pow(10, -28.44):-99."' \
		cmd='addcol LU_flux_w1_err "ADFLUX_IVAR_W1>0 ? (((1/pow(ADFLUX_IVAR_W1, 0.5))/MW_TRANSMISSION_W1))*pow(10, -28.44):-99."' \
		cmd='addcol LU_flux_w2_err "ADFLUX_IVAR_W2>0 ? (((1/pow(ADFLUX_IVAR_W2, 0.5))/MW_TRANSMISSION_W2))*pow(10, -28.44):-99."' \
		cmd='addcol LU_flux_w3_err "ADFLUX_IVAR_W3>0 ? (((1/pow(ADFLUX_IVAR_W3, 0.5))/MW_TRANSMISSION_W3))*pow(10, -28.44):-99."' \
		cmd='addcol LU_flux_w4_err "ADFLUX_IVAR_W4>0 ? (((1/pow(ADFLUX_IVAR_W4, 0.5))/MW_TRANSMISSION_W4))*pow(10, -28.44):-99."'
