%_HSCall.fits: %.fits %_GALEX.fits %_LS.fits %_UKIDSS.fits %_VHS.fits %_ALLWISE_sum.fits %_GALEX_UL.fits %_HSC.fits
	# this is the same as above, except
	# 1) requires both LS & HSC detections -- could remove LS but then we don't have WISE?
	# 2) uses extended or not from HSC
	stilts tmatchn nin=8 out=$@ \
		in1=$*.fits suffix1= values1=id \
		in2=$*_GALEX.fits suffix2=_GALEX values2=id \
		in3=$*_LS.fits suffix3=_LS values3=id \
		in4=$*_UKIDSS.fits suffix4=_UKIDSS values4=id \
		in5=$*_VHS.fits suffix5=_VHS values5=id \
		in6=$*_ALLWISE_sum.fits suffix6=_ALLWISE values6=id \
		in7=$*_GALEX_UL.fits suffix7=_GALEXUL values7=id \
		in8=$*_HSC.fits suffix8=_HSC values8=id \
		fixcols=all matcher=exact \
		ocmd='addcol pointlike "!(extended_HSC)&&'${ALLOW_POINTLIKE}'"' \
		ocmd='addcol inMzLSBASS "DEC>32&&RA>90&&RA<300"' \
		ocmd='addcol goodfitsLS "(fitbits_LS & (1 | 4 | 8)) == 0 && (maskbits_LS & (1024 | 2048)) == 0"' \
		ocmd='addcol goodfitsg "goodfitsLS && (maskbits_LS & 4) == 0 && fracin_g_LS>0.5"' \
		ocmd='addcol goodfitsr "goodfitsLS && (maskbits_LS & 8) == 0 && fracin_r_LS>0.5"' \
		ocmd='addcol goodfitsz "goodfitsLS && (maskbits_LS & 16) == 0 && fracin_z_LS>0.5"' \
		ocmd='addcol goodfitsi "goodfitsLS && (maskbits_LS & 16384) == 0 && fracin_i_LS>0.5"' \
		ocmd='addcol isolatedLS "(max(fracflux_g_LS,fracflux_r_LS,fracflux_z_LS,fracflux_w1_LS,fracflux_w2_LS)<0.1&&max(fracflux_w3_LS,fracflux_w4_LS)<10)"' \
		ocmd='addcol W34_blended "max(fracflux_w1_LS,fracflux_w2_LS)>0.1||max(fracflux_w4_LS,fracflux_w3_LS)>1"' \
		ocmd='addcol FUV Fflux_real_LU_GALEX*1e26' \
		ocmd='addcol NUV "e_Nflux_real_LU_GALEX>0 ? Nflux_real_LU_GALEX*1e26 : 3 * NUV_fluxlim_GALEXUL"' \
		ocmd='addcol FUV_err "e_Fflux_real_LU_GALEX*1e26"' \
		ocmd='addcol NUV_err "e_Nflux_real_LU_GALEX>0 ? e_Nflux_real_LU_GALEX*1e26 : -NUV"' \
		ocmd='addcol decam_g "(!inMzLSBASS && goodfitsg?LU_flux_g_LS*1e26:-99)"' \
		ocmd='addcol decam_r "(!inMzLSBASS && goodfitsr?LU_flux_r_LS*1e26:-99)"' \
		ocmd='addcol decam_i "(goodfitsi?LU_flux_i_LS*1e26:-99)"' \
		ocmd='addcol decam_z "(!inMzLSBASS && goodfitsz?LU_flux_z_LS*1e26:-99)"' \
		ocmd='addcol decam_g_err "(!inMzLSBASS && goodfitsg?LU_flux_g_err_LS*1e26/fracin_g_LS:-99)"' \
		ocmd='addcol decam_r_err "(!inMzLSBASS && goodfitsr?LU_flux_r_err_LS*1e26/fracin_r_LS:-99)"' \
		ocmd='addcol decam_i_err "(goodfitsi?LU_flux_i_err_LS*1e26/fracin_i_LS:-99)"' \
		ocmd='addcol decam_z_err "(!inMzLSBASS && goodfitsz?LU_flux_z_err_LS*1e26/fracin_z_LS:-99)"' \
		ocmd='addcol 90prime_r "(inMzLSBASS && goodfitsr?LU_flux_r_LS*1e26:-99)"' \
		ocmd='addcol 90prime_g "(inMzLSBASS && goodfitsg?LU_flux_g_LS*1e26:-99)"' \
		ocmd='addcol zd_mosaic "(inMzLSBASS && goodfitsz?LU_flux_z_LS*1e26:-99)"' \
		ocmd='addcol 90prime_r_err "(inMzLSBASS && goodfitsr?LU_flux_r_err_LS*1e26/fracin_r_LS:-99)"' \
		ocmd='addcol zd_mosaic_err "(inMzLSBASS && goodfitsz?LU_flux_z_err_LS*1e26/fracin_z_LS:-99)"' \
		ocmd='addcol 90prime_g_err "(inMzLSBASS && goodfitsg?LU_flux_g_err_LS*1e26/fracin_g_LS:-99)"' \
		icmd8='addcol psf_g_usable "!g_pixelflags_saturated&&!g_pixelflags_saturatedcenter&&g_inputcount_value>0&&!g_inputcount_flag&&!g_localbackground_flag&&!g_pixelflags&&!g_psfflux_flag"' \
		icmd8='addcol psf_r_usable "!r_pixelflags_saturated&&!r_pixelflags_saturatedcenter&&r_inputcount_value>0&&!r_inputcount_flag&&!r_localbackground_flag&&!r_pixelflags&&!r_psfflux_flag"' \
		icmd8='addcol psf_i_usable "!i_pixelflags_saturated&&!i_pixelflags_saturatedcenter&&i_inputcount_value>0&&!i_inputcount_flag&&!i_localbackground_flag&&!i_pixelflags&&!i_psfflux_flag"' \
		icmd8='addcol psf_z_usable "!z_pixelflags_saturated&&!z_pixelflags_saturatedcenter&&z_inputcount_value>0&&!z_inputcount_flag&&!z_localbackground_flag&&!z_pixelflags&&!z_psfflux_flag"' \
		icmd8='addcol psf_y_usable "!y_pixelflags_saturated&&!y_pixelflags_saturatedcenter&&y_inputcount_value>0&&!y_inputcount_flag&&!y_localbackground_flag&&!y_pixelflags&&!y_psfflux_flag"' \
		ocmd='addcol HSC_g "(psf_g_usable_HSC&&!g_psfflux_flag_HSC&&!g_apertureflux_57_flag_HSC?LU_flux_g_HSC*1e26:-99)"' \
		ocmd='addcol HSC_r "(psf_r_usable_HSC&&!r_psfflux_flag_HSC&&!r_apertureflux_57_flag_HSC?LU_flux_r_HSC*1e26:-99)"' \
		ocmd='addcol HSC_i "(psf_i_usable_HSC&&!i_psfflux_flag_HSC&&!i_apertureflux_57_flag_HSC?LU_flux_i_HSC*1e26:-99)"' \
		ocmd='addcol HSC_z "(psf_z_usable_HSC&&!z_psfflux_flag_HSC&&!z_apertureflux_57_flag_HSC?LU_flux_z_HSC*1e26:-99)"' \
		ocmd='addcol HSC_y "(psf_y_usable_HSC&&!y_psfflux_flag_HSC&&!y_apertureflux_57_flag_HSC?LU_flux_y_HSC*1e26:-99)"' \
		ocmd='addcol HSC_g_err "(psf_g_usable_HSC&&!g_psfflux_flag_HSC&&!g_apertureflux_57_flag_HSC?LU_flux_g_err_HSC*1e26:-99)"' \
		ocmd='addcol HSC_r_err "(psf_r_usable_HSC&&!r_psfflux_flag_HSC&&!r_apertureflux_57_flag_HSC?LU_flux_r_err_HSC*1e26:-99)"' \
		ocmd='addcol HSC_i_err "(psf_i_usable_HSC&&!i_psfflux_flag_HSC&&!i_apertureflux_57_flag_HSC?LU_flux_i_err_HSC*1e26:-99)"' \
		ocmd='addcol HSC_z_err "(psf_z_usable_HSC&&!z_psfflux_flag_HSC&&!z_apertureflux_57_flag_HSC?LU_flux_z_err_HSC*1e26:-99)"' \
		ocmd='addcol HSC_y_err "(psf_y_usable_HSC&&!y_psfflux_flag_HSC&&!y_apertureflux_57_flag_HSC?LU_flux_y_err_HSC*1e26:-99)"' \
		ocmd='addcol UV_Y "Yerrbits_VHS==0?(pointlike?UV_Yapc4flux_VHS:UV_Yap6flux_VHS)*1e26:-99"' \
		ocmd='addcol UV_J "Jerrbits_VHS==0?(pointlike?UV_Japc4flux_VHS:UV_Jap6flux_VHS)*1e26:-99"' \
		ocmd='addcol UV_H "Herrbits_VHS==0?(pointlike?UV_Hapc4flux_VHS:UV_Hap6flux_VHS)*1e26:-99"' \
		ocmd='addcol UV_K "Kserrbits_VHS==0?(pointlike?UV_Ksapc4flux_VHS:UV_Ksap6flux_VHS)*1e26:-99"' \
		ocmd='addcol UV_Y_err "Yerrbits_VHS==0?(pointlike?UV_Yapc4flux_err_VHS:UV_Yap6flux_err_VHS)*1e26:-99"' \
		ocmd='addcol UV_J_err "Jerrbits_VHS==0?(pointlike?UV_Japc4flux_err_VHS:UV_Jap6flux_err_VHS)*1e26:-99"' \
		ocmd='addcol UV_H_err "Herrbits_VHS==0?(pointlike?UV_Hapc4flux_err_VHS:UV_Hap6flux_err_VHS)*1e26:-99"' \
		ocmd='addcol UV_K_err "Kserrbits_VHS==0?(pointlike?UV_Ksapc4flux_err_VHS:UV_Ksap6flux_err_VHS)*1e26:-99"' \
		ocmd='addcol WFCAM_Y "Yerrbits_UKIDSS==0?(pointlike?WFCAM_Yapc4flux_UKIDSS:WFCAM_Yap6flux_UKIDSS)*1e26:-99"' \
		ocmd='addcol WFCAM_J "Jerrbits_UKIDSS==0?(pointlike?WFCAM_Japc4flux_UKIDSS:WFCAM_Jap6flux_UKIDSS)*1e26:-99"' \
		ocmd='addcol WFCAM_H "Herrbits_UKIDSS==0?(pointlike?WFCAM_Hapc4flux_UKIDSS:WFCAM_Hap6flux_UKIDSS)*1e26:-99"' \
		ocmd='addcol WFCAM_Ks "Kerrbits_UKIDSS==0?(pointlike?WFCAM_Kapc4flux_UKIDSS:WFCAM_Kap6flux_UKIDSS)*1e26:-99"' \
		ocmd='addcol WFCAM_Y_err "Yerrbits_UKIDSS==0?(pointlike?WFCAM_Yapc4flux_err_UKIDSS:WFCAM_Yap6flux_err_UKIDSS)*1e26:-99"' \
		ocmd='addcol WFCAM_J_err "Jerrbits_UKIDSS==0?(pointlike?WFCAM_Japc4flux_err_UKIDSS:WFCAM_Jap6flux_err_UKIDSS)*1e26:-99"' \
		ocmd='addcol WFCAM_H_err "Herrbits_UKIDSS==0?(pointlike?WFCAM_Hapc4flux_err_UKIDSS:WFCAM_Hap6flux_err_UKIDSS)*1e26:-99"' \
		ocmd='addcol WFCAM_Ks_err "Kerrbits_UKIDSS==0?(pointlike?WFCAM_Kapc4flux_err_UKIDSS:WFCAM_Kap6flux_err_UKIDSS)*1e26:-99"' \
		ocmd='addcol WISE1_origin "isolatedLS&&goodfitsW1 ? \"LS10\" : ( (fracflux_w1_LS > 1 || WISE1_ALLWISE < LU_flux_w1_LS * 1e26 * (1 + fracflux_w1_LS)) ? \"AllWISE\" : \"LS10UL\")"' \
		ocmd='addcol WISE1        "isolatedLS&&goodfitsW1 ? LU_flux_w1_LS * 1e26 : ((fracflux_w1_LS > 1 || WISE1_ALLWISE < LU_flux_w1_LS * 1e26 * (1 + fracflux_w1_LS)) ? WISE1_ALLWISE : (LU_flux_w1_LS + LU_flux_w1_err_LS) * 1e26 * (1 + fracflux_w1_LS))"' \
		ocmd='addcol WISE1_err "isolatedLS&&goodfitsW1 ? LU_flux_w1_err_LS * 1e26 : ((fracflux_w1_LS > 1 || WISE1_ALLWISE < LU_flux_w1_LS * 1e26 * (1 + fracflux_w1_LS)) ? -WISE1_ALLWISE : -LU_flux_w1_err_LS * 1e26 * (1 + fracflux_w1_LS))"' \
		ocmd='addcol WISE2_origin "isolatedLS&&goodfitsW2 ? \"LS10\" : ((fracflux_w2_LS > 1 || WISE2_ALLWISE < LU_flux_w2_LS * 1e26 * (1 + fracflux_w2_LS)) ? \"AllWISE\" : \"LS10UL\")"' \
		ocmd='addcol WISE2 "isolatedLS&&goodfitsW2 ? LU_flux_w2_LS * 1e26 : ((fracflux_w2_LS > 1 || WISE2_ALLWISE < LU_flux_w2_LS * 1e26 * (1 + fracflux_w2_LS)) ? WISE2_ALLWISE : (LU_flux_w2_LS + LU_flux_w2_err_LS) * 1e26 * (1 + fracflux_w2_LS))"' \
		ocmd='addcol WISE2_err "isolatedLS&&goodfitsW2 ? LU_flux_w2_err_LS * 1e26 : ((fracflux_w2_LS > 1 || WISE2_ALLWISE < LU_flux_w2_LS * 1e26 * (1 + fracflux_w2_LS)) ? -WISE2_ALLWISE : -LU_flux_w2_err_LS * 1e26 * (1 + fracflux_w2_LS))"' \
		ocmd='addcol WISE34_origin "isolatedLS&&goodfitsW1&&goodfitsW2&&!W34_blended ? \"LS10\" : \"AllWISE\""' \
		ocmd='addcol WISE3 "isolatedLS&&goodfitsW1&&goodfitsW2&&!W34_blended ? LU_flux_w3_LS*1e26 : -99"' \
		ocmd='addcol WISE3_err "isolatedLS&&goodfitsW1&&goodfitsW2&&!W34_blended ? LU_flux_w3_err_LS*1e26 : -99"' \
		ocmd='addcol WISE4 "isolatedLS&&goodfitsW1&&goodfitsW2&&!W34_blended ? LU_flux_w4_LS*1e26 : -99"' \
		ocmd='addcol WISE4_err "isolatedLS&&goodfitsW1&&goodfitsW2&&!W34_blended ? LU_flux_w4_err_LS*1e26 : -99"' \
