%_VHS.fits: %_VHS_noirlab.fits
	# for point sources: aperture corrected fluxes within 2.8'' diameter aperture (apermag4)
	# for extended sources: not-aperture corrected fluxes within 5.6'' diameter aperture (apermagnoapercorr6)
	# filter-specific offsets is the conversion from Vega to AB
	# 48.6 is to convert to erg/s/cm^2/Hz
	stilts tpipe in=$^ out=$@ \
		cmd='addcol UV_Yapc4flux "Yapermag4>0 ? (pow(10, -(Yapermag4-ay+0.634+48.60)/(2.5))): -99."' \
		cmd='addcol UV_Yapc4flux_err "Yapermag4err>0 ? (pow(10, -(Yapermag4-ay-Yapermag4err+0.634+48.60)/(2.5)))-UV_Yapc4flux: -99."' \
		cmd='addcol UV_Japc4flux "Japermag4>0 ? (pow(10, -(Japermag4-aj+0.938+48.60)/(2.5))): -99."' \
		cmd='addcol UV_Japc4flux_err "Japermag4err>0 ? (pow(10, -(Japermag4-aj-Japermag4err+0.938+48.60)/(2.5)))-UV_Japc4flux: -99."' \
		cmd='addcol UV_Hapc4flux "Hapermag4>0 ? (pow(10, -(Hapermag4-ah+1.379+48.60)/(2.5))): -99."' \
		cmd='addcol UV_Hapc4flux_err "Hapermag4err>0 ? (pow(10, -(Hapermag4-ah-Hapermag4err+1.379+48.60)/(2.5)))-UV_Hapc4flux: -99."' \
		cmd='addcol UV_Ksapc4flux "Ksapermag4>0 ? (pow(10, -(Ksapermag4-aks+1.9+48.60)/(2.5))): -99."' \
		cmd='addcol UV_Ksapc4flux_err "Ksapermag4err>0 ? (pow(10, -(Ksapermag4-aks-Ksapermag4err+1.9+48.60)/(2.5)))-UV_Ksapc4flux: -99."' \
		\
		cmd='addcol UV_Yap6flux "Yapermagnoapercorr6>0 ? (pow(10, -(Yapermagnoapercorr6-ay+0.634+48.60)/(2.5))): -99."' \
		cmd='addcol UV_Yap6flux_err "Yapermag6err>0 ? (pow(10, -(Yapermagnoapercorr6-ay-Yapermag6err+0.634+48.60)/(2.5)))-UV_Yap6flux: -99."' \
		cmd='addcol UV_Jap6flux "Japermagnoapercorr6>0 ? (pow(10, -(Japermagnoapercorr6-aj+0.938+48.60)/(2.5))): -99."' \
		cmd='addcol UV_Jap6flux_err "Japermag6err>0 ? (pow(10, -(Japermagnoapercorr6-aj-Japermag6err+0.938+48.60)/(2.5)))-UV_Jap6flux: -99."' \
		cmd='addcol UV_Hap6flux "Hapermagnoapercorr6>0 ? (pow(10, -(Hapermagnoapercorr6-ah+1.379+48.60)/(2.5))): -99."' \
		cmd='addcol UV_Hap6flux_err "Hapermag6err>0 ? (pow(10, -(Hapermagnoapercorr6-ah-Hapermag6err+1.379+48.60)/(2.5)))-UV_Hap6flux: -99."' \
		cmd='addcol UV_Ksap6flux "Ksapermagnoapercorr6>0 ? (pow(10, -(Ksapermagnoapercorr6-aks+1.9+48.60)/(2.5))): -99."' \
		cmd='addcol UV_Ksap6flux_err "Ksapermag6err>0 ? (pow(10, -(Ksapermagnoapercorr6-aks-Ksapermag6err+1.9+48.60)/(2.5)))-UV_Ksap6flux: -99."' \
		cmd='addcol UV_Yapc6flux "Yapermagnoapercorr6>0 ? (pow(10, -(Yapermagnoapercorr6-ay+0.634+48.60)/(2.5))): -99."' \
		cmd='addcol UV_Yapc6flux_err "Yapermag6err>0 ? (pow(10, -(Yapermagnoapercorr6-ay-Yapermag6err+0.634+48.60)/(2.5)))-UV_Yapc6flux: -99."'
