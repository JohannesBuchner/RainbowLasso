%_UKIDSS.fits: %_UKIDSS_noirlab.fits
	# for point sources: aperture fluxes with diameter 2.8'' (apermag4)
	# for extended sources: aperture fluxes with diameter 5.6'' (apercorr6)
	#    !!! there are ONLY aperture corrected fluxes provided!
	# filter-specific offsets is the conversion from Vega to AB
	# 48.6 is to convert to erg/s/cm^2/Hz
	# AB offsets from http://wsa.roe.ac.uk/pubs/hewett-ukidss.pdf
	stilts tpipe in=$^ out=$@ \
		cmd='addcol WFCAM_Yapc4flux "Yapermag4>0 ? (pow(10, -(Yapermag4-ay+0.634+48.60)/(2.5))): -99."' \
		cmd='addcol WFCAM_Yapc4flux_err "Yapermag4err>0 ? (pow(10, -(Yapermag4-ay-Yapermag4err+0.634+48.60)/(2.5)))-WFCAM_Yapc4flux: -99."' \
		cmd='addcol WFCAM_Japc4flux "Japermag4>0 ? (pow(10, -(Japermag4-aj+0.938+48.60)/(2.5))): -99."' \
		cmd='addcol WFCAM_Japc4flux_err "Japermag4err>0 ? (pow(10, -(Japermag4-aj-Japermag4err+0.938+48.60)/(2.5)))-WFCAM_Japc4flux: -99."' \
		cmd='addcol WFCAM_Hapc4flux "Hapermag4>0 ? (pow(10, -(Hapermag4-ah+1.379+48.60)/(2.5))): -99."' \
		cmd='addcol WFCAM_Hapc4flux_err "Hapermag4err>0 ? (pow(10, -(Hapermag4-ah-Hapermag4err+1.379+48.60)/(2.5)))-WFCAM_Hapc4flux: -99."' \
		cmd='addcol WFCAM_Kapc4flux "Kapermag4>0 ? (pow(10, -(Kapermag4-ak+1.9+48.60)/(2.5))): -99."' \
		cmd='addcol WFCAM_Kapc4flux_err "Kapermag4err>0 ? (pow(10, -(Kapermag4-ak-Kapermag4err+1.9+48.60)/(2.5)))-WFCAM_Kapc4flux: -99."' \
		\
		cmd='addcol WFCAM_Yap6flux "Yapermag6>0 ? (pow(10, -(Yapermag6-ay+0.634+48.60)/(2.5))): -99."' \
		cmd='addcol WFCAM_Yap6flux_err "Yapermag6err>0 ? (pow(10, -(Yapermag6-ay-Yapermag6err+0.634+48.60)/(2.5)))-WFCAM_Yap6flux: -99."' \
		cmd='addcol WFCAM_Jap6flux "Japermag6>0 ? (pow(10, -(Japermag6-aj+0.938+48.60)/(2.5))): -99."' \
		cmd='addcol WFCAM_Jap6flux_err "Japermag6err>0 ? (pow(10, -(Japermag6-aj-Japermag6err+0.938+48.60)/(2.5)))-WFCAM_Jap6flux: -99."' \
		cmd='addcol WFCAM_Hap6flux "Hapermag6>0 ? (pow(10, -(Hapermag6-ah+1.379+48.60)/(2.5))): -99."' \
		cmd='addcol WFCAM_Hap6flux_err "Hapermag6err>0 ? (pow(10, -(Hapermag6-ah-Hapermag6err+1.379+48.60)/(2.5)))-WFCAM_Hap6flux: -99."' \
		cmd='addcol WFCAM_Kap6flux "Kapermag6>0 ? (pow(10, -(Kapermag6-ak+1.9+48.60)/(2.5))): -99."' \
		cmd='addcol WFCAM_Kap6flux_err "Kapermag6err>0 ? (pow(10, -(Kapermag6-ak-Kapermag6err+1.9+48.60)/(2.5)))-WFCAM_Kap6flux: -99."' \
		cmd='addcol WFCAM_Yapc6flux "Yapermag6>0 ? (pow(10, -(Yapermag6-ay+0.634+48.60)/(2.5))): -99."' \
		cmd='addcol WFCAM_Yapc6flux_err "Yapermag6err>0 ? (pow(10, -(Yapermag6-ay-Yapermag6err+0.634+48.60)/(2.5)))-WFCAM_Yapc6flux: -99."'
