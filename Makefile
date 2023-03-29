.PHONY: all
.SUFFIXES: # no built-in rules
.SECONDARY: # do not delete intermediate products

all:
	echo "run: make myfilename_clean.fits myfilename_all.fits"

# p: Y, J, H, Ks e_
# ap3: Y, J, H, Ks e_
# apc3: Y, J, H, Ks e_
# flux_u g r i z eflux
# flux_y j h k w1 w2
%_coords.csv: %.fits
	stilts tpipe in=$^ cmd='keepcols "ID RA DEC"' out=$@

%_LS_aper.fits: fetchLS.py %.fits
	python3 $^ $@

%_LS.fits: %_LS_aper.fits
	# adflux_* = adaptive flux column, depending on source type
	# for extended sources: 5'' aperture = OPT:apflux_*_7 and IR:apflux_*_2
	# for point sources: model flux
	# Milky way extinction correction
	# 28.44 is to convert to erg/s/cm^2/Hz
	stilts tpipe in=$^ out=$@ \
		cmd='addcol adflux_g "type==\"PSF\" ? flux_g : apflux_g_7"' \
		cmd='addcol adflux_r "type==\"PSF\" ? flux_r : apflux_r_7"' \
		cmd='addcol adflux_z "type==\"PSF\" ? flux_z : apflux_z_7"' \
		cmd='addcol adflux_W1 "type==\"PSF\" ? flux_W1 : apflux_W1_2"' \
		cmd='addcol adflux_W2 "type==\"PSF\" ? flux_W2 : apflux_W2_2"' \
		cmd='addcol adflux_W3 "type==\"PSF\" ? flux_W3 : apflux_W3_2"' \
		cmd='addcol adflux_W4 "type==\"PSF\" ? flux_W4 : apflux_W4_2"' \
		cmd='addcol adflux_ivar_g  "type==\"PSF\" ? flux_ivar_g : apflux_ivar_g_7"' \
		cmd='addcol adflux_ivar_r  "type==\"PSF\" ? flux_ivar_r : apflux_ivar_r_7"' \
		cmd='addcol adflux_ivar_z  "type==\"PSF\" ? flux_ivar_z : apflux_ivar_z_7"' \
		cmd='addcol adflux_ivar_W1 "type==\"PSF\" ? flux_ivar_W1 : apflux_ivar_W1_2"' \
		cmd='addcol adflux_ivar_W2 "type==\"PSF\" ? flux_ivar_W2 : apflux_ivar_W2_2"' \
		cmd='addcol adflux_ivar_W3 "type==\"PSF\" ? flux_ivar_W3 : apflux_ivar_W3_2"' \
		cmd='addcol adflux_ivar_W4 "type==\"PSF\" ? flux_ivar_W4 : apflux_ivar_W4_2"' \
		cmd='addcol LU_flux_g "FLUX_G>-10000000000L ? ((ADFLUX_G/MW_TRANSMISSION_G))*pow(10, -28.44):-99."' \
		cmd='addcol LU_flux_r "FLUX_R>-10000000000L ? ((ADFLUX_R/MW_TRANSMISSION_R))*pow(10, -28.44):-99."' \
		cmd='addcol LU_flux_z "FLUX_Z>-10000000000L ? ((ADFLUX_Z/MW_TRANSMISSION_Z))*pow(10, -28.44):-99."' \
		cmd='addcol LU_flux_w1 "FLUX_W1>-10000000000L ? ((ADFLUX_W1/MW_TRANSMISSION_W1))*pow(10, -28.44):-99."' \
		cmd='addcol LU_flux_w2 "FLUX_W2>-10000000000L ? ((ADFLUX_W2/MW_TRANSMISSION_W2))*pow(10, -28.44):-99."' \
		cmd='addcol LU_flux_w3 "FLUX_W3>-10000000000L ? ((ADFLUX_W3/MW_TRANSMISSION_W3))*pow(10, -28.44):-99."' \
		cmd='addcol LU_flux_w4 "FLUX_W4>-10000000000L ? ((ADFLUX_W4/MW_TRANSMISSION_W4))*pow(10, -28.44):-99."' \
		cmd='addcol LU_flux_g_err "ADFLUX_IVAR_G>0 ? (((1/pow(ADFLUX_IVAR_G, 0.5))/MW_TRANSMISSION_G))*pow(10, -28.44):-99."' \
		cmd='addcol LU_flux_r_err "ADFLUX_IVAR_R>0 ? (((1/pow(ADFLUX_IVAR_R, 0.5))/MW_TRANSMISSION_R))*pow(10, -28.44):-99."' \
		cmd='addcol LU_flux_z_err "ADFLUX_IVAR_Z>0 ? (((1/pow(ADFLUX_IVAR_Z, 0.5))/MW_TRANSMISSION_Z))*pow(10, -28.44):-99."' \
		cmd='addcol LU_flux_w1_err "ADFLUX_IVAR_W1>0 ? (((1/pow(ADFLUX_IVAR_W1, 0.5))/MW_TRANSMISSION_W1))*pow(10, -28.44):-99."' \
		cmd='addcol LU_flux_w2_err "ADFLUX_IVAR_W2>0 ? (((1/pow(ADFLUX_IVAR_W2, 0.5))/MW_TRANSMISSION_W2))*pow(10, -28.44):-99."' \
		cmd='addcol LU_flux_w3_err "ADFLUX_IVAR_W3>0 ? (((1/pow(ADFLUX_IVAR_W3, 0.5))/MW_TRANSMISSION_W3))*pow(10, -28.44):-99."' \
		cmd='addcol LU_flux_w4_err "ADFLUX_IVAR_W4>0 ? (((1/pow(ADFLUX_IVAR_W4, 0.5))/MW_TRANSMISSION_W4))*pow(10, -28.44):-99."'

%_GALEX.fits: %_coords.csv
	# NUV measurements are only valid down to NUVmag<20.8 AB mag NFlux>17.4
	# FUV measurements are only valid down to FUVmag<19.9 AB mag FFlux>39.7
	stilts cdsskymatch cdstable=II/335/galex_ais  \
		find=each in=$< ra=RA dec=DEC \
		radius=1 out=$@ \
		ocmd='addcol EBmV "E(B-V)"' \
		ocmd='addcol FUV "FUVmag>0&&FUVmag<20.8 ? FUVmag-(8.286*$22): -99."' \
		ocmd='addcol NUV "NUVmag>0&&NUVmag<19.9 ? NUVmag-(8.621*$22): -99."' \
		ocmd='addcol e_FUV "FUVmag>0&&FUVmag<20.8 ?e_FUVmag:-99"' \
		ocmd='addcol e_NUV "NUVmag>0&&NUVmag<19.9?e_NUVmag:-99"' \
		ocmd='addcol Fflux_LU_A "Fflux>39.7 ? Fflux*(1.26727*pow(10, -17)): -99."' \
		ocmd='addcol Nflux_LU_A "Nflux>17.4 ? Nflux*(5.59444*pow(10, -18)): -99."' \
		ocmd='addcol e_Fflux_LU_A "e_Fflux>0 ? e_Fflux*(1.26727*pow(10, -17)): -99."' \
		ocmd='addcol e_Nflux_LU_A "e_Nflux>0 ? e_Nflux*(5.59444*pow(10, -18)): -99."' \
		ocmd='addcol Fflux_LU "Fflux>39.7 ? Fflux*(pow(10, -29)): -99."' \
		ocmd='addcol Nflux_LU "Nflux>17.4 ? Nflux*(pow(10, -29)): -99."' \
		ocmd='addcol e_Fflux_LU "e_Fflux>0 ? e_Fflux*(pow(10, -29)): -99."' \
		ocmd='addcol e_Nflux_LU "e_Nflux>0 ? e_Nflux*(pow(10, -29)): -99."'  \
		ocmd='addcol A_v "EBmV>0 ? (3.1*EBmV): -99."' \
		ocmd='addcol A_lam_fuv "A_v>0 ? A_v*(-0.3459882441104826 + 9.173790513029425/3.1): -99."' \
		ocmd='addcol A_lam_nuv "A_v>0 ? A_v*(0.6743929174392934 + 1.2480067831962458/3.1): -99."' \
		ocmd='addcol Fflux_real_LU_A "Fflux_LU_A>0 ? Fflux_LU_A*(pow(10, (0.4*A_lam_fuv))): -99."' \
		ocmd='addcol Nflux_real_LU_A "Nflux_LU_A>0 ? Nflux_LU_A*(pow(10, (0.4*A_lam_nuv))): -99."' \
		ocmd='addcol e_Fflux_real_LU_A "e_Fflux_LU_A>0 ? e_Fflux_LU_A*(pow(10, (0.4*A_lam_fuv))): -99."' \
		ocmd='addcol e_Nflux_real_LU_A "e_Nflux_LU_A>0 ? e_Nflux_LU_A*(pow(10, (0.4*A_lam_nuv))): -99."' \
		ocmd='addcol Fflux_real_LU "Fflux_real_LU_A>0 ? Fflux_real_LU_A*(3.34*pow(10, -19))*(pow(1538.6, 2)): -99."' \
		ocmd='addcol Nflux_real_LU "Nflux_real_LU_A>0 ? Nflux_real_LU_A*(3.34*pow(10, -19))*(pow(2315.7, 2)): -99."' \
		ocmd='addcol e_Fflux_real_LU "e_Fflux_real_LU_A>0 ? e_Fflux_real_LU_A*(3.34*pow(10, -19))*(pow(1538.6, 2)): -99."' \
		ocmd='addcol e_Nflux_real_LU "e_Nflux_real_LU_A>0 ? e_Nflux_real_LU_A*(3.34*pow(10, -19))*(pow(2315.7, 2)): -99."' \


%_ALLWISE.fits: %_coords.csv
	# match to AllWISE
	# convert Vega magnitudes to mJy fluxes
	stilts cdsskymatch cdstable=II/328/allwise  \
		find=each in=$< ra=RA dec=DEC \
		radius=4 out=$@ \
		ocmd='addcol WISE4 "W4mag>0 ? (pow(10, -(W4mag+6.620+48.60)/(2.5)))*1e26: -99."' \
		ocmd='addcol WISE4_err "e_W4mag>0 ? (pow(10, -(W4mag-e_W4mag+6.620+48.60)/(2.5)))*1e26-WISE4: -99."' \
		ocmd='addcol WISE3 "W3mag>0 ? (pow(10, -(W3mag+5.174+48.60)/(2.5)))*1e26: -99."' \
		ocmd='addcol WISE3_err "e_W3mag>0 ? (pow(10, -(W3mag-e_W3mag+5.174+48.60)/(2.5)))*1e26-WISE3: -99."' \
		ocmd='addcol WISE2 "W2mag>0 ? (pow(10, -(W2mag+3.339+48.60)/(2.5)))*1e26: -99."' \
		ocmd='addcol WISE2_err "e_W2mag>0 ? (pow(10, -(W2mag-e_W2mag+3.339+48.60)/(2.5)))*1e26-WISE2: -99."' \
		ocmd='addcol WISE1 "W1mag>0 ? (pow(10, -(W1mag+2.699+48.60)/(2.5)))*1e26: -99."' \
		ocmd='addcol WISE1_err "e_W1mag>0 ? (pow(10, -(W1mag-e_W1mag+2.699+48.60)/(2.5)))*1e26-WISE1: -99."' 


%_VHS.csv: %_coords.csv
	@echo
	@echo "Steps to manually create $@:"
	@echo
	@echo " 1. go to https://datalab.noirlab.edu/xmatch.php"
	@echo " 2. In the Table Management tab, upload the file $<. Give it a name."
	@echo " 3. In the Xmatch tab, select the catalog file you just uploaded, and the ra dec columns"
	@echo " 4. pick as 2nd table "
	@echo "      - database 'vhs_dr5'. Wait. "
	@echo "      - select as table 'vhs_dr5.vhs_cat_v3'. Wait. "
	@echo "      - select ra2000 and dec2000 as coordinate columns."
	@echo "      - select 'All' for which columns to select."
	@echo " 5. For Result Table Name set 'vhs'"
	@echo " 6. Select radius 1 arcsec, 'find all matches'"
	@echo " 7. tick 'Download results to your computer only'"
	@echo " 8. click submit, which will download a csv file."
	@echo " 9. save file to ${PWD}/$@"
	exit 1

%_UKIDSS.csv: %_coords.csv
	@echo
	@echo "Steps to manually create $@:"
	@echo
	@echo " 1. go to https://datalab.noirlab.edu/xmatch.php"
	@echo " 2. In the Table Management tab, upload the file $<. Give it a name."
	@echo " 3. In the Xmatch tab, select the catalog file you just uploaded, and the ra dec columns"
	@echo " 4. pick as 2nd table "
	@echo "      - database 'ukidss_dr11plus'. Wait. "
	@echo "      - select as table 'ukidss_dr11plus.lassource'. Wait. "
	@echo " 5. For Result Table Name set 'ukidss'"
	@echo " 6. Select radius 1 arcsec, 'find all matches'"
	@echo " 7. tick 'Download results to your computer only'"
	@echo " 8. click submit, which will download a csv file."
	@echo " 9. save file to ${PWD}/$@"
	@exit 1

%_VHS.fits: %_VHS.csv
	# for point sources: aperture corrected fluxes within 2.8'' (apermag4)
	# for extended sources: not-aperture corrected fluxes within 5.6'' (apermagnoapercorr6)
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

%_UKIDSS.fits: %_UKIDSS.csv
	# for point sources: fluxes within 2.8'' (apermag4)
	# for extended sources: fluxes within 5.6'' (apercorr6)
	#    !!! there are ONLY aperture corrected fluxes provided!
	# filter-specific offsets is the conversion from Vega to AB
	# 48.6 is to convert to erg/s/cm^2/Hz
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


%_all.fits: %.fits %_GALEX.fits %_LS.fits %_UKIDSS.fits %_VHS.fits %_ALLWISE.fits
	# merge everything together and use sensible column names
	# keep only WISE fluxes when ALLWISE also has a detection there
	# and if there are no blending issues
	# check LS9 fitbits for issues. bits 1, 2, 3, 7, 9, 12 are acceptable.
	# otherwise delete
	stilts tmatchn nin=6 out=$@ \
		in1=$*.fits suffix1= values1=id \
		in2=$*_GALEX.fits suffix2=_GALEX values2=id \
		in3=$*_LS.fits suffix3=_LS values3=id \
		in4=$*_UKIDSS.fits suffix4=_UKIDSS values4=t1_id \
		in5=$*_VHS.fits suffix5=_VHS values5=t1_id \
		in6=$*_ALLWISE.fits suffix6=_ALLWISE values6=id_in \
		fixcols=all matcher=exact \
		ocmd='addcol pointlike "!(type_LS != \"PSF\")"' \
		ocmd='addcol FUV Fflux_real_LU_GALEX*1e26' \
		ocmd='addcol NUV Nflux_real_LU_GALEX*1e26' \
		ocmd='addcol FUV_err e_Fflux_real_LU_GALEX*1e26' \
		ocmd='addcol NUV_err e_Nflux_real_LU_GALEX*1e26' \
		ocmd='addcol decam_g LU_flux_g_LS*1e26' \
		ocmd='addcol decam_r LU_flux_r_LS*1e26' \
		ocmd='addcol decam_z LU_flux_z_LS*1e26' \
		ocmd='addcol decam_g_err LU_flux_g_err_LS*1e26' \
		ocmd='addcol decam_r_err LU_flux_r_err_LS*1e26' \
		ocmd='addcol decam_z_err LU_flux_z_err_LS*1e26' \
		ocmd='addcol UV_Y "(pointlike?UV_Yapc4flux_VHS:UV_Yap6flux_VHS)*1e26"' \
		ocmd='addcol UV_J "(pointlike?UV_Japc4flux_VHS:UV_Jap6flux_VHS)*1e26"' \
		ocmd='addcol UV_H "(pointlike?UV_Hapc4flux_VHS:UV_Hap6flux_VHS)*1e26"' \
		ocmd='addcol UV_K "(pointlike?UV_Ksapc4flux_VHS:UV_Ksap6flux_VHS)*1e26"' \
		ocmd='addcol UV_Y_err "(pointlike?UV_Yapc4flux_err_VHS:UV_Yap6flux_err_VHS)*1e26"' \
		ocmd='addcol UV_J_err "(pointlike?UV_Japc4flux_err_VHS:UV_Jap6flux_err_VHS)*1e26"' \
		ocmd='addcol UV_H_err "(pointlike?UV_Hapc4flux_err_VHS:UV_Hap6flux_err_VHS)*1e26"' \
		ocmd='addcol UV_K_err "(pointlike?UV_Ksapc4flux_err_VHS:UV_Ksap6flux_err_VHS)*1e26"' \
		ocmd='addcol WFCAM_Y "(pointlike?WFCAM_Yapc4flux_UKIDSS:WFCAM_Yap6flux_UKIDSS)*1e26"' \
		ocmd='addcol WFCAM_J "(pointlike?WFCAM_Japc4flux_UKIDSS:WFCAM_Jap6flux_UKIDSS)*1e26"' \
		ocmd='addcol WFCAM_H "(pointlike?WFCAM_Hapc4flux_UKIDSS:WFCAM_Hap6flux_UKIDSS)*1e26"' \
		ocmd='addcol WFCAM_Ks "(pointlike?WFCAM_Kapc4flux_UKIDSS:WFCAM_Kap6flux_UKIDSS)*1e26"' \
		ocmd='addcol WFCAM_Y_err "(pointlike?WFCAM_Yapc4flux_err_UKIDSS:WFCAM_Yap6flux_err_UKIDSS)*1e26"' \
		ocmd='addcol WFCAM_J_err "(pointlike?WFCAM_Japc4flux_err_UKIDSS:WFCAM_Jap6flux_err_UKIDSS)*1e26"' \
		ocmd='addcol WFCAM_H_err "(pointlike?WFCAM_Hapc4flux_err_UKIDSS:WFCAM_Hap6flux_err_UKIDSS)*1e26"' \
		ocmd='addcol WFCAM_Ks_err "(pointlike?WFCAM_Kapc4flux_err_UKIDSS:WFCAM_Kap6flux_err_UKIDSS)*1e26"' \
		ocmd='addcol goodfits_LS "((fitbits_LS & (16383 ^ (2 | 4 | 8 | 128 | 512 | 4096))) == 0)"' \
		ocmd='addcol isolated_LS "(max(fracflux_g_LS,fracflux_z_LS,fracflux_w2_LS)<0.1&&max(fracflux_w4_LS,fracflux_w3_LS)<10)"' \
		ocmd='addcol W34_blended "max(fracflux_w1_LS,fracflux_w2_LS)>0.1||max(fracflux_w4_LS,fracflux_w3_LS)>1"' \
		ocmd='addcol WISE1 "(WISE1_ERR_ALLWISE>0?LU_flux_w1_LS*1e26:-99)"' \
		ocmd='addcol WISE1_err "(WISE1_ERR_ALLWISE>0?LU_flux_w1_err_LS*1e26:-99)"' \
		ocmd='addcol WISE2 "(WISE2_ERR_ALLWISE>0?LU_flux_w2_LS*1e26:-99)"' \
		ocmd='addcol WISE2_err "(WISE2_ERR_ALLWISE>0?LU_flux_w2_err_LS*1e26:-99)"' \
		ocmd='addcol WISE3 "(WISE3_ERR_ALLWISE>0&&!W34_blended?LU_flux_w3_LS*1e26:-99)"' \
		ocmd='addcol WISE3_err "(WISE3_ERR_ALLWISE>0&&!W34_blended?LU_flux_w3_err_LS*1e26:-99)"' \
		ocmd='addcol WISE4 "(WISE4_ERR_ALLWISE>0&&!W34_blended?LU_flux_w4_LS*1e26:-99)"' \
		ocmd='addcol WISE4_err "(WISE4_ERR_ALLWISE>0&&!W34_blended?LU_flux_w4_err_LS*1e26:-99)"' \

%_lite.fits: %_all.fits
	stilts tpipe in=$^ out=$@ cmd='delcols "adflux*_LS LU_flux*_LS *_GALEX *_VHS *_UKIDSS _*ALLWISE"'
	stilts tpipe in=$@ omode=stats

%.fits_errors.pdf %.fits_fluxes.pdf: %.fits
	python3 checkphotometry.py $<
