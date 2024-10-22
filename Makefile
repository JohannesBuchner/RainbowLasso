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

%_LS_staraper.fits: fetchLSstar.py %.fits
	python3 $^ $@

%_LS_extaper_in.fits: %_LS_staraper.fits %_LS_aper.fits
	stilts tmatch2 out=$@ fixcols=all matcher=exact \
		in1=$*_LS_aper.fits suffix1= values1=id \
		in2=$*_LS_staraper.fits suffix2=_STAR values2=id

%_LS_extaper.fits: addextflux.py %_LS_aper.fits
	python3 $^ $@

%_all_extflux.fits: %_all.fits %_LS_extaper.fits
	stilts tmatch2 out=$@ matcher=exact find=best join=all1 fixcols=dups \
		in1=$*_all.fits values1=id \
		in2=$*_LS_extaper.fits values2=id  icmd2='keepcols "id apfluxext*_7 MW_TRANSMISSION_*"' \
		ocmd='addcol prior_GALflux_decam_g "(fracflux_g_LS<0.1&&decam_g_err>0)?(apfluxext_g_7/MW_TRANSMISSION_G)*pow(10, -28.44)*1e26:-99"' \
		ocmd='addcol prior_GALflux_decam_r "(fracflux_r_LS<0.1&&decam_r_err>0)?(apfluxext_r_7/MW_TRANSMISSION_R)*pow(10, -28.44)*1e26:-99"' \
		ocmd='addcol prior_GALflux_decam_z "(fracflux_z_LS<0.1&&decam_z_err>0)?(apfluxext_z_7/MW_TRANSMISSION_Z)*pow(10, -28.44)*1e26:-99"' \
		ocmd='addcol prior_GALflux_decam_g_errlo "(fracflux_g_LS<0.1&&decam_g_err>0)?(apfluxext_err_g_7/MW_TRANSMISSION_G)*pow(10, -28.44)*1e26:-99"' \
		ocmd='addcol prior_GALflux_decam_r_errlo "(fracflux_r_LS<0.1&&decam_r_err>0)?(apfluxext_err_r_7/MW_TRANSMISSION_R)*pow(10, -28.44)*1e26:-99"' \
		ocmd='addcol prior_GALflux_decam_z_errlo "(fracflux_z_LS<0.1&&decam_z_err>0)?(apfluxext_err_z_7/MW_TRANSMISSION_R)*pow(10, -28.44)*1e26:-99"' \
		ocmd='addcol prior_GALflux_decam_g_errhi "1e10"' \
		ocmd='addcol prior_GALflux_decam_r_errhi "1e10"' \
		ocmd='addcol prior_GALflux_decam_z_errhi "1e10"' \
		ocmd='delcols "apfluxext*_7 MW_TRANSMISSION_* id_2"' \
		ocmd='colmeta -name id "id_1"'

%_LS.fits: %_LS_aper.fits
	# adflux_* = adaptive flux column, depending on source type
	# for extended sources: 5'' aperture = OPT:apflux_*_7 and IR:apflux_*_2
	# for point sources: model flux
	# Milky way extinction correction
	# 28.44 is to convert to erg/s/cm^2/Hz
	stilts tpipe in=$^ out=$@ \
		cmd='addcol adflux_g "type==\"PSF\" ? flux_g : apflux_g_7"' \
		cmd='addcol adflux_r "type==\"PSF\" ? flux_r : apflux_r_7"' \
		cmd='addcol adflux_i "type==\"PSF\" ? flux_i : apflux_i_7"' \
		cmd='addcol adflux_z "type==\"PSF\" ? flux_z : apflux_z_7"' \
		cmd='addcol adflux_W1 "type==\"PSF\" ? flux_W1 : apflux_W1_2"' \
		cmd='addcol adflux_W2 "type==\"PSF\" ? flux_W2 : apflux_W2_2"' \
		cmd='addcol adflux_W3 "type==\"PSF\" ? flux_W3 : apflux_W3_2"' \
		cmd='addcol adflux_W4 "type==\"PSF\" ? flux_W4 : apflux_W4_2"' \
		cmd='addcol adflux_ivar_g  "type==\"PSF\" ? flux_ivar_g : apflux_ivar_g_7"' \
		cmd='addcol adflux_ivar_r  "type==\"PSF\" ? flux_ivar_r : apflux_ivar_r_7"' \
		cmd='addcol adflux_ivar_i  "type==\"PSF\" ? flux_ivar_i : apflux_ivar_i_7"' \
		cmd='addcol adflux_ivar_z  "type==\"PSF\" ? flux_ivar_z : apflux_ivar_z_7"' \
		cmd='addcol adflux_ivar_W1 "type==\"PSF\" ? flux_ivar_W1 : apflux_ivar_W1_2"' \
		cmd='addcol adflux_ivar_W2 "type==\"PSF\" ? flux_ivar_W2 : apflux_ivar_W2_2"' \
		cmd='addcol adflux_ivar_W3 "type==\"PSF\" ? flux_ivar_W3 : apflux_ivar_W3_2"' \
		cmd='addcol adflux_ivar_W4 "type==\"PSF\" ? flux_ivar_W4 : apflux_ivar_W4_2"' \
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

%_HSC_aper.fits: fetchHSC.py %.fits
	python3 $^ $@

%_HSC.fits: %_HSC_aper.fits %.fits
	# adflux_* = adaptive flux column, depending on source type
	# for extended sources: 5'' aperture = OPT:apflux_*_7 and IR:apflux_*_2
	# for point sources: model flux
	# Milky way extinction correction
	# 28.44 is to convert to erg/s/cm^2/Hz
	stilts tmatch2 out=$@ matcher=sky find=best params=1 fixcols=all \
		in1=$*.fits values1="RA DEC" icmd1='keepcols "id RA DEC"' suffix1=_in \
		in2=$< values2="RA DEC" suffix2= \
		ocmd='addcol extended "g_extendedness_value==1||r_extendedness_value==1||i_extendedness_value==1||z_extendedness_value==1||y_extendedness_value==1"' \
		ocmd='addcol adflux_g "extended ? g_psfflux_flux : g_apertureflux_57_flux"' \
		ocmd='addcol adflux_r "extended ? r_psfflux_flux : r_apertureflux_57_flux"' \
		ocmd='addcol adflux_i "extended ? i_psfflux_flux : i_apertureflux_57_flux"' \
		ocmd='addcol adflux_z "extended ? z_psfflux_flux : z_apertureflux_57_flux"' \
		ocmd='addcol adflux_y "extended ? y_psfflux_flux : y_apertureflux_57_flux"' \
		ocmd='addcol adflux_g_err "extended ? g_psfflux_fluxerr : g_apertureflux_57_fluxerr"' \
		ocmd='addcol adflux_r_err "extended ? r_psfflux_fluxerr : r_apertureflux_57_fluxerr"' \
		ocmd='addcol adflux_i_err "extended ? i_psfflux_fluxerr : i_apertureflux_57_fluxerr"' \
		ocmd='addcol adflux_z_err "extended ? z_psfflux_fluxerr : z_apertureflux_57_fluxerr"' \
		ocmd='addcol adflux_y_err "extended ? y_psfflux_fluxerr : y_apertureflux_57_fluxerr"' \
		ocmd='addcol LU_flux_g "ADFLUX_G>-10000000000L ? ((ADFLUX_G/1000*(pow(10, -(-a_g)/(2.5)))))*pow(10, -28.44):-99."' \
		ocmd='addcol LU_flux_r "ADFLUX_R>-10000000000L ? ((ADFLUX_R/1000*(pow(10, -(-a_r)/(2.5)))))*pow(10, -28.44):-99."' \
		ocmd='addcol LU_flux_i "ADFLUX_I>-10000000000L ? ((ADFLUX_I/1000*(pow(10, -(-a_i)/(2.5)))))*pow(10, -28.44):-99."' \
		ocmd='addcol LU_flux_z "ADFLUX_Z>-10000000000L ? ((ADFLUX_Z/1000*(pow(10, -(-a_z)/(2.5)))))*pow(10, -28.44):-99."' \
		ocmd='addcol LU_flux_y "ADFLUX_Y>-10000000000L ? ((ADFLUX_Y/1000*(pow(10, -(-a_y)/(2.5)))))*pow(10, -28.44):-99."' \
		ocmd='addcol LU_flux_g_err "ADFLUX_G_err>-10000000000L ? ((ADFLUX_G_err/1000*(pow(10, -(-a_g)/(2.5)))))*pow(10, -28.44):-99."' \
		ocmd='addcol LU_flux_r_err "ADFLUX_R_err>-10000000000L ? ((ADFLUX_R_err/1000*(pow(10, -(-a_r)/(2.5)))))*pow(10, -28.44):-99."' \
		ocmd='addcol LU_flux_i_err "ADFLUX_I_err>-10000000000L ? ((ADFLUX_I_err/1000*(pow(10, -(-a_i)/(2.5)))))*pow(10, -28.44):-99."' \
		ocmd='addcol LU_flux_z_err "ADFLUX_Z_err>-10000000000L ? ((ADFLUX_Z_err/1000*(pow(10, -(-a_z)/(2.5)))))*pow(10, -28.44):-99."' \
		ocmd='addcol LU_flux_y_err "ADFLUX_Y_err>-10000000000L ? ((ADFLUX_Y_err/1000*(pow(10, -(-a_y)/(2.5)))))*pow(10, -28.44):-99."' \
		ocmd='delcols "RA_in DEC_in"; colmeta -name id id_in'

%_GALEX.fits: %_coords.csv
	# NUV measurements are only valid down to NUVmag<20.8 AB mag NFlux>17.4
	# FUV measurements are only valid down to FUVmag<19.9 AB mag FFlux>39.7
	# _LU_A: convert energy fluxes [uJy] to photon fluxes in lephare units (LU) (factors of 1e-17/18, per Angstrom?) (?)
	# _LU: convert units to lephare-units
	#    some factors used here:
	#    (1.26727*pow(10, -17)) = 1e-29 / (3.34*pow(10, -19))*(pow(1538.6, 2))
	#    (5.59444*pow(10, -18)) = 1e-29 / (3.34*pow(10, -19))*(pow(2315.7, 2))
	# _real_LU_A: correct by milky way absorption
	# _real_LU: convert photon fluxes back to energy fluxes conversion (?)
	# 
	# no upper limits are provided for GALEX, because
	#   only no-entry cases when matching to this catalog would be suitable
	#   in the milky way, the upper limit would need to be corrected up, and thus be drastically weaker (higher)
	#   GALEX is very shallow already, so those upper limits are probably not so informative
	stilts cdsskymatch cdstable=II/335/galex_ais  \
		find=each in=$< ra=RA dec=DEC \
		radius=1 out=$@ \
		ocmd='addcol EBmV "E(B-V)"' \
		ocmd='addcol A_v "EBmV>0 ? (3.1*EBmV): -99."' \
		ocmd='addcol A_lam_fuv "A_v>0 ? A_v*(-0.3459882441104826 + 9.173790513029425/3.1): -99."' \
		ocmd='addcol A_lam_nuv "A_v>0 ? A_v*(0.6743929174392934 + 1.2480067831962458/3.1): -99."' \
		\
		ocmd='addcol FUV "FUVmag>0&&FUVmag<20.8 ? FUVmag-(8.286*EBmV): -99."' \
		ocmd='addcol NUV "NUVmag>0&&NUVmag<19.9 ? NUVmag-(8.621*EBmV): -99."' \
		ocmd='addcol e_FUV "FUVmag>0&&FUVmag<20.8 ?e_FUVmag:-99"' \
		ocmd='addcol e_NUV "NUVmag>0&&NUVmag<19.9?e_NUVmag:-99"' \
		\
		ocmd='addcol Fgood_extraction "e_Fflux>0&&Fafl==0&&Fexf==0"' \
		ocmd='addcol Ngood_extraction "e_Nflux>0&&Nafl==0&&Nexf==0"' \
		ocmd='addcol Fflux_LU "Fflux>39.7&&Fgood_extraction ? Fflux*(pow(10, -29)): -99."' \
		ocmd='addcol Nflux_LU "Nflux>17.4&&Ngood_extraction ? Nflux*(pow(10, -29)): -99."' \
		ocmd='addcol e_Fflux_LU "Fflux>39.7&&Fgood_extraction ? e_Fflux*(pow(10, -29)): -99."' \
		ocmd='addcol e_Nflux_LU "Nflux>17.4&&Ngood_extraction ? e_Nflux*(pow(10, -29)): -99."'  \
		ocmd='addcol Fflux_LU_A "Fflux>39.7&&Fgood_extraction ? Fflux*(1.26727*pow(10, -17)): -99."' \
		ocmd='addcol Nflux_LU_A "Nflux>17.4&&Ngood_extraction ? Nflux*(5.59444*pow(10, -18)): -99."' \
		ocmd='addcol e_Fflux_LU_A "Fflux>39.7&&Fgood_extraction ? e_Fflux*(1.26727*pow(10, -17)): -99."' \
		ocmd='addcol e_Nflux_LU_A "Nflux>17.4&&Ngood_extraction ? e_Nflux*(5.59444*pow(10, -18)): -99."' \
		\
		ocmd='addcol Fflux_real_LU_A "Fflux_LU_A>0 ? Fflux_LU_A*(pow(10, (0.4*A_lam_fuv))): -99."' \
		ocmd='addcol Nflux_real_LU_A "Nflux_LU_A>0 ? Nflux_LU_A*(pow(10, (0.4*A_lam_nuv))): -99."' \
		ocmd='addcol e_Fflux_real_LU_A "e_Fflux_LU_A>0 ? e_Fflux_LU_A*(pow(10, (0.4*A_lam_fuv))): -99."' \
		ocmd='addcol e_Nflux_real_LU_A "e_Nflux_LU_A>0 ? e_Nflux_LU_A*(pow(10, (0.4*A_lam_nuv))): -99."' \
		ocmd='addcol Fflux_real_LU "Fflux_real_LU_A>0 ? Fflux_real_LU_A*(3.34*pow(10, -19))*(pow(1538.6, 2)): -99."' \
		ocmd='addcol Nflux_real_LU "Nflux_real_LU_A>0 ? Nflux_real_LU_A*(3.34*pow(10, -19))*(pow(2315.7, 2)): -99."' \
		ocmd='addcol e_Fflux_real_LU "e_Fflux_real_LU_A>0 ? e_Fflux_real_LU_A*(3.34*pow(10, -19))*(pow(1538.6, 2)): -99."' \
		ocmd='addcol e_Nflux_real_LU "e_Nflux_real_LU_A>0 ? e_Nflux_real_LU_A*(3.34*pow(10, -19))*(pow(2315.7, 2)): -99."'

galex_ais_ctrs_ebv.fits: galexebv.py galex_ais_ctrs.fits
	python3 $<

%_GALEX_UL.fits: %_coords.csv galex_ais_ctrs_ebv.fits
	# search within 0.56 degrees
	stilts tmatch2 out=$@ matcher=sky find=best params=1980 \
		in1=$< values1="RA Dec" \
		in2=galex_ais_ctrs_ebv.fits values2="RAfdeg DEfdeg" \
		ocmd='addcol NUV_fluxlim "0.02 * pow(10, EBV_min * 1.35)"'
		#ocmd='addcol NUV_maglim "20.8 - 8.621 * EBV_min"' \
		#ocmd='addcol NUV_fluxlim "3.631 * pow(10, (NUV_maglim - 22.5) / 2.5)"'


%_ALLWISE_all.fits: %_coords.csv
	# match to AllWISE
	# convert Vega magnitudes to mJy fluxes
	# AllWISE is providing detections based on WISE data alone
	#  and is helpful for providing upper limits, by finding the brightest
	stilts cdsskymatch cdstable=II/328/allwise  \
		find=all in=$< ra=RA dec=DEC \
		radius=20 out=$@ \
		ocmd='addcol WISE1 "W1mag>0 ? (pow(10, -(W1mag+2.699+48.60)/(2.5)))*1e26: -99."' \
		ocmd='addcol WISE1_err "e_W1mag>0 ? (pow(10, -(W1mag-e_W1mag+2.699+48.60)/(2.5)))*1e26-WISE1: -99."' \
		ocmd='addcol WISE2 "W2mag>0 ? (pow(10, -(W2mag+3.339+48.60)/(2.5)))*1e26: -99."' \
		ocmd='addcol WISE2_err "e_W2mag>0 ? (pow(10, -(W2mag-e_W2mag+3.339+48.60)/(2.5)))*1e26-WISE2: -99."' \
		ocmd='addcol WISE3 "W3mag>0 ? (pow(10, -(W3mag+5.174+48.60)/(2.5)))*1e26: -99."' \
		ocmd='addcol WISE3_err "e_W3mag>0 ? (pow(10, -(W3mag-e_W3mag+5.174+48.60)/(2.5)))*1e26-WISE3: -99."' \
		ocmd='addcol WISE4 "W4mag>0 ? (pow(10, -(W4mag+6.620+48.60)/(2.5)))*1e26: -99."' \
		ocmd='addcol WISE4_err "e_W4mag>0 ? (pow(10, -(W4mag-e_W4mag+6.620+48.60)/(2.5)))*1e26-WISE4: -99."' \

%_ALLWISE_sum.fits: groupsumWISE.py %_coords.csv %_ALLWISE_all.fits
	python3 $^ $@

%_ALLWISE.fits: %_coords.csv
	# match to AllWISE
	# convert Vega magnitudes to mJy fluxes
	# AllWISE is providing detections based on WISE data alone
	#  and is helpful for providing upper limits
	stilts cdsskymatch cdstable=II/328/allwise  \
		find=each in=$< ra=RA dec=DEC \
		radius=4 out=$@ \
		ocmd='addcol WISE1 "W1mag>0 ? (pow(10, -(W1mag+2.699+48.60)/(2.5)))*1e26: -99."' \
		ocmd='addcol WISE1_err "e_W1mag>0 ? (pow(10, -(W1mag-e_W1mag+2.699+48.60)/(2.5)))*1e26-WISE1: -99."' 
		ocmd='addcol WISE2 "W2mag>0 ? (pow(10, -(W2mag+3.339+48.60)/(2.5)))*1e26: -99."' \
		ocmd='addcol WISE2_err "e_W2mag>0 ? (pow(10, -(W2mag-e_W2mag+3.339+48.60)/(2.5)))*1e26-WISE2: -99."' \
		ocmd='addcol WISE3 "W3mag>0 ? (pow(10, -(W3mag+5.174+48.60)/(2.5)))*1e26: -99."' \
		ocmd='addcol WISE3_err "e_W3mag>0 ? (pow(10, -(W3mag-e_W3mag+5.174+48.60)/(2.5)))*1e26-WISE3: -99."' \
		ocmd='addcol WISE4 "W4mag>0 ? (pow(10, -(W4mag+6.620+48.60)/(2.5)))*1e26: -99."' \
		ocmd='addcol WISE4_err "e_W4mag>0 ? (pow(10, -(W4mag-e_W4mag+6.620+48.60)/(2.5)))*1e26-WISE4: -99."' \


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

%_VHS_manual.fits: %_VHS.csv
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

%_VHS_noirlab.fits: fetchVHS.py %.fits
	python3 $^ $@

%_VHS.fits: %_VHS_noirlab.fits
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

%_UKIDSS_noirlab.fits: fetchUKIDSS.py %.fits
	python3 $^ $@

%_UKIDSS_manual.fits: %_UKIDSS.csv
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

%_UKIDSS.fits: %_UKIDSS_noirlab.fits
	# for point sources: fluxes within 2.8'' (apermag4)
	# for extended sources: fluxes within 5.6'' (apercorr6)
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


%_all_manual.fits: %.fits %_GALEX.fits %_LS.fits %_UKIDSS_manual.fits %_VHS_manual.fits %_ALLWISE.fits
	# merge everything together and use sensible column names
	# keep only WISE fluxes when ALLWISE also has a detection there
	# and if there are no blending issues
	# check LS10 fitbits for issues. bits 1, 2, 3, 7, 9, 12 are acceptable.
	# otherwise delete
	# take GALEX upper limits from eFEDS   0.02mJy = 20.8 mag; depth = 0.02mJy * pow(10, -3 + E(B-V)*1.35)
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
		ocmd='addcol FUV_err "e_Fflux_real_LU_GALEX*1e26"' \
		ocmd='addcol NUV "e_Nflux_real_LU_GALEX>0 ? Nflux_real_LU_GALEX*1e26 : GALEX"' \
		ocmd='addcol NUV_err "e_Nflux_real_LU_GALEX>0 ? e_Nflux_real_LU_GALEX*1e26 : -0.02"' \
		ocmd='addcol decam_g "(fracflux_g_LS>0.9?LU_flux_g_LS*1e26:-99)"' \
		ocmd='addcol decam_r "(fracflux_r_LS>0.9?LU_flux_r_LS*1e26:-99)"' \
		ocmd='addcol decam_i "(fracflux_i_LS>0.9?LU_flux_i_LS*1e26:-99)"' \
		ocmd='addcol decam_z "(fracflux_z_LS>0.9?LU_flux_z_LS*1e26:-99)"' \
		ocmd='addcol decam_g_err "(fracflux_g_LS>0.9?LU_flux_g_err_LS*1e26:-99)"' \
		ocmd='addcol decam_r_err "(fracflux_r_LS>0.9?LU_flux_r_err_LS*1e26:-99)"' \
		ocmd='addcol decam_i_err "(fracflux_i_LS>0.9?LU_flux_i_err_LS*1e26:-99)"' \
		ocmd='addcol decam_z_err "(fracflux_z_LS>0.9?LU_flux_z_err_LS*1e26:-99)"' \
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
		ocmd='addcol goodfitsLS "((fitbits_LS & (1 | 4 | 8)) == 0)"' \
		ocmd='addcol isolatedLS "(max(fracflux_g_LS,fracflux_r_LS,fracflux_z_LS,fracflux_w1_LS,fracflux_w2_LS)<0.1&&max(fracflux_w3_LS,fracflux_w4_LS)<10)"' \
		ocmd='addcol W34_blended "max(fracflux_w1_LS,fracflux_w2_LS)>0.1||max(fracflux_w4_LS,fracflux_w3_LS)>1"' \
		ocmd='addcol WISE1 "(isolated_LS ? LU_flux_w1_LS * 1e26 : (WISE1_ERR_ALLWISE > 0 ? WISE1_ALLWISE : -99))"' \
		ocmd='addcol WISE1_err "(isolated_LS ? LU_flux_w1_err_LS * 1e26 : (WISE1_ERR_ALLWISE > 0 ? -WISE1_ERR_ALLWISE : -99))"' \
		ocmd='addcol WISE2 "(isolated_LS ? LU_flux_w2_LS * 1e26 : (WISE2_ERR_ALLWISE > 0 ? WISE2_ALLWISE : -99))"' \
		ocmd='addcol WISE2_err "(isolated_LS ? LU_flux_w2_err_LS * 1e26 : (WISE2_ERR_ALLWISE > 0 ? -WISE2_ERR_ALLWISE : -99))"' \
		ocmd='addcol WISE3 "(isolated_LS&&!W34_blended ? LU_flux_w3_LS * 1e26 : (WISE3_ERR_ALLWISE > 0 ? WISE3_ALLWISE : -99))"' \
		ocmd='addcol WISE3_err "(isolated_LS&&!W34_blended ? LU_flux_w3_err_LS * 1e26 : (WISE3_ERR_ALLWISE > 0 ? -WISE3_ERR_ALLWISE : -99))"' \
		ocmd='addcol WISE4 "(isolated_LS&&!W34_blended ? LU_flux_w4_LS * 1e26 : (WISE4_ERR_ALLWISE > 0 ? WISE4_ALLWISE : -99))"' \
		ocmd='addcol WISE4_err "(isolated_LS&&!W34_blended ? LU_flux_w4_err_LS * 1e26 : (WISE4_ERR_ALLWISE > 0 ? -WISE4_ERR_ALLWISE : -99))"' \

%_all.fits: %.fits %_GALEX.fits %_LS.fits %_UKIDSS.fits %_VHS.fits %_ALLWISE_sum.fits %_GALEX_UL.fits
	# merge everything together and use sensible column names
	# keep only WISE fluxes when ALLWISE also has a detection there
	# and if there are no blending issues
	# check LS9 fitbits for issues. bits 1 | 4 | 8 | 128 | 512 | 4096 are problematic.
	#    otherwise delete bands
	# for check fracin_ columns, and unly use fluxes if large
	#     correct up error by these values, because we may miss 10% of the flux
	#     but we do not know about the exact galaxy morphology here compared to the coverage
	# for WISE, check if W3 or W4 are blended
	#    if they are, use 3sigma-flux * (3 * fracflux + 1) to estimate a conservative 3 sigma total flux of all sources
	#    and use that as a (3 sigma) upper limit
	stilts tmatchn nin=7 out=$@ \
		in1=$*.fits suffix1= values1=id \
		in2=$*_GALEX.fits suffix2=_GALEX values2=id \
		in3=$*_LS.fits suffix3=_LS values3=id \
		in4=$*_UKIDSS.fits suffix4=_UKIDSS values4=id \
		in5=$*_VHS.fits suffix5=_VHS values5=id \
		in6=$*_ALLWISE_sum.fits suffix6=_ALLWISE values6=id \
		in7=$*_GALEX_UL.fits suffix7=_GALEXUL values7=id \
		fixcols=all matcher=exact \
		ocmd='addcol pointlike "!(type_LS != \"PSF\")"' \
		ocmd='addcol inMzLSBASS "DEC>32&&RA>90&&RA<300"' \
		ocmd='addcol goodfitsLS "(fitbits_LS & (1 | 4 | 8)) == 0"' \
		ocmd='addcol isolatedLS "(max(fracflux_g_LS,fracflux_r_LS,fracflux_z_LS,fracflux_w1_LS,fracflux_w2_LS)<0.1&&max(fracflux_w3_LS,fracflux_w4_LS)<10)"' \
		ocmd='addcol W34_blended "max(fracflux_w1_LS,fracflux_w2_LS)>0.1||max(fracflux_w4_LS,fracflux_w3_LS)>1"' \
		ocmd='addcol FUV Fflux_real_LU_GALEX*1e26' \
		ocmd='addcol NUV "e_Nflux_real_LU_GALEX>0 ? Nflux_real_LU_GALEX*1e26 : 3 * NUV_fluxlim_GALEXUL"' \
		ocmd='addcol FUV_err "e_Fflux_real_LU_GALEX*1e26"' \
		ocmd='addcol NUV_err "e_Nflux_real_LU_GALEX>0 ? e_Nflux_real_LU_GALEX*1e26 : -NUV"' \
		ocmd='addcol decam_g "(!inMzLSBASS && fracin_g_LS>0.5?LU_flux_g_LS*1e26:-99)"' \
		ocmd='addcol decam_r "(!inMzLSBASS && fracin_r_LS>0.5?LU_flux_r_LS*1e26:-99)"' \
		ocmd='addcol decam_i "(fracin_i_LS>0.5?LU_flux_i_LS*1e26:-99)"' \
		ocmd='addcol decam_z "(!inMzLSBASS && fracin_z_LS>0.5?LU_flux_z_LS*1e26:-99)"' \
		ocmd='addcol decam_g_err "(!inMzLSBASS && fracin_g_LS>0.5?LU_flux_g_err_LS*1e26/fracin_g_LS:-99)"' \
		ocmd='addcol decam_r_err "(!inMzLSBASS && fracin_r_LS>0.5?LU_flux_r_err_LS*1e26/fracin_r_LS:-99)"' \
		ocmd='addcol decam_i_err "(fracin_i_LS>0.5?LU_flux_i_err_LS*1e26/fracin_i_LS:-99)"' \
		ocmd='addcol decam_z_err "(!inMzLSBASS && fracin_z_LS>0.5?LU_flux_z_err_LS*1e26/fracin_z_LS:-99)"' \
		ocmd='addcol 90prime_r "(inMzLSBASS && fracin_r_LS>0.5?LU_flux_r_LS*1e26:-99)"' \
		ocmd='addcol 90prime_g "(inMzLSBASS && fracin_g_LS>0.5?LU_flux_g_LS*1e26:-99)"' \
		ocmd='addcol zd_mosaic "(inMzLSBASS && fracin_z_LS>0.5?LU_flux_z_LS*1e26:-99)"' \
		ocmd='addcol 90prime_r_err "(inMzLSBASS && fracin_r_LS>0.5?LU_flux_r_err_LS*1e26/fracin_r_LS:-99)"' \
		ocmd='addcol zd_mosaic_err "(inMzLSBASS && fracin_z_LS>0.5?LU_flux_z_err_LS*1e26/fracin_z_LS:-99)"' \
		ocmd='addcol 90prime_g_err "(inMzLSBASS && fracin_g_LS>0.5?LU_flux_g_err_LS*1e26/fracin_g_LS:-99)"' \
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
		ocmd='addcol WISE1_origin "isolatedLS ? \"LS10\" : ( (fracflux_w1_LS > 1 || WISE1_ALLWISE < LU_flux_w1_LS * 1e26 * (1 + fracflux_w1_LS)) ? \"AllWISE\" : \"LS10UL\")"' \
		ocmd='addcol WISE1        "isolatedLS ? LU_flux_w1_LS * 1e26 : ((fracflux_w1_LS > 1 || WISE1_ALLWISE < LU_flux_w1_LS * 1e26 * (1 + fracflux_w1_LS)) ? WISE1_ALLWISE : LU_flux_w1_LS * 1e26)"' \
		ocmd='addcol WISE1_err "isolatedLS ? LU_flux_w1_err_LS * 1e26 : -WISE1"' \
		ocmd='addcol WISE2_origin "isolatedLS ? \"LS10\" : ( (fracflux_w2_LS > 1 || WISE2_ALLWISE < LU_flux_w2_LS * 1e26 * (1 + fracflux_w2_LS)) ? \"AllWISE\" : \"LS10UL\")"' \
		ocmd='addcol WISE2 "isolatedLS ? LU_flux_w2_LS * 1e26 : ((fracflux_w2_LS > 1 || WISE2_ALLWISE < LU_flux_w2_LS * 1e26 * (1 + fracflux_w2_LS)) ? WISE2_ALLWISE : LU_flux_w2_LS * 1e26)"' \
		ocmd='addcol WISE2_err "isolatedLS ? LU_flux_w2_err_LS * 1e26 : -WISE2"' \
		ocmd='addcol WISE34_origin "isolatedLS&&!W34_blended ? \"LS10\" : \"AllWISE\""' \
		ocmd='addcol WISE3 "isolatedLS&&!W34_blended ? LU_flux_w3_LS*1e26 : WISE3_ALLWISE"' \
		ocmd='addcol WISE3_err "isolatedLS&&!W34_blended ? LU_flux_w3_err_LS*1e26 : -WISE3"' \
		ocmd='addcol WISE4 "isolatedLS&&!W34_blended ? LU_flux_w4_LS*1e26 : WISE4_ALLWISE"' \
		ocmd='addcol WISE4_err "isolatedLS&&!W34_blended ? LU_flux_w4_err_LS*1e26 : -WISE4"' \

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
		ocmd='addcol pointlike "!(extended_HSC)"' \
		ocmd='addcol inMzLSBASS "DEC>32&&RA>90&&RA<300"' \
		ocmd='addcol goodfitsLS "(fitbits_LS & (1 | 4 | 8)) == 0"' \
		ocmd='addcol isolatedLS "(max(fracflux_g_LS,fracflux_r_LS,fracflux_z_LS,fracflux_w1_LS,fracflux_w2_LS)<0.1&&max(fracflux_w3_LS,fracflux_w4_LS)<10)"' \
		ocmd='addcol W34_blended "max(fracflux_w1_LS,fracflux_w2_LS)>0.1||max(fracflux_w4_LS,fracflux_w3_LS)>1"' \
		ocmd='addcol FUV Fflux_real_LU_GALEX*1e26' \
		ocmd='addcol NUV "e_Nflux_real_LU_GALEX>0 ? Nflux_real_LU_GALEX*1e26 : 3 * NUV_fluxlim_GALEXUL"' \
		ocmd='addcol FUV_err "e_Fflux_real_LU_GALEX*1e26"' \
		ocmd='addcol NUV_err "e_Nflux_real_LU_GALEX>0 ? e_Nflux_real_LU_GALEX*1e26 : -NUV"' \
		ocmd='addcol decam_g "(!inMzLSBASS && fracin_g_LS>0.5?LU_flux_g_LS*1e26:-99)"' \
		ocmd='addcol decam_r "(!inMzLSBASS && fracin_r_LS>0.5?LU_flux_r_LS*1e26:-99)"' \
		ocmd='addcol decam_i "(fracin_i_LS>0.5?LU_flux_i_LS*1e26:-99)"' \
		ocmd='addcol decam_z "(!inMzLSBASS && fracin_z_LS>0.5?LU_flux_z_LS*1e26:-99)"' \
		ocmd='addcol decam_g_err "(!inMzLSBASS && fracin_g_LS>0.5?LU_flux_g_err_LS*1e26/fracin_g_LS:-99)"' \
		ocmd='addcol decam_r_err "(!inMzLSBASS && fracin_r_LS>0.5?LU_flux_r_err_LS*1e26/fracin_r_LS:-99)"' \
		ocmd='addcol decam_i_err "(fracin_i_LS>0.5?LU_flux_i_err_LS*1e26/fracin_i_LS:-99)"' \
		ocmd='addcol decam_z_err "(!inMzLSBASS && fracin_z_LS>0.5?LU_flux_z_err_LS*1e26/fracin_z_LS:-99)"' \
		ocmd='addcol 90prime_r "(inMzLSBASS && fracin_r_LS>0.5?LU_flux_r_LS*1e26:-99)"' \
		ocmd='addcol 90prime_g "(inMzLSBASS && fracin_g_LS>0.5?LU_flux_g_LS*1e26:-99)"' \
		ocmd='addcol zd_mosaic "(inMzLSBASS && fracin_z_LS>0.5?LU_flux_z_LS*1e26:-99)"' \
		ocmd='addcol 90prime_r_err "(inMzLSBASS && fracin_r_LS>0.5?LU_flux_r_err_LS*1e26/fracin_r_LS:-99)"' \
		ocmd='addcol zd_mosaic_err "(inMzLSBASS && fracin_z_LS>0.5?LU_flux_z_err_LS*1e26/fracin_z_LS:-99)"' \
		ocmd='addcol 90prime_g_err "(inMzLSBASS && fracin_g_LS>0.5?LU_flux_g_err_LS*1e26/fracin_g_LS:-99)"' \
		ocmd='addcol HSC_g "(!g_psfflux_flag_HSC&&!g_apertureflux_57_flag_HSC?LU_flux_g_HSC*1e26:-99)"' \
		ocmd='addcol HSC_r "(!r_psfflux_flag_HSC&&!r_apertureflux_57_flag_HSC?LU_flux_r_HSC*1e26:-99)"' \
		ocmd='addcol HSC_i "(!i_psfflux_flag_HSC&&!i_apertureflux_57_flag_HSC?LU_flux_i_HSC*1e26:-99)"' \
		ocmd='addcol HSC_z "(!z_psfflux_flag_HSC&&!z_apertureflux_57_flag_HSC?LU_flux_z_HSC*1e26:-99)"' \
		ocmd='addcol HSC_y "(!y_psfflux_flag_HSC&&!y_apertureflux_57_flag_HSC?LU_flux_y_HSC*1e26:-99)"' \
		ocmd='addcol HSC_g_err "(!g_psfflux_flag_HSC&&!g_apertureflux_57_flag_HSC?LU_flux_g_err_HSC*1e26:-99)"' \
		ocmd='addcol HSC_r_err "(!r_psfflux_flag_HSC&&!r_apertureflux_57_flag_HSC?LU_flux_r_err_HSC*1e26:-99)"' \
		ocmd='addcol HSC_i_err "(!i_psfflux_flag_HSC&&!i_apertureflux_57_flag_HSC?LU_flux_i_err_HSC*1e26:-99)"' \
		ocmd='addcol HSC_z_err "(!z_psfflux_flag_HSC&&!z_apertureflux_57_flag_HSC?LU_flux_z_err_HSC*1e26:-99)"' \
		ocmd='addcol HSC_y_err "(!y_psfflux_flag_HSC&&!y_apertureflux_57_flag_HSC?LU_flux_y_err_HSC*1e26:-99)"' \
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
		ocmd='addcol WISE1_origin "isolatedLS ? \"LS10\" : ( (fracflux_w1_LS > 1 || WISE1_ALLWISE < LU_flux_w1_LS * 1e26 * (1 + fracflux_w1_LS)) ? \"AllWISE\" : \"LS10UL\")"' \
		ocmd='addcol WISE1        "isolatedLS ? LU_flux_w1_LS * 1e26 : ((fracflux_w1_LS > 1 || WISE1_ALLWISE < LU_flux_w1_LS * 1e26 * (1 + fracflux_w1_LS)) ? WISE1_ALLWISE : LU_flux_w1_LS * 1e26)"' \
		ocmd='addcol WISE1_err "isolatedLS ? LU_flux_w1_err_LS * 1e26 : -WISE1"' \
		ocmd='addcol WISE2_origin "isolatedLS ? \"LS10\" : ( (fracflux_w2_LS > 1 || WISE2_ALLWISE < LU_flux_w2_LS * 1e26 * (1 + fracflux_w2_LS)) ? \"AllWISE\" : \"LS10UL\")"' \
		ocmd='addcol WISE2 "isolatedLS ? LU_flux_w2_LS * 1e26 : ((fracflux_w2_LS > 1 || WISE2_ALLWISE < LU_flux_w2_LS * 1e26 * (1 + fracflux_w2_LS)) ? WISE2_ALLWISE : LU_flux_w2_LS * 1e26)"' \
		ocmd='addcol WISE2_err "isolatedLS ? LU_flux_w2_err_LS * 1e26 : -WISE2"' \
		ocmd='addcol WISE34_origin "isolatedLS&&!W34_blended ? \"LS10\" : \"AllWISE\""' \
		ocmd='addcol WISE3 "isolatedLS&&!W34_blended ? LU_flux_w3_LS*1e26 : WISE3_ALLWISE"' \
		ocmd='addcol WISE3_err "isolatedLS&&!W34_blended ? LU_flux_w3_err_LS*1e26 : -WISE3"' \
		ocmd='addcol WISE4 "isolatedLS&&!W34_blended ? LU_flux_w4_LS*1e26 : WISE4_ALLWISE"' \
		ocmd='addcol WISE4_err "isolatedLS&&!W34_blended ? LU_flux_w4_err_LS*1e26 : -WISE4"' \

%_lite.fits: %.fits
	stilts tpipe in=$^ out=$@ cmd='delcols "adflux*_LS *_LS LU_flux*_LS *_GALEX *_VHS *_UKIDSS *_ALLWISE *_GALEXUL *_HSC"'
	stilts tpipe in=$@ omode=stats

%.fits_errors.pdf %.fits_fluxes.pdf: %.fits
	python3 checkphotometry.py $<

.SECONDARY:
.PRECIOUS:

