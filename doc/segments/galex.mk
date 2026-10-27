%_GALEX.fits: %_coords.csv
	# NUV measurements are only valid down to NUVmag<20.8 AB mag NFlux>17.4
	# FUV measurements are only valid down to FUVmag<19.9 AB mag FFlux>39.7
	# _LU_A: convert energy fluxes [uJy] to photon fluxes in lephare units (LU) (factors of 1e-17/18, per Angstrom?) (?)
	# _LU: convert units to lephare-units (not used further)
	#    some factors used here:
	#    (1.26727*pow(10, -17)) = 1e-29 / (3.34*pow(10, -19))*(pow(1538.6, 2))
	#    (5.59444*pow(10, -18)) = 1e-29 / (3.34*pow(10, -19))*(pow(2315.7, 2))
	# _real_LU_A: correct _LU_A by milky way absorption
	# _real_LU: convert photon fluxes back to energy fluxes conversion (?)
	# conversion from Vega to AB: 18.82, 20.08 (https://asd.gsfc.nasa.gov/archive/galex/FAQ/counts_background.html)
	# no upper limits are provided for GALEX, because
	#   only no-entry cases when matching to this catalog would be suitable
	#   in the milky way, the upper limit would need to be corrected up, and thus be drastically weaker (higher)
	#   GALEX is very shallow already, so those upper limits are probably not so informative
	# FUV_4 is an aperture with 4 pixel radius, pixel size is 1.5'' -> 6'' radius: very large!
	# FWHM of GALEX is very large, best to treat all objects at z>>0 as point sources
	stilts cdsskymatch cdstable=II/335/galex_ais  \
		find=each in=$< ra=RA dec=DEC \
		radius=4 out=$@ \
		ocmd='fixcolnames' \
		ocmd='colmeta -name EBmV "E_B_V_"' \
		ocmd='addcol A_v "EBmV>0 ? (3.1*EBmV): -99."' \
		ocmd='addcol A_lam_fuv "A_v>0 ? A_v*(-0.3459882441104826 + 9.173790513029425/3.1): -99."' \
		ocmd='addcol A_lam_nuv "A_v>0 ? A_v*(0.6743929174392934 + 1.2480067831962458/3.1): -99."' \
		\
		ocmd='addcol FUV "FUVmag>0&&FUVmag<19.9 ? FUVmag-(8.286*EBmV): -99."' \
		ocmd='addcol NUV "NUVmag>0&&NUVmag<20.8 ? NUVmag-(8.621*EBmV): -99."' \
		ocmd='addcol e_FUV "FUVmag>0&&FUVmag<19.9 ?e_FUVmag:-99"' \
		ocmd='addcol e_NUV "NUVmag>0&&NUVmag<20.8 ?e_NUVmag:-99"' \
		\
		ocmd='addcol FUV_aper "FUVmag>0&&FUVmag<19.9 ? FUV_4+18.82-(8.286*EBmV): -99."' \
		ocmd='addcol NUV_aper "NUVmag>0&&NUVmag<20.8 ? NUV_4+20.08-(8.621*EBmV): -99."' \
		ocmd='addcol e_FUV_aper "FUVmag>0&&FUVmag<19.9 ? e_FUV_4:-99"' \
		ocmd='addcol e_NUV_aper "NUVmag>0&&NUVmag<20.8 ? e_NUV_4:-99"' \
		\
		ocmd='addcol Fflux_aper "e_FUV_aper>0 ? pow(10, -0.4*FUV_aper)*1e6*3631 : -99."' \
		ocmd='addcol Nflux_aper "e_NUV_aper>0 ? pow(10, -0.4*NUV_aper)*1e6*3631 : -99."' \
		ocmd='addcol e_Fflux_aper "e_FUV_aper>0 ? pow(10, -0.4*(FUV_aper-e_FUV_aper))*1e6*3631 - Fflux_aper : -99"' \
		ocmd='addcol e_Nflux_aper "e_NUV_aper>0 ? pow(10, -0.4*(NUV_aper-e_NUV_aper))*1e6*3631 - Nflux_aper : -99"' \
		\
		ocmd='addcol Fgood_extraction "e_Fflux>0&&Fafl==0&&Fexf==0"' \
		ocmd='addcol Ngood_extraction "e_Nflux>0&&Nafl==0&&Nexf==0"' \
		ocmd='addcol Fflux_real_LU "Fflux>39.7&&Fgood_extraction ? Fflux*pow(10, (0.4*A_lam_fuv))*1e-29: -99."' \
		ocmd='addcol Nflux_real_LU "Nflux>17.4&&Ngood_extraction ? Nflux*pow(10, (0.4*A_lam_nuv))*1e-29: -99."' \
		ocmd='addcol e_Fflux_real_LU "Fflux>39.7&&Fgood_extraction ? e_Fflux*pow(10, (0.4*A_lam_fuv))*1e-29: -99."' \
		ocmd='addcol e_Nflux_real_LU "Nflux>17.4&&Ngood_extraction ? e_Nflux*pow(10, (0.4*A_lam_nuv))*1e-29: -99."' \
		ocmd='addcol Fflux_aper_real_LU "Fflux_aper>39.7&&Fgood_extraction ? Fflux_aper*pow(10, (0.4*A_lam_fuv))*1e-29: -99."' \
		ocmd='addcol Nflux_aper_real_LU "Nflux_aper>17.4&&Ngood_extraction ? Nflux_aper*pow(10, (0.4*A_lam_nuv))*1e-29: -99."' \
		ocmd='addcol e_Fflux_aper_real_LU "Fflux_aper>39.7&&Fgood_extraction ? e_Fflux_aper*pow(10, (0.4*A_lam_fuv))*1e-29: -99."' \
		ocmd='addcol e_Nflux_aper_real_LU "Nflux_aper>17.4&&Ngood_extraction ? e_Nflux_aper*pow(10, (0.4*A_lam_nuv))*1e-29: -99."' \
