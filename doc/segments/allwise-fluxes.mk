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
