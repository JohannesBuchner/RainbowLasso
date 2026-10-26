%_VHS_noirlab.fits: %_coords.csv
	datalab mydb_import --table=rblvhs --data=$^ --append=False
	datalab query --sql="select id, vhs.* from vhs_dr5.vhs_cat_v3 as vhs, mydb://rblvhs WHERE q3c_join(ra, dec, vhs.ra2000, vhs.dec2000, 1./60/60)" ${DLLONGQARGS} --out=mydb://rblvhsout
	datalab query --sql="select * from mydb://rblvhsout" ${DLLONGQARGS} --fmt=fits --out=$@
	stilts tpipe in=$@ ifmt=fits omode=count
	datalab mydb_drop --table=rblvhs || true
	datalab mydb_drop --table=rblvhsout || true
