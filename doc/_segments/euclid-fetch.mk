%_EuclidQ1_noirlab.fits: %_coords.csv
	datalab mydb_import --table=rbeuclidq1 --data=$^ --append=False
	datalab query --sql="select id, euclid.* from euclid_q1.object as euclid, mydb://rbeuclidq1 as input WHERE q3c_join(input.ra, input.dec, euclid.right_ascension, euclid.declination, 1./60/60)" ${DLLONGQARGS} --out=mydb://rbeuclidq1out
	datalab query --sql="select * from mydb://rbeuclidq1out" ${DLLONGQARGS} --fmt=fits --out=$@
	stilts tpipe in=$@ ifmt=fits omode=count
	datalab mydb_drop --table=rbeuclidq1 || true
	datalab mydb_drop --table=rbeuclidq1out || true
