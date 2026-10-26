%_UKIDSS_noirlab.fits: %_coords.csv
	datalab mydb_import --table=rblukidss --data=$^ --append=False
	datalab query --sql="select id, ukidss.* from ukidss_dr11plus.lassource as ukidss, mydb://rblukidss as input WHERE q3c_join(input.ra, input.dec, ukidss.ra, ukidss.dec, 1./60/60)" ${DLLONGQARGS} --out=mydb://rblukidssout
	datalab query --sql="select * from mydb://rblukidssout" ${DLLONGQARGS} --fmt=fits --out=$@
	stilts tpipe in=$@ ifmt=fits omode=count
	datalab mydb_drop --table=rblukidss || true
	datalab mydb_drop --table=rblukidssout || true
