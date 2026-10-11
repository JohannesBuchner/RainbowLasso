"""Smoke test for the lite catalogue filtering.

Applies the %_lite.fits column filtering (see Makefile) to the committed
selection catalogue with astropy, writes the result and checks it. This
avoids the stilts and NOIRLab (astro-datalab) dependencies of the full
data pipeline, so `make test` can run anywhere.
"""
import fnmatch
import os
import tempfile

from astropy.table import Table

# same column patterns as the stilts delcols call in the %_lite.fits rule
LITE_DELCOLS = (
	'adflux*_LS *_LS LU_flux*_LS *_GALEX *_VHS *_UKIDSS '
	'*_ALLWISE *_GALEXUL *_HSC_extaper *_HSC psf*_SDSS aper*_SDSS '
	'id*_SDSS ra*_SDSS dec*_SDSS'
).split()


def main():
	table = Table.read('dr16QWX_selection.fits')
	nrows = len(table)
	columns = list(table.colnames)

	drop = [c for c in columns
		if any(fnmatch.fnmatch(c, p) for p in LITE_DELCOLS)]
	keep = [c for c in columns if c not in drop]
	for core in ('id', 'RA', 'DEC'):
		assert core in keep, 'core column %s must survive the lite filtering' % core
	assert len(keep) > 0, 'lite filtering must not remove all columns'

	table.remove_columns(drop)
	out = os.path.join(tempfile.mkdtemp(), 'dr16QWX_selection_lite_test.fits')
	table.write(out, overwrite=True)

	check = Table.read(out)
	assert len(check) == nrows, 'row count must be preserved'
	assert len(check.columns) == len(keep), \
		'column count must match the filtered set'

	print('input : %d rows, %d columns' % (nrows, len(columns)))
	print('lite  : %d rows, %d columns (dropped %d)' % (nrows, len(keep), len(drop)))
	print('OK: lite filtering round-trips through astropy')


if __name__ == '__main__':
	main()
