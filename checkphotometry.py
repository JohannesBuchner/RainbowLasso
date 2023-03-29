import matplotlib.pyplot as plt
from astropy.table import Table
import sys
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages

t = Table.read(sys.argv[1], memmap=True)
print(t.colnames)
fluxcols = [col for col in t.colnames if col + '_err' in t.colnames]
print("flux columns:", fluxcols)
sequential_pairs = list(zip(fluxcols[:-1], fluxcols[1:]))
similar_pairs = []
for col in fluxcols:
	if '_' in col:
		band = col.split('_')[1]
		similar_pairs += [(col, col2) for col2 in fluxcols if col2.endswith('_' + band) and col2 != col]


with PdfPages(sys.argv[1] + '_fluxes.pdf') as pdf_fluxes, PdfPages(sys.argv[1] + '_errors.pdf') as pdf_errs:
	for last_col, col in similar_pairs + sequential_pairs:
		print(" %s vs %s ..." % (col, last_col))
		last_vals = t[last_col]
		last_errs = t[last_col + '_err']
		vals = t[col]
		errs = t[col + '_err']
		if not np.any(vals > 0):
			continue

		if last_vals is not None:
			both_defined = np.logical_and(errs > 0, last_errs > 0)
			if both_defined.any():
				plt.figure()
				plt.scatter(last_vals[both_defined], vals[both_defined])
				ylo, yhi = vals[both_defined].min(), vals[both_defined].max()
				plt.plot([ylo, yhi], [ylo, yhi], ls='--', color='k')
				plt.xlabel(last_col)
				plt.ylabel(col)
				plt.xscale('log')
				plt.yscale('log')
				plt.tight_layout()
				pdf_fluxes.savefig()
				plt.close()

		plt.figure()
		plt.scatter(errs, vals / errs, label=col)
		if last_vals is not None:
			plt.scatter(last_errs, last_vals / last_errs, c='lightgray', alpha=0.3, label=last_col)
		#plt.ylim(0, None)
		plt.xlabel("Error (%s)" % (col))
		plt.ylabel("Value / Error (%s)" % (col))
		plt.xscale('log')
		plt.yscale('log')
		if (errs < 0).any():
			plt.text(0.98, 0,  
				'%.2f%% have err<0\n[%s-%s]' % (
				(errs<0).mean()*100, errs[errs<0].min(), errs[errs<0].max()), 
				va='bottom', ha='right', transform=plt.gca().transAxes)
		plt.legend(loc='best')
		pdf_errs.savefig()
		plt.close()

		last_col, last_vals, last_errs = col, vals, errs
