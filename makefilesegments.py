"""
Extracts Makefile rule segments into files for the sphinx documentation.

Each documentation page for a survey/catalog includes the relevant Makefile
rule (target + recipe) via sphinx literalinclude. This script parses the
Makefile and writes one file per segment into doc/segments/.

Usage::

	python3 makefilesegments.py                 # extract all segments
	python3 makefilesegments.py galex legacysurvey-fetch   # extract a subset
	python3 makefilesegments.py --list          # list available segments

The segment files are regenerated whenever the Makefile changes.
"""
import argparse
import os
import re
import sys

# mapping of documentation segment name -> Makefile target name
SEGMENTS = {
	'coords': '%_coords.csv',
	'galex': '%_GALEX.fits',
	'galex-coverage': '%_GALEX_UL.fits',
	'galexebv': 'galex_ais_ctrs_ebv.fits',
	'legacysurvey-fetch': '%_LS_noirlab.fits',
	'legacysurvey-fluxes': '%_LS.fits',
	'legacysurvey-extflux': '%_LS_extaper.fits',
	'legacysurvey-extflux-merge': '%_all_extflux.fits',
	'vhs-fetch': '%_VHS_noirlab.fits',
	'vhs-fluxes': '%_VHS.fits',
	'ukidss-fetch': '%_UKIDSS_noirlab.fits',
	'ukidss-fluxes': '%_UKIDSS.fits',
	'allwise-fetch': '%_ALLWISE_all.fits',
	'allwise-sum': '%_ALLWISE_sum.fits',
	'allwise-fluxes': '%_ALLWISE.fits',
	'final': '%_all.fits',
	'final-lite': '%_lite.fits',
	'hsc-fetch': '%_HSC_aper.fits',
	'hsc-fluxes': '%_HSC.fits',
	'hscall': '%_HSCall.fits',
	'hsc-extflux': '%_HSCall_extflux.fits',
	'sdss-fetch': '%_SDSS.fits',
	'sdss-fluxes': '%_withSDSS.fits',
	'euclid-fetch': '%_EuclidQ1_noirlab.fits',
	'euclid-prior': '%_all_extflux_withEuclid.fits',
}


def parse_rules(makefile_path):
	"""Parses the Makefile into rule blocks (target line + recipe).

	Returns a dict of target name -> full rule text (target line,
	prerequisites and recipe, including comment lines within the rule).
	"""
	rules = {}
	current = []
	for line in open(makefile_path):
		line = line.rstrip('\n')
		if line and not line[0].isspace():
			# line starting at column 0: begins a new rule or variable
			store_rule(rules, current)
			current = [line]
		else:
			current.append(line)
	store_rule(rules, current)
	return rules


def store_rule(rules, block):
	"""Adds a rule block to the rules dict, if it defines a target."""
	block = strip_empty_lines(block)
	if not block:
		return
	m = re.match(r'^([^:#\s]+)\s*:', block[0])
	if m:
		rules[m.group(1)] = '\n'.join(block)


def strip_empty_lines(block):
	"""Removes trailing empty lines from a rule block."""
	while block and not block[-1].strip():
		block = block[:-1]
	return block


def extract_segment(rules, target):
	"""Returns the rule text for the given target, or None if not found."""
	if target in rules:
		return rules[target]
	return None


def main():
	"""Writes the requested Makefile segments into doc/segments/."""
	parser = argparse.ArgumentParser(
		description="Extract Makefile rule segments for the sphinx docs.")
	parser.add_argument('segments', nargs='*', default=None,
		help="segment names to extract (default: all)")
	parser.add_argument('--makefile', default='Makefile',
		help="Makefile to parse (default: Makefile)")
	parser.add_argument('--outdir', default=os.path.join('doc', 'segments'),
		help="output directory (default: doc/segments)")
	parser.add_argument('--list', action='store_true',
		help="list available segments and exit")
	args = parser.parse_args()

	if args.list:
		for name in sorted(SEGMENTS):
			print('%-30s %s' % (name, SEGMENTS[name]))
		return

	rules = parse_rules(args.makefile)
	names = args.segments if args.segments else sorted(SEGMENTS)
	unknown = [name for name in names if name not in SEGMENTS]
	if unknown:
		print("unknown segments: %s (see --list)" % ', '.join(unknown),
			file=sys.stderr)
		sys.exit(1)

	os.makedirs(args.outdir, exist_ok=True)
	for name in names:
		target = SEGMENTS[name]
		text = extract_segment(rules, target)
		if text is None:
			print("target %s (%s) not found in %s" %
				(target, name, args.makefile), file=sys.stderr)
			sys.exit(1)
		outpath = os.path.join(args.outdir, name + '.mk')
		with open(outpath, 'w') as f:
			f.write(text + '\n')
		print("wrote %s  (target %s)" % (outpath, target))


if __name__ == '__main__':
	main()
