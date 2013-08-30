#!/usr/bin/env python

import argparse # Requires Python 2.7+
import numpy as np

# Parse command line arguments
description = """Strip headers and unwanted columns from WMAP TT Power Spectrum data ready for LAPD
"""

parser = argparse.ArgumentParser(description=description)
parser.add_argument(
    'incat', type=str,
    help='ASCII catalogue containing TT power spectra, with a column layout following the format '+
    'of wmap_tt_spectrum_7yr_v4p1.txt')
parser.add_argument(
    'outcat', type=str,
    help='Filename for ASCII output after stripping header and unwanted 3rd and 4th columns to '
    'generate an input file ready for the LAPD program')
args = parser.parse_args()

# Read in numeric data (ignoring headers) and write out the three columns to the desired outcat
data = np.loadtxt(args.incat)
np.savetxt(args.outcat, np.array((data[:, 0], data[:, 1], data[:, 2])).T, fmt='%7d %12.4f %12.4f')
# done!
