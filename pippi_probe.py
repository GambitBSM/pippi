
#############################################################
# pippi: parse it, plot it
# ------------------------
# Simple HDF5 probing program for pippi.
#
# Author: Pat Scott (p.scott@imperial.ac.uk)
# Originally developed: Oct 2015
#############################################################

from pippi_utils import *

# Define probe-specific pip file entries
hdf5_cols = dataObject('assign_hdf5_label_to_column',string_dictionary)
keys = keys+[hdf5_cols]

#input:   filename = the name of the pip file
def probe(filename):
  #Parse pip file
  getIniData(filename,keys)
  # Open main chain and read in contents
  mainArray = getChainData(mainChain.value, hdf5_assignments=hdf5_cols, probe_only=True)
