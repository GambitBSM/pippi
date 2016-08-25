#############################################################
# pippi: parse it, plot it
# ------------------------
# Utilities module for pippi.
#
# Author: Pat Scott (patscott@physics.mcgill.ca)
# Originally developed: March 2012
# Modified: Christoph Weniger, 2015
#############################################################

import sys
import subprocess
import os.path
import re
import datetime
import numpy as np
from pippi_colours import *

def mapToRefLike(otherlike,array):
  #Map from one permitted likelihood form to another
  return [permittedLikes[otherlike](member) for member in array]


def safe_dir(dirname):
  #Check that dirname exists, if not then make it
  if not os.path.exists(dirname):
    os.makedirs(os.path.abspath(dirname))


def getIniData(filename,keys,savekeys=None,savedir=None):
  #Read initialisation data from pip file filename into keys

  if (len(filename) != 1):
    usage()
    return
  filename = filename[0]

  #Try to open pip file
  pipfile = safe_open(filename)
  #Get the contents
  original_parse_options = pipfile.readlines()
  #Shut it
  pipfile.close
  #Remove all comments
  for i,line in enumerate(original_parse_options): original_parse_options[i] = re.sub(r';.*$', '', line)
  #Stitch together multi-line lines
  parse_options = list([])
  oldline = ""
  for i,line in enumerate(original_parse_options):
    newline = line.rstrip()
    oldline += newline
    if (oldline != "" and oldline[-1] == "\\"):
      oldline = oldline[:-1]
    else:
      parse_options.append(oldline+"\n")
      oldline = ""

  #Find relevant bits in pipfile
  for i,key in enumerate(keys):
    lines = filter(key.seek,parse_options)
    if (len(lines) == 0):
      sys.exit('Error: field '+key.pipFileKey+' required for requested operation not found in '+filename+'.  Quitting...\n')
    if (len(lines) > 1):
      sys.exit('Error: field '+key.pipFileKey+' required for requested operation duplicated in '+filename+'.  Quitting...\n')
    # Convert key entries to correct internal data format
    key.convert(re.sub(r'.*=', '', lines[0]))

  # Save requested keys
  if savekeys is not None:

    # Make sure saving is actually possible
    if not mainChain in keys:
      print '\n  Warning: saving of keys not possible because mainChain is undefined.\n  Skipping save...'
      return

    # Open the file keys will be saved to
    if savedir is None or not savedir in keys or savedir.value is None:
      outfile = smart_open(re.sub(r"\..?.?.?$", '', mainChain.value)+'_savedkeys.pip','w')
    else:
      # Make sure savedir exists, make it if not
      safe_dir(savedir.value)
      outfile = smart_open(savedir.value+'/'+re.sub(r'.*/|\..?.?.?$', '', mainChain.value)+'_savedkeys.pip','w')
    outfile.write(';This file of saved keys created by pippi '\
                +pippiVersion+' on '+datetime.datetime.now().strftime('%c')+'\n')

    # Parse the pip file again and save the requested keys
    for key in savekeys:
      lines = filter(key.seek,parse_options)
      # Save keys verbatim to the savedkeys pip file
      if savekeys is not None and key in savekeys: outfile.write(lines[0])

  # Clean up
  if savekeys is not None: outfile.close


def getChainData(filename, hdf5_assignments=None, labels=None, silent=False, probe_only=False):
  # Open a chain file and read it into memory

  # Regular ASCII chain
  if filename.split(":")[0][-5:] != '.hdf5':
    #Try to open chain file
    chainfile = safe_open(filename)

    #Read each line into a new list, and save that within a bigger list
    data=[]
    for line in chainfile:
      lines = line.split()
      if lines != [] and line[0] not in ['#',';','!',]:
        if not any(x == "none" for x in lines): # For now we just exclude all points with any invalid entries. 
          data.append(lines)

    #Close the chainfile and indicate success
    chainfile.close
    if not silent:
      print
      print '  Read chain '+filename
      print

    #Turn the whole lot into a numpy array of doubles
    return np.array(data, dtype=np.float64)

  # HDF5 file
  else:
    filename, groupname = filename.split(":")
    if not silent:
      print
      print "  Reading HDF5 chain file"
      print "    filename:", filename
      print "    group:", groupname
      print

    # Parse group entry
    groups = groupname.split('/')
    if groups[0] != "" or groupname == "":
      raise ValueError("Group identifier should start with '/'")
    else:
      groups = groups[1:]

    # Open HDF5 file
    try:
      import h5py
      f = h5py.File(filename,'r')
    except IOError:
      print "ERROR while reading hdf5 file.  Check filename and file format."
      quit()

    # Get relevant group entries and column names
    entries = f
    if groups[0] != "":
      for key in groups:
        try:
          entries = entries[key]
        except KeyError:
          print "ERROR: requested group \""+key+"\" does not exist in hdf5 file."
          quit()
    column_names = filter(lambda x: x[-8:] != "_isvalid", list(entries))

    if not probe_only:
      # Get data                                                                                          
      data = []                                                                                           
      data_isvalid = []                                                                                   
      for column_name in column_names:                                                                    
        try:                                                                                              
          data.append(np.array(entries[column_name], dtype=np.float64))                                   
        except AttributeError:                                                                            
          print "ERROR: \""+column_name+"\" in group \""+groupname+"\" is not convertible to a float."    
          print "Probably you gave the wrong group in your pip file."                                     
          quit()                                                                                          
        data_isvalid.append(np.array(entries[column_name+"_isvalid"], dtype=np.float64))                  
      data = np.array(data, dtype=np.float64)                                                             
      data_isvalid = np.array(data_isvalid, dtype=np.float64)
      # Print the raw number of samples in the hdf5 file                                                  
      print "  Total samples: ", data[0].size                                                             

    # Reorganize MPIrank, pointID and other requested entries for convenience.
    indices = []
    index_count = 0
    primary_column_names = ['MPIrank', 'pointID']
    sorted_column_names = primary_column_names + [x for x in column_names if x not in primary_column_names]
    for column_name in sorted_column_names:
      if hdf5_assignments.value is not None:
        while index_count in hdf5_assignments.value:
          try_append(indices, column_names, hdf5_assignments.value[index_count])
          index_count += 1
      if hdf5_assignments.value is None or column_name not in hdf5_assignments.value:
        try_append(indices, column_names, column_name)
        index_count += 1
    column_names = np.array(column_names)[indices]

    # Print probed contents and split
    if probe_only:
      for i, column_name in enumerate(column_names):
        print "   ", i, ":", column_name
      print
      quit()

    data = data[indices]
    data_isvalid = data_isvalid[indices]

    # Filter out valid points, according to the likelihood only -- and only if called with labels provided
    if (labels):
      # Based on the likelihood entry only
      likelihood_index = [value for key, value in labels.value.iteritems() if key in permittedLikes]
      likelihood_index = likelihood_index[0]
      cut = (data_isvalid[likelihood_index] == 1)
      # Based on *all* entries.
      #cut = (data_isvalid.prod(axis=0) == 1)
      data = data[:,cut]
      data_isvalid = data_isvalid[:,cut]
      print "  Total valid samples: ", data[0].size
      print "  Fraction of samples deemed valid: %.4f"%(1.0*sum(cut)/len(cut))

    # Print list of contents for convenience
    if not silent:
      for i, column_name in enumerate(column_names):
        print "   ",i, ":", column_name
        print "        mean: %.2e  min: %.2e  max %.2e"%(data[i].mean(), data[i].min(), data[i].max())
        print "        Fraction of valid points where this is invalid: %.4f"%(1.0-data_isvalid[i].mean())
      print

    return np.array(data.T, dtype=np.float64)


def try_append(indices, cols, x):
  try:
    indices.append(cols.index(x))
  except:
    print "ERROR: hdf5 file does not contain a field titled \""+x+"\"."
    quit()


def has_multiplicity(labels, hdf5_cols):
  is_hdf5 = mainChain.value.split(":")[0][-5:] == '.hdf5'
  if is_hdf5 and hdf5_cols:
    return any(x in hdf5_cols.value for x in permittedMults)
  else:
    return any(x in labels.value for x in permittedMults)


def usage():
  #Print pippi usage information
  print
  print 'You must be new here (or have fat fingers).  You can use pippi to'
  print
  print '  merge two or more chains:'
  print '    pippi merge <chain1> <chain2> ... <chainx>'
  print
  print '  post-process a chain:'
  print '    pippi pare <chain> <name_of_python_module> <name_of_function_in_module>'
  print
  print '  parse a chain using options in iniFile.pip:'
  print '    pippi parse iniFile.pip'
  print
  print '  write plotting scipts for a chain using options in iniFile.pip:'
  print '    pippi script iniFile.pip'
  print
  print '  run plotting scipts for a chain using options in iniFile.pip:'
  print '    pippi plot iniFile.pip'
  print
  print '  print an hdf5 file\'s computed column indices, using options in iniFile.pip:'
  print '    pippi probe iniFile.pip'
  print
  print '  parse, script and plot in one go:'
  print '    pippi iniFile.pip'
  print


def safe_open(filename):
  #Try to open input file
  try:
    infile = open(filename,"r")
    return infile
  except IOError:
    #Crash if file does not exist
    sys.exit('\n  Sorry, file '+filename+' not found or read-protected.  Quitting...\n')


def smart_open(filename,mode):
  #Try to open file for output
  try:
    outfile = open(filename,mode)
    return outfile
  except IOError:
    #Crash if file cannot be opened the way requested
    sys.exit('\n  Sorry, file '+filename+' cannot be opened for writing.\nCheck disk space and folder permissions.  Quitting...\n')


class dataObject:
  #Class for pip file entries, containing their values, keys and expected data types

  pipFileKey = value = dataType = ''

  def __init__(self,pipEntry,dataType):
    #Save key and expected data type for new pip file field
    self.pipFileKey = pipEntry
    self.conversion = dataType

  def seek(self,string):
    #Look for field's key in string
    if self.pipFileKey in string: return True
    return False

  def convert(self,string):
    #Convert string to appropriate data format for pip file field
    string = string.strip()
    try:
      if string == '':
        self.value = None
      else:
        self.value = self.conversion(string)
    except:
      sys.exit('Error: invalid data format in field '+self.pipFileKey+'. Quitting...\n')


#Conversion functions for parsing pip file entries

def integer(x):
  x = re.sub(r"[':;,]", '', x)
  if len(x.split()) != 1: raise Exception
  return int(x)

def floater(x):
  x = re.sub(r"[':;,]", '', x)
  if len(x.split()) != 1: raise Exception
  return float(x)

def string(x):
  if len(re.findall("'", x)) != 2: raise Exception
  return re.sub(r"^.*?'|'.*?$", '', x)

def safe_string(x):
  if len(re.findall("'", x)) != 2 or ';' in x: raise Exception
  return re.sub(r"^.*?'|'.*?$", '', x)

def boolean(x):
  x = re.sub(r"[':;,]", '', x)
  if len(x.split()) != 1: raise Exception
  if x.lower().strip() in ("yes", "true", "t", "1"):
    return True
  elif x.lower().strip() in ("no", "false", "f", "0"):
    return False
  else: raise Exception

def int_list(x): return single_list(x,integer)
def float_list(x): return single_list(x,floater)
def single_list(x,singleType):
  x = re.sub(r"[':;,]", ' ', x)
  returnVal = []
  for entry in x.split():
    returnVal.append(singleType(entry))
  return returnVal

def intuple_list(x): return tuple_list(x,integer)
def floatuple_list(x): return tuple_list(x,floater)
def tuple_list(x,tupleType):
  if len(re.findall("\{", x)) != len(re.findall("\}", x)): raise Exception
  x = re.findall("\{.+?\}", x)
  for i, pair in enumerate(x):
    pair = re.sub(r"[':;,\}\{]", ' ', pair).split()
    if len(pair) < 2: raise Exception
    x[i] = [tupleType(j) for j in pair]
  return x

def string_list(x):
  returnVal = []
  if len(re.findall("'", x))%2 != 0: raise Exception
  x = re.findall("'.+?'", x)
  for i, entry in enumerate(x):
    returnVal.append(string(entry))
  return returnVal

def float_dictionary(x):
  returnVal = {}
  x = re.findall(".+?:.+?[\s,;]+?|.+?:.+?$", x)
  for i, pair in enumerate(x):
    pair = re.sub("[\s,;]+$", '', pair).split(':')
    pair = [integer(pair[0]), floater(pair[1])]
    returnVal[pair[0]] = pair[1]
  return returnVal

def floatuple_dictionary(x):
  returnVal = {}
  x = re.findall(".+?:{.+?}[\s,;]+?|.+?:{.+?}$", x)
  for i, pair in enumerate(x):
    pair = re.sub("[\s,;]+$", '', pair).split(':')
    pair = [integer(pair[0]), floatuple_list(pair[1])[0]]
    returnVal[pair[0]] = pair[1]
  return returnVal

def string_dictionary(x):
  returnVal = {}
  if len(re.findall("'", x))%2 != 0: raise Exception
  x = re.findall("(.+?:'.+?'[\s,;]*|'.+?':.+?[\s,;]*)", x)
  for i, pair in enumerate(x):
    pair = re.sub("[\s,;]+$", '', pair).split(':')
    for j, single in enumerate(pair):
      if single[0] == '\'':
        pair[j] = string(single)
      else:
        pair[j] = integer(single)
    returnVal[pair[0]] = pair[1]
    returnVal[pair[1]] = pair[0]
  return returnVal

def int_pair_string_dictionary(x):
  returnVal = {}
  if len(re.findall("'", x))%2 != 0: raise Exception
  x = re.findall(".+?:'.+?'[\s,;]*", x)
  for i, pair in enumerate(x):
    pair = re.sub("[\s,;]+$", '', pair).split(':')
    for j, single in enumerate(pair):
      if single[0] == '\'':
        pair[j] = string(single)
      else:
        pair[j] = intuple_list(single)[0]
    if pair[0][0] in returnVal:
      returnVal[pair[0][0]][pair[0][1]] = pair[1]
    else:
      returnVal[pair[0][0]] = {pair[0][1]:pair[1]}
  return returnVal

def internal(x):
  return permittedInternals[re.sub(r"[':;,]", ' ', x).lower()]

#Global constants and simple one-liners
pippiVersion = 'v1.0'

def times1(x): return x
def half(x): return x*0.5
def negative(x): return -x
def negativehalf(x): return -x*0.5
def negln(x): return -np.log(x)
permittedLikes = {'-lnlike':times1,'lnlike':negative,'-2lnlike':half,'chi2':half,'2lnlike':negativehalf,'like':negln}
refLike = '-lnlike'

permittedMults = ['mult','mult.','multiplicity','multiplic.','mtpcty','Posterior']
refMult = permittedMults[0]

permittedPriors = ['prior','priors','pri.']
refPrior = permittedPriors[0]

mcmc = 'mcmc'
multinest = 'multinest'
other = 'other'
permittedInternals = {'mcmc':mcmc, 'multinest':multinest, 'other':other}
if permittedSchemes is not None: permittedInternals.update(permittedSchemes)
allowedIntMethods = ['bilinear', 'spline']

mainChain = dataObject('main_chain',safe_string)
secChain = dataObject('comparison_chain',safe_string)
doPosterior = dataObject('do_posterior_pdf',boolean)
doProfile = dataObject('do_profile_like',boolean)
oneDplots = dataObject('oneD_plot_quantities',int_list)
twoDplots = dataObject('twoD_plot_quantities',intuple_list)
contours = dataObject('contour_levels',float_list)
keys = [mainChain,secChain,doPosterior,doProfile,oneDplots,twoDplots,contours]
