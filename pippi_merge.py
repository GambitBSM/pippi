
#############################################################
# pippi: parse it, plot it
# ------------------------
# Merge program for pippi.
#
# Author: Pat Scott (patscott@physics.mcgill.ca)
# Originally developed: March 2012
#############################################################

from pippi_utils import *

def merge(filenames):
  #Function to merge multiple chains and spit them out to stdout (for 
  #e.g. pipping into a file, or some other program).  Basically
  #just 'cat' with a little bit of error checking.

  if (len(filenames) == 0):
    usage()
    return

  firstLine = True

  for filename in filenames:

    #Try to open input file
    infile = safe_open(filename)

    #Read the first line
    line = infile.readline()

    #Skip comments at the beginning of the file
    while (line[0] == '#'):
      line = infile.readline()    
    
    #Work out the number of columns in the first valid line 
    columns = len(line.split())

    while (columns != 0):
      try:
        if (firstLine):
          #Save the number of columns if this is the first line of the first chain
          columnsPredicted = columns
          firstLine = False
        elif (columns != columnsPredicted):
          #Crash if a later chain or line has a different number of columns to the first one
          sys.exit('Error: chains do not match (number of columns differ).  Quitting...')      
        #Output the current line to stdout and get the next one
        print line.rstrip('\n')
        #Read the next line
        line = infile.readline()
        #Work out the number of columns in the next line 
        columns = len(line.split())
      except IOError:
        break

    #Shut the chain file and move on to the next
    infile.close  

