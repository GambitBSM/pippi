
#############################################################
# pippi: parse it, plot it
# ------------------------
# Colour scheme module for pippi.  Add your own at the end of
# this file.
#
# Author: Pat Scott (patscott@physics.mcgill.ca)
# Originally developed: March 2012
#############################################################

import re

permittedSchemes = {}

class colourScheme:
  # Class for pippi plotting colour schemes

  # Name of colour scheme
  name = ''

  # Default values for colours
  mainPostColour1D = 'Blue'
  mainProfColour1D = 'Red'
  mainPostContourColour2D = 'Black'
  mainProfContourColour2D = 'Black'
  comparisonPostColour1D = 'Grey'
  comparisonProfColour1D = 'Grey'
  comparisonPostContourColour2D = 'Grey'
  comparisonProfContourColour2D = 'Grey'
  baseProfColourMap = '#fff--#fff(contour2)--#f45(contour1)--#612'
  basePostColourMap = '#fff--#fff(contour2)--#88f(contour1)--#229'

  # Default values for 1D plotting styles
  fillTransparency1D = '0.85'
  main1DLineStyle = 'Solid'
  comparison1DLineStyle = 'Solid' #alt: 'Dots', 'Dashes'
  lineWidth1D = '0.9'

  # Default values for 2D contour plotting styles
  mainContourStyle = 'Solid'
  comparisonContourStyle = 'Solid'
  lineWidth2D = '0.9'

  # Default text and axis colours
  legendTextColour1D = 'Black'
  keyTextColour1D = 'Black'
  axisColour1D = 'Black'
  legendTextColour2D = 'Black'
  keyTextColour2D = 'Black'
  axisColour2D = 'Black'

  # Default markers and their colours
  referenceMarkerInner = 'Times'
  referenceMarkerInnerScale = 1.2
  referenceMarkerInnerColour = 'Yellow'
  referenceMarkerOuter = 'TimesOpen'
  referenceMarkerOuterScale = 1.2
  referenceMarkerOuterColour = 'Grey'

  mainBestFitMarker = 'Star'
  mainBestFitMarkerScale = 1.5
  mainBestFitColour1D = '#300'
  mainBestFitColour2D = '#300'

  mainPostMeanMarker = 'Bullet'
  mainPostMeanMarkerScale = 0.9
  mainPostMeanColour1D = '#004'
  mainPostMeanColour2D = '#004'

  comparisonBestFitMarker = 'Star'
  comparisonBestFitMarkerScale = 1.5
  comparisonBestFitColour = 'Grey'

  comparisonPostMeanMarker = 'Bullet'
  comparisonPostMeanMarkerScale = 0.9
  comparisonPostMeanColour = 'Grey'

  def __init__(self,name):
    global permittedSchemes
    name = name.lower()
    self.name = name
    if permittedSchemes is None:
      permittedSchemes = {name:self}
    else:
      permittedSchemes[name] = self

  def colourMap(self,contours,kind):
    #Construct colourmap from base colour map and contour levels
    if kind == 'post':
      localColourMap = self.basePostColourMap
    elif kind == 'like':
      localColourMap = self.baseProfColourMap
    else:
      sys.exit('    Error: unrecognised type of colourmap requested.\n    Quitting...\n')
    for i, contour in enumerate(contours):
      localColourMap = re.sub(r'contour'+str(i+1), contour, localColourMap) 
    return localColourMap
    
# basic colour scheme
basic = colourScheme('basic')

# iceCube colour scheme
iceCube = colourScheme('iceCube')
iceCube.baseProfColourMap = '#fff--#fff(contour2)--#292(contour1)--#f55(contour1)--#000'
iceCube.basePostColourMap = '#fff--#fff(contour2)--#29d(contour1)--#f55(contour1)--#000'
iceCube.mainBestFitColour1D = 'Black'
iceCube.mainPostMeanColour1D = 'Black'
iceCube.mainBestFitColour2D = 'Black'
iceCube.mainPostMeanColour2D = 'Black'

# iceCube3sig colour scheme
iceCube = colourScheme('iceCube3sig')
iceCube.baseProfColourMap = '#fff--#fff(contour3)--#292(contour2)--#fff(contour2)--#929(contour1)--#f55(contour1)--#000'
iceCube.basePostColourMap = '#fff--#fff(contour3)--#29d(contour2)--#fff(contour2)--#929(contour1)--#f55(contour1)--#000'
iceCube.mainBestFitColour1D = 'Black'
iceCube.mainPostMeanColour1D = 'Black'
iceCube.mainBestFitColour2D = 'Black'
iceCube.mainPostMeanColour2D = 'Black'

# SBClassic colour scheme
SBClassic = colourScheme('SBClassic')
SBClassic.baseProfColourMap = '#fff--#fff(contour2)--#2f2(contour1)--#f33(0.5)--#000'
SBClassic.basePostColourMap = '#fff--#fff(contour2)--#95d(contour1)--#f33(0.5)--#000'
SBClassic.mainBestFitColour1D = 'Black'
SBClassic.mainPostMeanColour1D = 'Black'
SBClassic.mainBestFitColour2D = 'Black'
SBClassic.mainPostMeanColour2D = 'Black'

# BlueGold colour scheme
BlueGold = colourScheme('BlueGold')
BlueGold.baseProfColourMap = '#fff--#fff(contour2)--#f44(contour2)--#f44(contour1)--#ece(contour1)--#ece'
BlueGold.basePostColourMap = '#fff--#fff(contour2)--#44f(contour2)--#44f(contour1)--#fc0(contour1)--#fc0'
BlueGold.mainPostContourColour2D = 'DarkBlue'
BlueGold.mainProfContourColour2D = 'Maroon'
BlueGold.mainBestFitColour = 'Black'
BlueGold.mainBestFitColour1D = 'Black'
BlueGold.mainPostMeanColour1D = 'Black'
BlueGold.mainBestFitColour2D = 'Black'
BlueGold.mainPostMeanColour2D = 'Black'

# nightOfTheAllanachs colour scheme
nightOfTheAllanachs = colourScheme('nightOfTheAllanachs')
nightOfTheAllanachs.basePostColourMap = '#000--#000(contour2)--#808(contour1)--#f33(0.5)--#ff0'
nightOfTheAllanachs.baseProfColourMap = '#000--#000(contour2)--#33f(contour1)--#0ff(0.5)--#ff0'
nightOfTheAllanachs.mainPostContourColour2D = 'White'
nightOfTheAllanachs.mainProfContourColour2D = 'White'
nightOfTheAllanachs.axisColour2D = 'White'
nightOfTheAllanachs.mainBestFitColour1D = 'Black'
nightOfTheAllanachs.mainPostMeanColour1D = 'Black'
nightOfTheAllanachs.mainBestFitColour2D = 'White'
nightOfTheAllanachs.mainPostMeanColour2D = 'White'
nightOfTheAllanachs.legendTextColour2D = 'White'
nightOfTheAllanachs.keyTextColour2D = 'White'


