
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
import copy

permittedSchemes = {}

def Blockshading(colour,line_code, fill_code):
  scheme = colourScheme('Blockshading_'+colour)
  scheme.backgroundColour = '#fff'
  scheme.baseProfColourMap = '#fff--#fff(contour1)--#'+fill_code+'(contour1)--#'+fill_code
  scheme.basePostColourMap = '#fff--#fff(contour1)--#'+fill_code+'(contour1)--#'+fill_code
  scheme.mainPostContourColour2D = '\'#'+line_code+'\''
  scheme.mainProfContourColour2D = '\'#'+line_code+'\''
  scheme.mainBestFitColour1D = '#'+line_code
  scheme.mainPostMeanColour1D = '#'+line_code
  scheme.mainBestFitColour2D = '#'+line_code
  scheme.mainPostMeanColour2D = '#'+line_code
  scheme.fillTransparency2D = '0.85' #so far this is not actually functional in ctioga2; presumably it will work in later versions.
  return scheme


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
  backgroundColour = '#fff'
  baseProfColourMap = '#fff--#fff(contour2)--#f45(contour1)--#612'
  basePostColourMap = '#fff--#fff(contour2)--#88f(contour1)--#229'

  # Default values for 1D plotting styles
  fillTransparency1D = '0.85'
  main1DLineStyle = 'Solid'
  comparison1DLineStyle = 'Solid' #alt: 'Dots', 'Dashes'
  lineWidth1D = '0.9'

  # Default values for 2D contour plotting styles
  fillTransparency2D = '1.0'
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
  referenceMarkerInner = 'Cross'
  referenceMarkerInnerScale = 0.7
  referenceMarkerInnerColour = 'Black'
  referenceMarkerOuter = 'Times'
  referenceMarkerOuterScale = 0.7
  referenceMarkerOuterColour = 'Yellow'

  mainBestFitMarker = 'Star'
  mainBestFitMarkerScale = 0.8
  mainBestFitColour1D = '#300'
  mainBestFitColour2D = '#300'
  mainBestFitColourOutline2D = 'Black'

  mainPostMeanMarker = 'Bullet'
  mainPostMeanMarkerScale = 0.6
  mainPostMeanColour1D = '#004'
  mainPostMeanColour2D = '#004'
  mainPostMeanColourOutline2D = 'Black'

  comparisonBestFitMarker = 'Star'
  comparisonBestFitMarkerScale = 0.8
  comparisonBestFitColour = 'Grey'

  comparisonPostMeanMarker = 'Bullet'
  comparisonPostMeanMarkerScale = 0.6
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
    elif kind == 'obs':
      localColourMap = self.baseObsColourMap
    else:
      sys.exit('    Error: unrecognised type of colourmap requested.\n    Quitting...\n')
    for i, contour in enumerate(contours):
      localColourMap = re.sub(r'contour'+str(i+1), contour, localColourMap)
    return localColourMap

# basic colour scheme
basic = colourScheme('basic')

# iceCube colour scheme
iceCube = colourScheme('iceCube')
iceCube.backgroundColour = '#fff'
iceCube.baseProfColourMap = '#fff--#fff(contour2)--#292(contour1)--#f55(contour1)--#000'
iceCube.basePostColourMap = '#fff--#fff(contour2)--#29d(contour1)--#f55(contour1)--#000'
iceCube.baseObsColourMap = 'hls:White(contour1)--Red(contour2)--Green(contour3)'
iceCube.mainBestFitColour1D = 'Black'
iceCube.mainPostMeanColour1D = 'Black'
iceCube.mainBestFitColour2D = 'Black'
iceCube.mainPostMeanColour2D = 'Black'

# iceCube79 colour scheme
iceCube79 = colourScheme('iceCube79')
iceCube79.backgroundColour = '#fff'
iceCube79.baseProfColourMap = '#fff--#fff(contour2)--#fab(contour1)--#f45'
iceCube79.basePostColourMap = '#fff--#fff(contour2)--#ddf(contour1)--#88f'
iceCube79.baseObsColourMap = 'hls:White(contour1)--Red(contour2)--Green(contour3)'
iceCube79.mainBestFitColour1D = 'Black'
iceCube79.mainPostMeanColour1D = 'Black'
iceCube79.mainBestFitColour2D = 'Black'
iceCube79.mainPostMeanColour2D = 'Black'
iceCube79.mainPostContourColour2D = 'Grey'
iceCube79.mainProfContourColour2D = 'Grey'
iceCube79.lineWidth2D = '1.5'

# iceCube3sig colour scheme
iceCube3sig = colourScheme('iceCube3sig')
iceCube3sig.backgroundColour = '#fff'
iceCube3sig.baseProfColourMap = '#fff--#fff(contour3)--#292(contour2)--#fff(contour2)--#929(contour1)--#f55(contour1)--#000'
iceCube3sig.basePostColourMap = '#fff--#fff(contour3)--#29d(contour2)--#fff(contour2)--#929(contour1)--#f55(contour1)--#000'
iceCube3sig.baseObsColourMap = 'hls:White(contour1)--Red(contour2)--Green(contour3)'
iceCube3sig.mainBestFitColour1D = 'Black'
iceCube3sig.mainPostMeanColour1D = 'Black'
iceCube3sig.mainBestFitColour2D = 'Black'
iceCube3sig.mainPostMeanColour2D = 'Black'

# SBClassic colour scheme
SBClassic = colourScheme('SBClassic')
SBClassic.backgroundColour = '#fff'
SBClassic.baseProfColourMap = '#fff--#fff(contour2)--#2f2(contour1)--#f33(0.5)--#000'
SBClassic.basePostColourMap = '#fff--#fff(contour2)--#95d(contour1)--#f33(0.5)--#000'
SBClassic.baseObsColourMap = 'hls:White(contour1)--Red(contour2)--Green(contour3)'
SBClassic.mainBestFitColour1D = 'Black'
SBClassic.mainPostMeanColour1D = 'Black'
SBClassic.mainBestFitColour2D = 'Black'
SBClassic.mainPostMeanColour2D = 'Black'

# BlueGold colour scheme
BlueGold = colourScheme('BlueGold')
BlueGold.backgroundColour = '#fff'
BlueGold.baseProfColourMap = '#fff--#fff(contour2)--#f44(contour2)--#f44(contour1)--#ece(contour1)--#ece'
BlueGold.basePostColourMap = '#fff--#fff(contour2)--#44f(contour2)--#44f(contour1)--#fc0(contour1)--#fc0'
BlueGold.baseObsColourMap = 'hls:White(contour1)--Red(contour2)--Green(contour3)'
BlueGold.mainPostContourColour2D = 'DarkBlue'
BlueGold.mainProfContourColour2D = 'Maroon'
BlueGold.mainBestFitColour = 'Black'
BlueGold.mainBestFitColour1D = 'Black'
BlueGold.mainPostMeanColour1D = 'Black'
BlueGold.mainBestFitColour2D = 'Black'
BlueGold.mainPostMeanColour2D = 'Black'

# nightOfTheAllanachs colour scheme
nightOfTheAllanachs = colourScheme('nightOfTheAllanachs')
nightOfTheAllanachs.backgroundColour = '#000'
nightOfTheAllanachs.basePostColourMap = '#000--#000(contour2)--#808(contour1)--#f33(0.5)--#ff0'
nightOfTheAllanachs.baseProfColourMap = '#000--#000(contour2)--#33f(contour1)--#0ff(0.5)--#ff0'
nightOfTheAllanachs.baseObsColourMap = 'Black(contour1)--Red(contour2)--Green(contour3)'
nightOfTheAllanachs.mainPostContourColour2D = 'White'
nightOfTheAllanachs.mainProfContourColour2D = 'White'
nightOfTheAllanachs.axisColour2D = 'White'
nightOfTheAllanachs.mainBestFitColour1D = 'Black'
nightOfTheAllanachs.mainPostMeanColour1D = 'Black'
nightOfTheAllanachs.mainBestFitColour2D = 'White'
nightOfTheAllanachs.mainPostMeanColour2D = 'White'
nightOfTheAllanachs.legendTextColour2D = 'White'
nightOfTheAllanachs.keyTextColour2D = 'White'
nightOfTheAllanachs.comparisonContourStyle = 'Dashes'
nightOfTheAllanachs.comparison1DLineStyle = 'Dashes'

# nightOfTheAllanachs2 colour scheme
nightOfTheAllanachs2 = colourScheme('nightOfTheAllanachs2')
nightOfTheAllanachs2.backgroundColour = '#000'
nightOfTheAllanachs2.basePostColourMap = '#000--#000(contour2)--#808(contour1)--#f33(0.5)--#ff0'
nightOfTheAllanachs2.baseProfColourMap = '#000--#000(contour2)--#33f(contour1)--#0ff(0.5)--#ff0'
nightOfTheAllanachs2.baseObsColourMap = 'Black(contour1)--Red(contour2)--#00FFFF(contour3)'
nightOfTheAllanachs2.mainPostContourColour2D = 'White'
nightOfTheAllanachs2.mainProfContourColour2D = 'White'
nightOfTheAllanachs2.axisColour2D = 'White'
nightOfTheAllanachs2.mainBestFitColour1D = 'Red'
nightOfTheAllanachs2.mainPostMeanColour1D = 'Blue'
nightOfTheAllanachs2.mainBestFitColour2D = 'White'
nightOfTheAllanachs2.mainBestFitColourOutline2D = 'Black'
nightOfTheAllanachs2.mainPostMeanColour2D = 'White'
nightOfTheAllanachs2.mainPostMeanColourOutline2D = 'Black'
nightOfTheAllanachs2.legendTextColour2D = 'White'
nightOfTheAllanachs2.keyTextColour2D = 'White'
nightOfTheAllanachs2.comparisonContourStyle = 'Dashes'
nightOfTheAllanachs2.comparison1DLineStyle = 'Dashes'

# nightOfTheAllanachs3 colour scheme
nightOfTheAllanachs3 = colourScheme('nightOfTheAllanachs3')
nightOfTheAllanachs3.backgroundcolour = '#000'
nightOfTheAllanachs3.basePostColourMap = '#000--#000(contour2)--#808(contour1)--#f33(0.5)--#ff0'
nightOfTheAllanachs3.baseProfColourMap = '#000--#000(contour2)--#33f(contour1)--#0ff(0.5)--#ff0'
nightOfTheAllanachs3.baseObsColourMap = 'Black(contour1)--Blue(contour2)--Orange(contour3)'
nightOfTheAllanachs3.mainPostContourColour2D = 'White'
nightOfTheAllanachs3.mainProfContourColour2D = 'White'
nightOfTheAllanachs3.axisColour2D = 'White'
nightOfTheAllanachs3.mainBestFitColour1D = 'Red'
nightOfTheAllanachs3.mainPostMeanColour1D = 'Blue'
nightOfTheAllanachs3.mainBestFitColour2D = 'White'
nightOfTheAllanachs3.mainBestFitColourOutline2D = 'Black'
nightOfTheAllanachs3.mainPostMeanColour2D = 'White'
nightOfTheAllanachs3.mainPostMeanColourOutline2D = 'Black'
nightOfTheAllanachs3.legendTextColour2D = 'White'
nightOfTheAllanachs3.keyTextColour2D = 'White'
nightOfTheAllanachs3.comparisonContourStyle = 'Dashes'
nightOfTheAllanachs3.comparison1DLineStyle = 'Dashes'

# Blockshading colour schemes
Blockshading_red = Blockshading("red", "800", "e00")
Blockshading_green = Blockshading("green", "080", "0e0")
Blockshading_blue = Blockshading("blue", "005", "44f")
Blockshading_pink = Blockshading("pink", "808", "e0e")
Blockshading_purple = Blockshading("purple", "303", "80e")
Blockshading_orange = Blockshading("orange", "840", "f90")
Blockshading_yellow = Blockshading("yellow", "870", "fe0")
Blockshading_cyan = Blockshading("cyan", "088", "3ee")

# Extended Basic scheme - red
ExtendedBasic_red = colourScheme('extendedBasic_red')
ExtendedBasic_red.basePostColourMap = '#fff--#fff(contour3)--#fab(contour2)--#f45(contour1)--#612'
ExtendedBasic_red.baseProfColourMap = '#fff--#fff(contour3)--#fab(contour2)--#f45(contour1)--#612'
ExtendedBasic_red.baseObsColourMap = '#fff(contour1)--#00f(contour2)--#0f0(contour3)'
ExtendedBasic_red.mainBestFitColour1D = '#300'
ExtendedBasic_red.mainBestFitColour2D = '#300'
ExtendedBasic_red.mainPostContourColour2D = 'Black'
ExtendedBasic_red.mainProfContourColour2D = 'Black'
ExtendedBasic_red.axisColour2D = 'Black'
ExtendedBasic_red.mainBestFitColour1D = 'Red'
ExtendedBasic_red.mainPostMeanColour1D = 'Red'
ExtendedBasic_red.mainBestFitColour2D = 'White'
ExtendedBasic_red.mainBestFitColourOutline2D = 'Black'
ExtendedBasic_red.mainPostMeanColour2D = 'White'
ExtendedBasic_red.mainPostMeanColourOutline2D = 'Black'
ExtendedBasic_red.legendTextColour2D = 'Black'
ExtendedBasic_red.keyTextColour2D = 'Black'
ExtendedBasic_red.comparisonContourStyle = 'Dashes'
ExtendedBasic_red.comparison1DLineStyle = 'Dashes'

# Extended Basic scheme - green
ExtendedBasic_green = colourScheme('extendedBasic_green')
ExtendedBasic_green.basePostColourMap = '#fff--#fff(contour3)--#ad9(contour2)--#3a5(contour1)--#162'
ExtendedBasic_green.baseProfColourMap = '#fff--#fff(contour3)--#ad9(contour2)--#3a5(contour1)--#162'
ExtendedBasic_green.baseObsColourMap = '#fff(contour1)--#00f(contour2)--#f00(contour3)'
ExtendedBasic_green.mainBestFitColour1D = '#030'
ExtendedBasic_green.mainBestFitColour2D = '#030'
ExtendedBasic_green.mainPostContourColour2D = 'Black'
ExtendedBasic_green.mainProfContourColour2D = 'Black'
ExtendedBasic_green.axisColour2D = 'Black'
ExtendedBasic_green.mainBestFitColour1D = 'Green'
ExtendedBasic_green.mainPostMeanColour1D = 'Green'
ExtendedBasic_green.mainBestFitColour2D = 'White'
ExtendedBasic_green.mainBestFitColourOutline2D = 'Black'
ExtendedBasic_green.mainPostMeanColour2D = 'White'
ExtendedBasic_green.mainPostMeanColourOutline2D = 'Black'
ExtendedBasic_green.legendTextColour2D = 'Black'
ExtendedBasic_green.keyTextColour2D = 'Black'
ExtendedBasic_green.comparisonContourStyle = 'Dashes'
ExtendedBasic_green.comparison1DLineStyle = 'Dashes'

# Extended Basic scheme - blue
ExtendedBasic_blue = colourScheme('extendedBasic_blue')
ExtendedBasic_blue.basePostColourMap = '#fff--#fff(contour3)--#9ce(contour2)--#38b(contour1)--#126'
ExtendedBasic_blue.baseProfColourMap = '#fff--#fff(contour3)--#9ce(contour2)--#38b(contour1)--#126'
ExtendedBasic_blue.baseObsColourMap = '#fff(contour1)--#f00(contour2)--#0f0(contour3)'
ExtendedBasic_blue.mainBestFitColour1D = '#003'
ExtendedBasic_blue.mainBestFitColour2D = '#003'
ExtendedBasic_blue.mainPostContourColour2D = 'Black'
ExtendedBasic_blue.mainProfContourColour2D = 'Black'
ExtendedBasic_blue.axisColour2D = 'Black'
ExtendedBasic_blue.mainBestFitColour1D = 'Blue'
ExtendedBasic_blue.mainPostMeanColour1D = 'Blue'
ExtendedBasic_blue.mainBestFitColour2D = 'White'
ExtendedBasic_blue.mainBestFitColourOutline2D = 'Black'
ExtendedBasic_blue.mainPostMeanColour2D = 'White'
ExtendedBasic_blue.mainPostMeanColourOutline2D = 'Black'
ExtendedBasic_blue.legendTextColour2D = 'Black'
ExtendedBasic_blue.keyTextColour2D = 'Black'
ExtendedBasic_blue.comparisonContourStyle = 'Dashes'
ExtendedBasic_blue.comparison1DLineStyle = 'Dashes'

# Extended Basic scheme - grey
ExtendedBasic_grey = colourScheme('extendedBasic_grey')
ExtendedBasic_grey.basePostColourMap = '#fff--#fff(contour3)--#aaa(contour2)--#555(contour1)--#222'
ExtendedBasic_grey.baseProfColourMap = '#fff--#fff(contour3)--#aaa(contour2)--#555(contour1)--#222'
ExtendedBasic_grey.baseObsColourMap = '#fff(contour1)--#999(contour2)--#333(contour3)'
ExtendedBasic_grey.mainBestFitColour1D = '#000'
ExtendedBasic_grey.mainBestFitColour2D = '#000'
ExtendedBasic_grey.mainPostContourColour2D = 'Black'
ExtendedBasic_grey.mainProfContourColour2D = 'Black'
ExtendedBasic_grey.axisColour2D = 'Black'
ExtendedBasic_grey.mainBestFitColour1D = 'Black'
ExtendedBasic_grey.mainPostMeanColour1D = 'Black'
ExtendedBasic_grey.mainBestFitColour2D = 'White'
ExtendedBasic_grey.mainBestFitColourOutline2D = 'Black'
ExtendedBasic_grey.mainPostMeanColour2D = 'White'
ExtendedBasic_grey.mainPostMeanColourOutline2D = 'Black'
ExtendedBasic_grey.legendTextColour2D = 'Black'
ExtendedBasic_grey.keyTextColour2D = 'Black'
ExtendedBasic_grey.comparisonContourStyle = 'Dashes'
ExtendedBasic_grey.comparison1DLineStyle = 'Dashes'

# Sharp Basic scheme - red
SharpBasic_red = colourScheme('sharpBasic_red')
SharpBasic_red.basePostColourMap = '#fff--#fff(contour3)--#fab(contour2)--#f89(contour2)--#d45(contour1)--#b34(contour1)--#612'
SharpBasic_red.baseProfColourMap = '#fff--#fff(contour3)--#fab(contour2)--#f89(contour2)--#d45(contour1)--#b34(contour1)--#612'
SharpBasic_red.baseObsColourMap = '#fff(contour1)--#7a9(contour2)--#365(contour3)'
SharpBasic_red.mainBestFitColour1D = '#300'
SharpBasic_red.mainBestFitColour2D = '#300'
SharpBasic_red.mainPostContourColour2D = 'Black'
SharpBasic_red.mainProfContourColour2D = 'Black'
SharpBasic_red.axisColour2D = 'Black'
SharpBasic_red.mainBestFitColour1D = 'Red'
SharpBasic_red.mainPostMeanColour1D = 'Red'
SharpBasic_red.mainBestFitColour2D = 'White'
SharpBasic_red.mainBestFitColourOutline2D = 'Black'
SharpBasic_red.mainPostMeanColour2D = 'White'
SharpBasic_red.mainPostMeanColourOutline2D = 'Black'
SharpBasic_red.legendTextColour2D = 'Black'
SharpBasic_red.keyTextColour2D = 'Black'
SharpBasic_red.comparisonContourStyle = 'Dashes'
SharpBasic_red.comparison1DLineStyle = 'Dashes'
SharpBasic_red.comparisonBestFitColour = 'Green'
SharpBasic_red.comparisonPostMeanColour = 'Green'
SharpBasic_red.comparisonPostContourColour2D = 'Green'
SharpBasic_red.comparisonProfContourColour2D = 'Green'


# Sharp Basic scheme - green
SharpBasic_green = colourScheme('sharpBasic_green')
SharpBasic_green.basePostColourMap = '#fff--#fff(contour3)--#ad9(contour2)--#7a9(contour2)--#385(contour1)--#365(contour1)--#142'
SharpBasic_green.baseProfColourMap = '#fff--#fff(contour3)--#ad9(contour2)--#7a9(contour2)--#385(contour1)--#365(contour1)--#142'
SharpBasic_green.baseObsColourMap = '#fff(contour1)--#f89(contour2)--#b34(contour3)'
SharpBasic_green.mainBestFitColour1D = '#030'
SharpBasic_green.mainBestFitColour2D = '#030'
SharpBasic_green.mainPostContourColour2D = 'Black'
SharpBasic_green.mainProfContourColour2D = 'Black'
SharpBasic_green.axisColour2D = 'Black'
SharpBasic_green.mainBestFitColour1D = 'Green'
SharpBasic_green.mainPostMeanColour1D = 'Green'
SharpBasic_green.mainBestFitColour2D = 'White'
SharpBasic_green.mainBestFitColourOutline2D = 'Black'
SharpBasic_green.mainPostMeanColour2D = 'White'
SharpBasic_green.mainPostMeanColourOutline2D = 'Black'
SharpBasic_green.legendTextColour2D = 'Black'
SharpBasic_green.keyTextColour2D = 'Black'
SharpBasic_green.comparisonContourStyle = 'Dashes'
SharpBasic_green.comparison1DLineStyle = 'Dashes'
SharpBasic_green.comparisonBestFitColour = 'Red'
SharpBasic_green.comparisonPostMeanColour = 'Red'
SharpBasic_green.comparisonPostContourColour2D = 'Red'
SharpBasic_green.comparisonProfContourColour2D = 'Red'

# Sharp Basic scheme - blue
SharpBasic_blue = colourScheme('sharpBasic_blue')
SharpBasic_blue.basePostColourMap = '#fff--#fff(contour3)--#9ce(contour2)--#6be(contour2)--#38b(contour1)--#169(contour1)--#126'
SharpBasic_blue.baseProfColourMap = '#fff--#fff(contour3)--#9ce(contour2)--#6be(contour2)--#38b(contour1)--#169(contour1)--#126'
SharpBasic_blue.baseObsColourMap = '#fff(contour1)--#f92(contour2)--#930(contour3)'
SharpBasic_blue.mainBestFitColour1D = '#003'
SharpBasic_blue.mainBestFitColour2D = '#003'
SharpBasic_blue.mainBestFitColour2D = 'White'
SharpBasic_blue.mainPostContourColour2D = 'Black'
SharpBasic_blue.mainProfContourColour2D = 'Black'
SharpBasic_blue.axisColour2D = 'Black'
SharpBasic_blue.mainBestFitColour1D = 'Blue'
SharpBasic_blue.mainPostMeanColour1D = 'Blue'
SharpBasic_blue.mainBestFitColour2D = 'White'
SharpBasic_blue.mainBestFitColourOutline2D = 'Black'
SharpBasic_blue.mainPostMeanColour2D = 'White'
SharpBasic_blue.mainPostMeanColourOutline2D = 'Black'
SharpBasic_blue.legendTextColour2D = 'Black'
SharpBasic_blue.keyTextColour2D = 'Black'
SharpBasic_blue.comparisonContourStyle = 'Solid'
SharpBasic_blue.comparison1DLineStyle = 'Solid'
#SharpBasic_blue.comparisonBestFitColour2D = 'Chocolate'
#SharpBasic_blue.comparisonBestFitcolourOutline2D = 'DarkChocolate'
#SharpBasic_blue.comparisonPostMeanColour2D = 'Chocolate'
#SharpBasic_blue.comparisonPostMeanColourOutline2D = 'DarkChocolate'
#SharpBasic_blue.comparisonPostContourColour2D = 'DarkChocolate'
#SharpBasic_blue.comparisonProfContourColour2D = 'DarkChocolate'

# Sharp Basic scheme - grey
SharpBasic_grey = colourScheme('sharpBasic_grey')
SharpBasic_grey.basePostColourMap = '#fff--#fff(contour3)--#ccc(contour2)--#aaa(contour2)--#777(contour1)--#555(contour1)--#222'
SharpBasic_grey.baseProfColourMap = '#fff--#fff(contour3)--#ccc(contour2)--#aaa(contour2)--#777(contour1)--#555(contour1)--#222'
SharpBasic_grey.baseObsColourMap = '#fff--#fff(contour1)--#ccc(contour2)--#aaa(contour2)--#777(contour3)--#555(contour3)--#222'
SharpBasic_grey.mainBestFitColour1D = '#000'
SharpBasic_grey.mainBestFitColour2D = '#000'
SharpBasic_grey.mainPostContourColour2D = 'Black'
SharpBasic_grey.mainProfContourColour2D = 'Black'
SharpBasic_grey.axisColour2D = 'Black'
SharpBasic_grey.mainBestFitColour1D = 'Black'
SharpBasic_grey.mainPostMeanColour1D = 'Black'
SharpBasic_grey.mainBestFitColour2D = 'White'
SharpBasic_grey.mainBestFitColourOutline2D = 'Black'
SharpBasic_grey.mainPostMeanColour2D = 'White'
SharpBasic_grey.mainPostMeanColourOutline2D = 'Black'
SharpBasic_grey.legendTextColour2D = 'Black'
SharpBasic_grey.keyTextColour2D = 'Black'
SharpBasic_grey.comparisonContourStyle = 'Dashes'
SharpBasic_grey.comparison1DLineStyle = 'Dashes'

# Orange-Red scheme
OrangeRed = colourScheme('orangeRed')
OrangeRed.baseProfColourMap = '#fff--#fff(contour3)--#dcd0b7(contour2)--#fdcc8a(contour2)--#da6b37(contour1)--#f9523f(contour1)--#612'
OrangeRed.mainBestFitColour1D = '#300'
OrangeRed.mainBestFitColour2D = '#300'

# Blue-Purple scheme
BluePurple = colourScheme('bluePurple')
BluePurple.baseProfColourMap = '#fff--#fff(contour3)--#cbd6d9(contour2)--#b3cde3(contour2)--#6a74a4(contour1)--#aa63bf(contour1)--#424'
BluePurple.mainBestFitColour1D = '#303'
BluePurple.mainBestFitColour2D = '#303'

# Green-Blue scheme
GreenBlue = colourScheme('greenBlue')
GreenBlue.baseProfColourMap = '#fff--#fff(contour3)--#d0d7c6(contour2)--#bae4bc(contour2)--#59aaa2(contour1)--#4daedf(contour1)--#126'
GreenBlue.mainBestFitColour1D = '#003'
GreenBlue.mainBestFitColour2D = '#003'

# Yellow-Green scheme
YellowGreen = colourScheme('yellowGreen')
YellowGreen.baseProfColourMap = '#fff--#fff(contour3)--#ddddaa(contour2)--#c2e699(contour2)--#56a457(contour1)--#45a665(contour1)--#142'
YellowGreen.mainBestFitColour1D = '#030'
YellowGreen.mainBestFitColour2D = '#030'
