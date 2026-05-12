# Python script to create plasticStrainVsYieldStress file for solids4foam
# plasticity laws
# Author: Andrew Whelan, UCD.
# Updated by Philip Cardiff, UCD.

import sys
import math

def hardeningLaw(x):
    y=320e6+688e6*x
    #y=(3**0.5)*505.92e6*(0.0038+x)**0.12
    return '('+str(x)+'    '+str(y)+')'

strainValues=[0.0001*(i) for i in range(10000)]

mainStr='( \n'
for i in strainValues:
    str1=hardeningLaw(i)
    mainStr=mainStr+str1+'\n'

mainStr=mainStr+')'
text_file = open("plasticStrainVsYieldStress", "w")
text_file.write(mainStr)
text_file.close()

