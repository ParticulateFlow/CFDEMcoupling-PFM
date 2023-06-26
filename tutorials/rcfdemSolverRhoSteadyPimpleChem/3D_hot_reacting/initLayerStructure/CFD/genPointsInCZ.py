#!/usr/bin/env python

import os
import shutil
import sys
import math

from operator import add


# ----------------------------------------------
# read set of all points and keep only those in
# predefined region
# ----------------------------------------------

def upperEdge(r):
    zu = 11
    zl = 5
    r0 = 3
    z = (zu-zl)*math.exp(-(r/r0)**2)+zl
    return z

def lowerEdge(r):
    z = upperEdge(r) - 2
    return z

with open("cellCenters") as f:
    content = f.readlines()
content = [x.rstrip(')\n') for x in content] 
content = [x.lstrip(' (') for x in content]

inPointcounter = 0
belowPointcounter = 0
fIn = open("system/pointsInCZ", "w")
fBelow = open("system/pointsBelowCZ", "w")

for i in range(0,len(content)):
    x,y,z=content[i].split()
    x=float(x)
    y=float(y)
    z=float(z)
    r=math.sqrt(x**2+y**2)
    if (z<upperEdge(r) and z>lowerEdge(r)):
        if (inPointcounter > 0):
            fIn.write('\n')
        fIn.write('\t\t(%f %f %f)' % (x, y, z))
        inPointcounter = inPointcounter + 1
    elif (z<=lowerEdge(r)):
        if (belowPointcounter > 0):
            fBelow.write('\n')
        fBelow.write('\t\t(%f %f %f)' % (x, y, z))
        belowPointcounter = belowPointcounter + 1

fIn.close()
fBelow.close()
