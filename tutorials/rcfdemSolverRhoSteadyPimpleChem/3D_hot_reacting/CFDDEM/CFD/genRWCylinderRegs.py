#!/usr/bin/env python

import os
import shutil
import sys
import math

from operator import add


# ----------------------------------------------
# main dimensions
# ----------------------------------------------

# diameter of BF, height of RWs, number of RWs
D = 12
h = 0
N = 32
RegDepth = 1.65 # 1.5
#RegRad = 0.3 #0.125
# get region radious via command line parameter
RegRad=float(sys.argv[1])

r = D/2

phi0 = 0 #2*math.pi/N*0.5

with open('system/RWRegs', 'w') as f:
    for i in range(0,N):
        phi = i*2*math.pi/N + phi0
        x1 = (r+0.5) * math.sin(phi)
        y1 = (r+0.5) * math.cos(phi)
        z1 = h

        x2 = (r-RegDepth) * math.sin(phi)
        y2 = (r-RegDepth) * math.cos(phi)
        z2 = h

        f.write('\t{\n\t\tname RWRegs;\n\t\ttype cellSet;\n')
        if (i==0):
            f.write('\t\taction new;\n')
        else:
            f.write('\t\taction add;\n')
        f.write('\t\tsource cylinderToCell;\n')
        f.write('\t\tsourceInfo\n\t\t{\n')
#        f.write('\t\t\tp1 (%f %f %f);\n' % (x1,y1,z1))
        f.write('\t\t\tp2 (%f %f %f);\n' % (x2,y2,z2))
        f.write('\t\t\tp1 (%f %f %f);\n' % (x1,y1,z1))
        f.write('\t\t\tradius %f;\n\t\t}\n\t}\n\n' % RegRad)

