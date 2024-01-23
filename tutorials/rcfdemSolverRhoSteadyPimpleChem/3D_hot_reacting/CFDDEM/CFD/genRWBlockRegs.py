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
z0 = 0
N = 32
RegDepth = 1.65 # 1.5
#RegRad = 0.3 #0.125
# get region radious via command line parameter
RegRad=float(sys.argv[1])

r = D/2

phi0 = 0

with open('system/RWRegs', 'w') as f:
    for i in range(0,N):
        phi = i*2*math.pi/N + phi0

        x1 = -RegRad
        y1 = r-RegDepth
        z1 = z0-RegRad

        xorigin = x1 * math.cos(phi) - y1 * math.sin(phi)
        yorigin = x1 * math.sin(phi) + y1 * math.cos(phi) 
        zorigin = z1

        x1 = 2*RegRad #RegRad
        y1 = 0  #r-RegDepth
        z1 = 0 #z0-RegRad

        xi = x1 * math.cos(phi) - y1 * math.sin(phi)
        yi = x1 * math.sin(phi) + y1 * math.cos(phi) 
        zi = z1

        x1 = 0 #-RegRad
        y1 = RegDepth + 0.05 #r+0.05
        z1 = 0 #z0-RegRad

        xj = x1 * math.cos(phi) - y1 * math.sin(phi)
        yj = x1 * math.sin(phi) + y1 * math.cos(phi)
        zj = z1

        x1 = 0 #-RegRad
        y1 = 0 #r-RegDepth
        z1 = 2*RegRad #z0+RegRad

        xk = x1 * math.cos(phi) - y1 * math.sin(phi)
        yk = x1 * math.sin(phi) + y1 * math.cos(phi)
        zk = z1




        f.write('\t{\n\t\tname RWRegs;\n\t\ttype cellSet;\n')
        if (i==0):
            f.write('\t\taction new;\n')
        else:
            f.write('\t\taction add;\n')
        f.write('\t\tsource rotatedBoxToCell;\n')
        f.write('\t\tsourceInfo\n\t\t{\n')
        f.write('\t\t\torigin (%f %f %f);\n' % (xorigin,yorigin,zorigin))
        f.write('\t\t\ti (%f %f %f);\n' % (xi,yi,zi))
        f.write('\t\t\tj (%f %f %f);\n' % (xj,yj,zj))
        f.write('\t\t\tk (%f %f %f);\n' % (xk,yk,zk))
        f.write('\t\t}\n\t}\n\n')

