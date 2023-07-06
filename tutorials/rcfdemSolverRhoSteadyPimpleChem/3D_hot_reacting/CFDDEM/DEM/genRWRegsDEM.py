#!/usr/bin/env python

import os
import shutil
import sys
import math

from operator import add

from random import seed
from random import randint
# seed random number generator
seed(1)


# ----------------------------------------------
# main dimensions
# ----------------------------------------------

# diameter of BF, height of RWs, number of RWs
r = 6    # wall radius
xcellorig = 0.37300
d1 = 0.5*xcellorig*math.sqrt(2)   # inner diameter of cylinder part; cells have been refined to have the orig side length
d2 = 0.5   # inner diameter of spherical cavity
D1 = d1*1.1 # outer diameter of cylinder part
D2 = d2*1.4 # outer diameter of spherical cavity
depth = 1.5 # distance of sphere center from wall


z0 = 0 # reference height, axis of cylinder
N = 32 # number of raceways
phi0 = 0 # reference angle

with open('in.RWRegions', 'w') as f:
    for i in range(0,N):
        phi = i*2*math.pi/N + phi0
        x1 = (r - depth) * math.sin(phi)
        y1 = (r - depth) * math.cos(phi)
        z1 = z0 - 0.5*d1 + 0.5*d2

        f.write('region RWSphere%d sphere %f %f %f %f\n' % (i,x1,y1,z1,0.5*d2))
        f.write('region RWSphereOuter%d sphere %f %f %f %f\n' % (i,x1,y1,z1,0.5*D2))
    f.write('\n\n')


    for i in range(0,N):
        phi = i*2*math.pi/N + phi0

        f.write('variable phi%d equal %f\n' % (i, phi))
        f.write('region RWCylinder%d cylinder x 0.0 %f %f %f %f ' % (i, z0, 0.5*d1, r-depth ,r))
        f.write('rotate v_phi%d 0 0 0 0 0 1\n' % i)
        f.write('region RWCylinderOuter%d cylinder x 0.0 %f %f %f %f ' % (i, z0, 0.5*D1, r-depth ,r))
        f.write('rotate v_phi%d 0 0 0 0 0 1\n' % i)
    f.write('\n\n')

# combine outer raceways into single region
    f.write('region RWOuter union %d ' % (2*N))
    for i in range(0,N):
        f.write('RWSphereOuter%d RWCylinderOuter%d ' % (i,i))
    f.write('\n\n')




# leave inner raceways separate for equal removal rates
#    for i in range(0,N):
#        f.write('region RWInner%d union 2 RWSphereInner%d RWCylinderInner%d\n' % (i,i,i))

#    f.write('\n\n')
#    for i in range(0,N):
#        f.write('delete_atoms region RWInner%d\n' % i)

#    for i in range(0,N):
#        f.write('delete_atoms region RWInner%d\n' % i)


# combine inner raceways into single region
    f.write('region RW union %d ' % (2*N))
    for i in range(0,N):
        f.write('RWSphere%d RWCylinder%d ' % (i,i))
    f.write('\n\n')
    f.write('delete_atoms region RW\n')
    f.write('group RW dynamic all region RW every 100')


with open('in.RWRemoval', 'w') as f:
    f.write('#########\n# removal\n#########\n\n')
    f.write('# remove material from inner RW regions to be mainly empty and additional material from corona around it\n\n')
    f.write('variable rCMassTot equal 32*${rCMass}\n')
    f.write('fix delRW all remove nevery ${NRW} massrate ${rCMassTot} style delete seed 1234567 region RW verbose no restart_read no\n')
    f.write('fix delRWOuter all remove nevery ${NRW} massrate ${rCMassTot} style delete seed 1234568 region RWOuter verbose no restart_read no')
#    for i in range(0,N):
#        f.write('fix delRW%d all remove nevery ${NRW} massrate ${rCMass} style delete seed %d region RWInner%d verbose no restart_read no\n' % (i,randint(0, 10000000),i))
 
