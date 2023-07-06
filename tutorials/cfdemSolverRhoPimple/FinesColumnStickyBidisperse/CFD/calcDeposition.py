import os
import shutil
import math

with open('postProcessing/volInt/0/volFieldValue.dat', 'r') as f:
    for line in f:
        if line.startswith('#'): 
            continue
        row = line.split()
    alphaTot = float(row[1])

with open('postProcessing/volIntLower/0/volFieldValue.dat', 'r') as f:
    for line in f:
        if line.startswith('#'): 
            continue
        row = line.split()
    alphaLower = float(row[1])

with open('postProcessing/volIntUpper/0/volFieldValue.dat', 'r') as f:
    for line in f:
        if line.startswith('#'): 
            continue
        row = line.split()
    alphaUpper = float(row[1])

with open('inletVelocity','r') as f:
    for line in f:
        row = line.split()
        v = float(row[1].rstrip(";"))

v_sf = abs(v)/1.25

vol=0.05**2*math.pi*0.3

alphaTot *= 1.0/vol
alphaLower *= 2.0/vol
alphaUpper *= 2.0/vol

with open("deposition.txt","w+") as f:
    f.write("# v_sf\talpha_tot\talpha_lower\talpha_upper\n")
    f.write("#--------------------------------------------------\n")
    f.write(str(v_sf)+"\t"+str(alphaTot)+"\t"+str(alphaLower)+"\t"+str(alphaUpper))
