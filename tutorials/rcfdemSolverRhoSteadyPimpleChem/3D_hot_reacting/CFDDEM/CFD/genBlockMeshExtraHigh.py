#!/usr/bin/env python

import os
import shutil
import sys
import math

from operator import add

exec(open('libMagicMesh.py').read())


# ----------------------------------------------
# main dimensions
# ----------------------------------------------

# reference diameter
D = 12

# relative diameter, absolute height coordinates; lowest height such that RW height at layer between cells
dVec = [1.0, 1.0, 1.23, 1.23, 0.7825, 0.7825, 0.7825]
hVec = [-1.12, 0.746, 3.913, 6.408, 24.035, 26.035, 28.0]
if (dVec.__len__() != hVec.__len__()):
    print('Inconsistent numbers of heights and diameters.')
    sys.exit()


# ----------------------------------------------
# discretisation
# ----------------------------------------------

dt = 0.35

nTan = 12
dz = 0.35  # square at RW level: dz = 12*math.pi/(8*nTan)



# ----------------------------------------------
# no user intervention beyond this point
# ----------------------------------------------

# magic screws
fc = 0.65;
mfc = 1 - ( 1 - (math.sin(math.pi/3)/math.sin(5*math.pi/12))) / (1-1/math.sqrt(2))

# magic constants
NN = 8





heights = list()
for j in range(0, hVec.__len__()):
	heights.append(hVec[j])
# end for

ro = list()
for j in range(0, dVec.__len__()):
	ro.append(0.5*D*dVec[j])
# end for

hIndices = list(range(0, heights.__len__(), 1))
cHeights = heights

rc = list()
for j in range(0, ro.__len__()):
	rc.append(fc*ro[j])
# end for

rcIn = list()
for j in range(0, ro.__len__()):
	rcIn.append(rc[j]/math.sqrt(2))
# end for

rcInt = list()
for j in range(0, ro.__len__()):
	dr = (1.0-mfc)*(rc[j] - rcIn[j])
	ri = rc[j]*(1.0-1.0/math.sqrt(2.0)) + dr/math.sqrt(2.0)
	rC = ri - rc[j] + dr
	#
	ec = math.sqrt(ri*ri + rC*rC - 2.0*ri*rC*math.cos(math.pi/24.0))
	rcInt.append(ec)
# end for






# number of cells
nR = max(1, int(nTan*math.pi/(4*fc)-3))

nT = list()
for j in range(0, ro.__len__()-1):
	nT.append(nTan)
# end for

nZ = list()
for j in range(0, heights.__len__()-1):
	nZ.append( max(1, int( (heights[j+1] - heights[j])/dz) ) )
# end for


# angles
dAngle = int(360/NN)
anglesDeg = list(range(0, 360, dAngle))
angles = list(range(0, NN))
iAngles = list(range(0, NN))
indices = list(range(0, NN, 1))


offsetAngle = (360.0/NN)

for i in range(0,anglesDeg.__len__()):
	angles[i] = float(anglesDeg[i] + offsetAngle)*math.pi/180.0
	iAngles[i] = (float(anglesDeg[i] + offsetAngle)+float(dAngle)/2.0)*math.pi/180.0
# end for


topLevel = hIndices.__len__()-1


# ----------------------------------------------
# 					Automagic
# ----------------------------------------------

# some evil global variable
blockCount = 0

# point coordinates
pC = {}

# interpolation-point coordinates
piC = {}


counter = 0
pTotal = {}

# compute coordinates
for j in hIndices:
	mKey = 'mmm'+str(j)
	pC[mKey] = [0.0, 0.0, cHeights[j]]
	pTotal[mKey] = counter # don't forget to count
	counter = inc(counter)

	for i in indices:	
		cKey = 'c'+str(i)+'h'+str(j)
		if (i%2 == 0):
			radius = rc[j]
		else:
			radius = rcIn[j] + mfc*(rc[j]-rcIn[j])
		# end if
		
#		pC[cKey] = [radius*math.sin(angles[i]), cHeights[j], radius*math.cos(angles[i])]
#		piC[cKey] = [rcInt[j]*math.sin(iAngles[i]), cHeights[j], rcInt[j]*math.cos(iAngles[i])]
# cycl. permutated
		pC[cKey] = [radius*math.cos(angles[i]), radius*math.sin(angles[i]), cHeights[j]]
		piC[cKey] = [rcInt[j]*math.cos(iAngles[i]), rcInt[j]*math.sin(iAngles[i]), cHeights[j]]
		pTotal[cKey] = counter
		counter = inc(counter)
		
		oKey = 'o'+str(i)+'h'+str(j)
#		pC[oKey] = [ro[j]*math.sin(angles[i]), heights[j], ro[j]*math.cos(angles[i])]
#		piC[oKey] = [ro[j]*math.sin(iAngles[i]), heights[j], ro[j]*math.cos(iAngles[i])]
		pC[oKey] = [ro[j]*math.cos(angles[i]), ro[j]*math.sin(angles[i]), heights[j]]
		piC[oKey] = [ro[j]*math.cos(iAngles[i]), ro[j]*math.sin(iAngles[i]), heights[j]]
		pTotal[oKey] = counter
		counter = inc(counter)
	# end for all angles
# end for all heights




# ----------------------------------------------
# set up IO
# ----------------------------------------------

filename = './system/blockMeshDict'
mywfile = open(filename,'w')



# ----------------------------------------------
# write blockMeshDict
beginBlockMeshDict(mywfile)


# ----------------------------------------------
# write points
mywfile.write("vertices\n")
mywfile.write(listStart)


# write vertices a-f
for j in hIndices:
	mKey = 'mmm'+str(j)
	mywfile.write("\t"+writeCoord(str(pC[mKey]))+"\t // "+mKey+" = "+str(pTotal[mKey])+"\n")

	for i in indices:
		cKey = 'c'+str(i)+'h'+str(j)
		mywfile.write("\t"+writeCoord(str(pC[cKey]))+"\t // "+cKey+" = "+str(pTotal[cKey])+"\n")
		
		oKey = 'o'+str(i)+'h'+str(j)
		mywfile.write("\t"+writeCoord(str(pC[oKey]))+"\t // "+oKey+" = "+str(pTotal[oKey])+"\n")
	# end for
# end for


mywfile.write(listEnd)


mywfile.write(newLine)
mywfile.write(newLine)


# ----------------------------------------------
# create blocks
mywfile.write("blocks\n")
mywfile.write(listStart)

grading = [1, 1, 1]

for j in hIndices[0:hIndices.__len__()-1]:
	for i in indices:
		# centre blocks
		if (i%2 == 0):
			lowerNodes = ['mmm'+str(j), 'c'+str((i-1)%NN)+'h'+str(j), 'c'+str((i)%NN)+'h'+str(j), 'c'+str((i+1)%NN)+'h'+str(j)]
			mywfile.write("    " + writeBlock(lowerNodes, [nT[j], nT[j], nZ[j]], grading) + "\n")
		# end if
		
		# outer blocks
		lowerNodes = ['c'+str(i)+'h'+str(j), 'o'+str(i)+'h'+str(j), 'o'+str((i+1)%NN)+'h'+str(j), 'c'+str((i+1)%NN)+'h'+str(j)]
		mywfile.write("    " + writeBlock(lowerNodes, [nR, nT[j], nZ[j]], grading) + "\n")
	# end for
# end for


mywfile.write(listEnd)
mywfile.write(newLine)
mywfile.write(newLine)


# ----------------------------------------------
# write edge definitions for curved edges
mywfile.write("edges\n")
mywfile.write(listStart)

for j in hIndices:
	for i in indices:
		# centre blocks
		cKey1 = 'c'+str(i%NN)+'h'+str(j)
		cKey2 = 'c'+str((i+1)%NN)+'h'+str(j)
		mywfile.write(writeEdge(cKey1, cKey2) + "\n")
		
		# outer blocks
		oKey1 = 'o'+str(i%NN)+'h'+str(j)
		oKey2 = 'o'+str((i+1)%NN)+'h'+str(j)
		mywfile.write(writeEdge(oKey1, oKey2) + "\n")
	# end for all angles
# end for all heights


mywfile.write(listEnd)
mywfile.write(newLine)
mywfile.write(newLine)


# ----------------------------------------------
# write patches
mywfile.write("boundary\n")
mywfile.write(listStart)

beginPatch(mywfile, "bottom", "patch")

for i in indices:
	# centre blocks
	if (i%2 == 0):
		key1 = 'mmm'+str(0)
		key2 = 'c'+str((i-1)%NN)+'h'+str(0)
		key3 = 'c'+str((i)%NN)+'h'+str(0)
		key4 = 'c'+str((i+1)%NN)+'h'+str(0)
		
		mywfile.write("    " + writeFullFace([key1, key4, key3, key2]) + "\n")
	# end if
	
	# outer blocks
	key1 = 'c'+str((i-1)%NN)+'h'+str(0)
	key2 = 'o'+str((i-1)%NN)+'h'+str(0)
	
	mywfile.write("    " + writeHFace(key2, key1) + "\n")
# end for
endPatch(mywfile)


beginPatch(mywfile, "top", "patch")
# write faces
for i in indices:
	# centre blocks
	if (i%2 == 0):
		key1 = 'mmm'+str(topLevel)
		key2 = 'c'+str((i-1)%NN)+'h'+str(topLevel)
		key3 = 'c'+str((i)%NN)+'h'+str(topLevel)
		key4 = 'c'+str((i+1)%NN)+'h'+str(topLevel)
		
		mywfile.write("    " + writeFullFace([key1, key2, key3, key4]) + "\n")
	# end if
	
	# outer blocks
	key1 = 'c'+str((i-1)%NN)+'h'+str(topLevel)
	key2 = 'o'+str((i-1)%NN)+'h'+str(topLevel)
	
	mywfile.write("    " + writeHFace(key1, key2) + "\n")
# end for
endPatch(mywfile)


beginPatch(mywfile, "walls", "wall")
# write faces
for j in hIndices[0:hIndices.__len__()-1]:
	for i in indices:
		key1 = 'o'+str((i)%NN)+'h'+str(j)
		key2 = 'o'+str((i+1)%NN)+'h'+str(j)
		mywfile.write("    " + writeVFace(key1, key2) + "\n")
	# end for
# end for


endPatch(mywfile)

mywfile.write(listEnd)
mywfile.write(newLine)

mywfile.write(seperator)

# finished
mywfile.close()

# ----------------------------------------------
#					T H E   E N D
# ----------------------------------------------

