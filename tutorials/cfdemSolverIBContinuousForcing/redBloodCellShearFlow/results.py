#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt

# Reading sphere data file
t,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,x9,y9,z9,x10,y10,z10 = np.genfromtxt('DEM/post/particle_position.txt', skip_header=0,usecols = (0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30),unpack=True)
t = t*1e-06

# Calculating the centre of mass of RBC
xcm = (x1+x2+x3+x4+x5+x6+x7+x8+x9+x10)/10
ycm = (y1+y2+y3+y4+y5+y6+y7+y8+y9+y10)/10
zcm = (z1+z2+z3+z4+z5+z6+z7+z8+z9+z10)/10

# Calculating the gyration tensor
G11 = ((x1-xcm)*(x1-xcm) + (x2-xcm)*(x2-xcm) + (x3-xcm)*(x3-xcm) + (x4-xcm)*(x4-xcm) + (x5-xcm)*(x5-xcm) + (x6-xcm)*(x6-xcm) + (x7-xcm)*(x7-xcm) + (x8-xcm)*(x8-xcm) + (x9-xcm)*(x9-xcm) + (x10-xcm)*(x10-xcm))/10
G12 = ((x1-xcm)*(y1-ycm) + (x2-xcm)*(y2-ycm) + (x3-xcm)*(y3-ycm) + (x4-xcm)*(y4-ycm) + (x5-xcm)*(y5-ycm) + (x6-xcm)*(y6-ycm) + (x7-xcm)*(y7-ycm) + (x8-xcm)*(y8-ycm) + (x9-xcm)*(y9-ycm) + (x10-xcm)*(y10-ycm))/10
G13 = ((x1-xcm)*(z1-zcm) + (x2-xcm)*(z2-zcm) + (x3-xcm)*(z3-zcm) + (x4-xcm)*(z4-zcm) + (x5-xcm)*(z5-zcm) + (x6-xcm)*(z6-zcm) + (x7-xcm)*(z7-zcm) + (x8-xcm)*(z8-zcm) + (x9-xcm)*(z9-zcm) + (x10-xcm)*(z10-zcm))/10
G21 = ((y1-ycm)*(x1-xcm) + (y2-ycm)*(x2-xcm) + (y3-ycm)*(x3-xcm) + (y4-ycm)*(x4-xcm) + (y5-ycm)*(x5-xcm) + (y6-ycm)*(x6-xcm) + (y7-ycm)*(x7-xcm) + (y8-ycm)*(x8-xcm) + (y9-ycm)*(x9-xcm) + (y10-ycm)*(x10-xcm))/10
G22 = ((y1-ycm)*(y1-ycm) + (y2-ycm)*(y2-ycm) + (y3-ycm)*(y3-ycm) + (y4-ycm)*(y4-ycm) + (y5-ycm)*(y5-ycm) + (y6-ycm)*(y6-ycm) + (y7-ycm)*(y7-ycm) + (y8-ycm)*(y8-ycm) + (y9-ycm)*(y9-ycm) + (y10-ycm)*(y10-ycm))/10
G23 = ((y1-ycm)*(z1-zcm) + (y2-ycm)*(z2-zcm) + (y3-ycm)*(z3-zcm) + (y4-ycm)*(z4-zcm) + (y5-ycm)*(z5-zcm) + (y6-ycm)*(z6-zcm) + (y7-ycm)*(z7-zcm) + (y8-ycm)*(z8-zcm) + (y9-ycm)*(z9-zcm) + (y10-ycm)*(z10-zcm))/10
G31 = ((z1-zcm)*(x1-xcm) + (z2-zcm)*(x2-xcm) + (z3-zcm)*(x3-xcm) + (z4-zcm)*(x4-xcm) + (z5-zcm)*(x5-xcm) + (z6-zcm)*(x6-xcm) + (z7-zcm)*(x7-xcm) + (z8-zcm)*(x8-xcm) + (z9-zcm)*(x9-xcm) + (z10-zcm)*(x10-xcm))/10
G32 = ((z1-zcm)*(y1-ycm) + (z2-zcm)*(y2-ycm) + (z3-zcm)*(y3-ycm) + (z4-zcm)*(y4-ycm) + (z5-zcm)*(y5-ycm) + (z6-zcm)*(y6-ycm) + (z7-zcm)*(y7-ycm) + (z8-zcm)*(y8-ycm) + (z9-zcm)*(y9-ycm) + (z10-zcm)*(y10-ycm))/10
G33 = ((z1-zcm)*(z1-zcm) + (z2-zcm)*(z2-zcm) + (z3-zcm)*(z3-zcm) + (z4-zcm)*(z4-zcm) + (z5-zcm)*(z5-zcm) + (z6-zcm)*(z6-zcm) + (z7-zcm)*(z7-zcm) + (z8-zcm)*(z8-zcm) + (z9-zcm)*(z9-zcm) + (z10-zcm)*(z10-zcm))/10

# Calculating the eigen values of the gyration tensor
min_eig = []

for i in range(len(t)):
	G = np.matrix([ [G11[i],G12[i],G13[i]],
			[G21[i],G22[i],G23[i]],
			[G31[i],G32[i],G33[i]] ])

	eigenvalues = np.linalg.eigvals(G)
	eig = np.min(eigenvalues)
	min_eig.append(eig)

plt.plot(t,min_eig,linewidth = 2)
plt.xlabel('Time (s)')
plt.ylabel('Min. eigen value of G')
plt.grid()
plt.show()

plt.plot(z1,y1)
plt.plot(z2,y2)
plt.plot(z3,y3)
plt.plot(z4,y4)
plt.plot(z5,y5)
plt.plot(z6,y6)
plt.plot(z7,y7)
plt.plot(z8,y8)
plt.plot(z9,y9)
plt.plot(z10,y10)
plt.show()

plt.plot(z1,x1)
plt.plot(z2,x2)
plt.plot(z3,x3)
plt.plot(z4,x4)
plt.plot(z5,x5)
plt.plot(z6,x6)
plt.plot(z7,x7)
plt.plot(z8,x8)
plt.plot(z9,x9)
plt.plot(z10,x10)
plt.show()
