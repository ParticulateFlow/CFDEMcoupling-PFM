voidfractionRec=particleCloud.recM().voidfraction()();
URec=particleCloud.recM().U()();
UsRec=particleCloud.recM().Us()();
phiRec=linearInterpolate(URec*voidfractionRec) & mesh.Sf();
