URec=particleCloud.recM().tU();
voidfractionRec=particleCloud.recM().tvoidfraction();
UsRec=particleCloud.recM().tUs();
phiRec=linearInterpolate(URec*voidfractionRec) & mesh.Sf();
