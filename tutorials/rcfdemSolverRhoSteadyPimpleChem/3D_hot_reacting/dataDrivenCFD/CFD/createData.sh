reconstructPar -latestTime -fields '(p rho T TMean U UMean partTempMean voidfractionMean voidfraction)'
rm -r VTK
foamToVTK -latestTime
rm -r [1-9]*
