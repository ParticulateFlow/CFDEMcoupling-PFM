cd db$1
cd CFD
cp -r orig.0 0
blockMesh
decomposePar -force
mpirun -np 4 cfdemSolverRhoPimple -parallel > createDB.log 2>&1

reconstructPar -time 2.5:
# remove uniform folder
find . -name \uniform -type d -exec rm -rf {} \;
# make database
mkdir dataBase$1

mv [2-5]* dataBase$1/

# shift times
cd dataBase$1
for t in *; do
        t0=$t
        break 1
done

for t in *; do
    tnew=$(awk "BEGIN {print $t-$t0}")
    mv $t $tnew
done

for t in *; do
    cd $t
    mv NuFieldMean NuField
    mv phiMean phi
    mv pMean p
    mv UsMean Us
    mv voidfractionMean voidfraction
    rm UMean
    rm rhoMean
    cd ..
done

cd ..
#cp -r constant dataBase$1/
#cp -r system dataBase$1/
