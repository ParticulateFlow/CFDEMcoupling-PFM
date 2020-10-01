python updateVelScaling.py

cd CFDDEM/CFD
reconstructPar -latestTime
latestTime=$(foamListTimes -noZero -latestTime)

symmetrizer -dict "symmetrizerPropertiesPost"

cd ../../dataDrivenCFD/CFD

typeset -i currIndex=$(cat "../../runIndex")
mv dataBase dataBase_${currIndex}
mkdir dataBase
mkdir dataBase/0

cp ../../CFDDEM/CFD/${latestTime}/voidfractionMean dataBase/0/voidfraction
cp ../../CFDDEM/CFD/${latestTime}/UsMean dataBase/0/Us
cp ../../CFDDEM/CFD/${latestTime}/UMean dataBase/0/UMean

cp -r dataBase/0 dataBase/1
cp -r dataBase/0 dataBase/2

latestTime=$(foamListTimes -noZero -latestTime -processor)

if [ -z "$latestTime" ]
then
  echo "No time folder beyond 0 found, assuming first run."
  decomposePar
else
  symmetrizer -dict "symmetrizerProperties"
  decomposePar -fields -latestTime
  for proc in processor*
    do
      cp $proc/${latestTime}/partTemp $proc/${latestTime}/partTempRef
    done
fi


let "newLatestTime=${latestTime}+18000"
echo "Updating endTime to ${newLatestTime}."
sed -i "/^endTime/c\endTime         ${newLatestTime};" system/controlDict

./decomposeParDB.sh

cp ../../CFDDEM/DEM/liggghts.restart ../../CFDDEM/DEM/liggghts.restart_${currIndex}
mv ../../CFDDEM/DEM/liggghts.restartCFDEM ../../CFDDEM/DEM/liggghts.restart
cp ../../CFDDEM/DEM/liggghts.restart ../DEM/liggghts.restart
cp ../../CFDDEM/DEM/initOre ../DEM/initOre
cp ../../CFDDEM/DEM/initCoke ../DEM/initCoke
cp ../../CFDDEM/DEM/initFineMix ../DEM/initFineMix

let "currIndex += 1"
echo ${currIndex}
echo ${currIndex} > ../../runIndex

cd ../DEM
mv post post_${currIndex}
mkdir post

cd ../..
