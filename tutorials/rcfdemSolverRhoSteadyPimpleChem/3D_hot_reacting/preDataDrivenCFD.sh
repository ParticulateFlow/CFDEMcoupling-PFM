python updateVelScaling.py

cd CFDDEM/CFD
reconstructPar -latestTime
latestTime=$(foamListTimes -noZero -latestTime)

symmetrizer -dict "symmetrizerPropertiesPost_1"
symmetrizer -dict "symmetrizerPropertiesPost_2"
symmetrizer -dict "symmetrizerPropertiesPost_3"
symmetrizer -dict "symmetrizerPropertiesPost_4"

cd ../../dataDrivenCFD/CFD

typeset -i currIndex=$(cat "../../runIndex")
mv postProcessing postProcessing_${currIndex}
mv dataBase dataBase_${currIndex}
mkdir dataBase
mkdir dataBase/0

cp ../../CFDDEM/CFD/${latestTime}/voidfractionMean dataBase/0/voidfraction
cp ../../CFDDEM/CFD/${latestTime}/UsMean dataBase/0/Us
cp ../../CFDDEM/CFD/${latestTime}/UMean dataBase/0/UMean

cp -r dataBase/0 dataBase/2
cp -r dataBase/0 dataBase/4

latestTime=$(foamListTimes -noZero -latestTime -processor)

if [ -z "$latestTime" ]
then
  echo "No time folder beyond 0 found, assuming first run."
  decomposePar
else
  symmetrizer -dict "symmetrizerProperties_1"
  symmetrizer -dict "symmetrizerProperties_2"
  symmetrizer -dict "symmetrizerProperties_3"
  symmetrizer -dict "symmetrizerProperties_4"
  decomposePar -fields -latestTime
  for proc in processor*
    do
      cp $proc/${latestTime}/partTemp $proc/${latestTime}/partTempRef
    done
fi


newLatestTime=$(awk -vlatestTime=${latestTime} 'BEGIN{printf "%.1f\n",latestTime+36000}')
echo "Updating endTime to ${newLatestTime}."
sed -i "/^endTime/c\endTime         ${newLatestTime};" system/controlDict

./decomposeParDB.sh

cp ../../CFDDEM/DEM/liggghts.restart ../../CFDDEM/DEM/liggghts.restart_${currIndex}
mv ../../CFDDEM/DEM/liggghts.restartCFDEM ../../CFDDEM/DEM/liggghts.restart

cp ../../CFDDEM/DEM/initOre ../DEM/initOre
cp ../../CFDDEM/DEM/initCoke ../DEM/initCoke
cp ../../CFDDEM/DEM/initFineMix ../DEM/initFineMix
cp ../../CFDDEM/DEM/initFullState ../DEM/initFullState

mv ../../CFDDEM/DEM/initOre ../../CFDDEM/DEM/initOre_${currIndex}
mv ../../CFDDEM/DEM/initCoke ../../CFDDEM/DEM/initCoke_${currIndex}
mv ../../CFDDEM/DEM/initFineMix ../../CFDDEM/DEM/initFineMix_${currIndex}
mv ../../CFDDEM/DEM/initFullState ../../CFDDEM/DEM/initFullState_${currIndex}
cp ../../CFDDEM/DEM/zInsLower ../../CFDDEM/DEM/zInsLower_${currIndex}

cd ../DEM
mv post post_${currIndex}
mkdir post

mv restart restart_${currIndex}
mkdir restart

let "currIndex += 1"
echo ${currIndex}
echo ${currIndex} > ../../runIndex

cd ../..
