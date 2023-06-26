cd dataDrivenCFD/CFD
reconstructPar -latestTime
latestTime=$(foamListTimes -noZero -latestTime)

if [ -z "$latestTime" ]
then
  echo "No time folder beyond 0 found, assuming first run."
  latestTime=0
  cp orig.0/partTemp ../../CFDDEM/CFD/0/partTempRef
  cp orig.0/partTemp ../../CFDDEM/CFD/0/partTemp
else
  cp ${latestTime}/partTemp ../../CFDDEM/CFD/0/partTempRef
  cp ${latestTime}/partTemp ../../CFDDEM/CFD/0/partTemp
fi

cd ../../CFDDEM/CFD
rm -r proc*
rm ABORT
latestTime2=$(foamListTimes -noZero -latestTime)
if [ -z "$latestTime2" ]
then
  echo "No time folder beyond 0 found in CFDDEM, assuming first run."
else
  echo "Removing latestTime folder ${latestTime2} in CFDDEM"
  rm -r ${latestTime2}
fi

symmetrizer -dict "symmetrizerPropertiesPre_1"
symmetrizer -dict "symmetrizerPropertiesPre_2"
symmetrizer -dict "symmetrizerPropertiesPre_3"
symmetrizer -dict "symmetrizerPropertiesPre_4"

decomposePar

if [ ${latestTime} -ne "0" ]
then
  typeset -i currIndex=$(cat "../../runIndex")
  mv postProcessing postProcessing_${currIndex}
  cd ../DEM
  mv post post_${currIndex}
  mv monitor monitor_${currIndex}
  mkdir post
  mkdir monitor
fi

cd ../..
