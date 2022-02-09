#!/bin/bash

#===================================================================#
# allrun sequential script
# Tim MJ Nijssen - September 2021
#===================================================================#

#--------------------------------------------------------------------------------#
#- run settings file
scriptPath="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
. $scriptPath/scriptSettings.sh $1
#--------------------------------------------------------------------------------#

#- find sequenceSettings file
if [ -z "$2" ]; then
    #- current path if not other specified
    sequenceSettingsPath=$casePath/sequenceSettings
else
    #- read input
    sequenceSettingsPath="$(dirname "$(readlink -f $2)/.")"
fi

#- other paths
boundaryFilesPath=$casePath/boundary_files

#- do all tasks before actual run
perl $scriptPath/updateControlDict.pl $casePath/CFD/system/controlDict $casePath/CFD/system/controlDict 0 ' ' # set startTime to 0 before decomposing
. $scriptPath/preRunAllPar.sh $casePath

echo "RunAllParSequence: duplicating LIGGGHTS restart file"
cp $casePath/DEM/post/restart/liggghts.restart $casePath/DEM/post/restart/liggghts.restartSequence # duplicate restart file to restartSequence

#- read sequenceSettings file line by line
while IFS= read -r stepLine || [[ -n "$stepLine" ]]; do
    if [ ! -z "$stepLine" ]; then                # skip empty lines
	if [ "${stepLine:0:1}" != "#" ]; then        # skip bash style comments
	    if [ "${stepLine:0:2}" != "//" ] ; then  # skip c style comments
	        stepWords=( $stepLine );
		stepName=${stepWords[0]}
		startTime=${stepWords[1]}
		endTime=${stepWords[2]}

		echo "runAllParSequence: running step $stepName from t=$startTime to t=$endTime"

		# update start and end times in controlDict
		perl $scriptPath/updateControlDict.pl $casePath/CFD/system/controlDict $casePath/CFD/system/controlDict $startTime $endTime

		# find and update boundary conditions
		for timeDir in $casePath/CFD/processor0/*/ ; do             # loop over time directories
		    timeName=$(basename $timeDir)                           # get time folder name

		    # check if numeric
		    re='^[0-9]+([.][0-9]+)?$'
		    if [[ $timeName =~ $re ]] ; then

			# find startTime folder
			if (( $(echo "$timeName == $startTime" |bc -l) )); then

			    for rFile in $boundaryFilesPath/$stepName/* ; do                          # loop over replacement files
				if [ "${rFile: -1}" != "~" ]; then                                        # skip temp files
				    
				    rFileName=$(basename $rFile)                                          # get file name
				    echo "runAllParSequence: updating boundary conditions in processor*/$timeName/$rFileName"
				    for procDir in $casePath/CFD/processor* ; do                          # loop over processor directories
				        oFile="$procDir/$timeName/$rFileName"                             # path to original file in CFD/processorX/time
				        perl $scriptPath/replaceBC.pl $oFile $rFile $oFile                # update boundary conditions
				    done
				fi
		            done
			fi
		    fi
		done
		#- run parallel CFD-DEM
		echo "runAllParSequence: Running CFD-DEM"
		. $scriptPath/runCFDDEMPar.sh $casePath </dev/null
	    fi
	fi
    fi
done < "$sequenceSettingsPath"    

#- do all tasks after actual run
perl $scriptPath/updateControlDict.pl $casePath/CFD/system/controlDict $casePath/CFD/system/controlDict 0 ' ' # set startTime to 0 before reconstructing
. $scriptPath/postRunAllPar.sh $casePath



