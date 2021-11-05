#!/bin/bash

#===================================================================#
# allrun conditional sequential script
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
    tappingSettingsPath=$casePath/tappingSettings
else
    #- read input
    tappingSettingsPath="$(dirname "$(readlink -f $2)/.")"
fi

#- other paths
boundaryFilesPath=$casePath/boundary_files

#- search for startTime setting
sequenceStartTime=0
while IFS= read -r line || [[ -n "$line" ]]; do
    if [ ! -z "$line" ]; then                      # skip empty lines
	if [ "${line:0:1}" != "#" ]; then              # skip bash style comments
	    if [ "${line:0:2}" != "//" ] ; then        # skip c style comments
		if [ "${line:0:9}" == "startTime" ] ; then # read startTime from file
		    timeWords=( $line );
		    timeSetting=${timeWords[1]}

		    # check if numeric
		    re='^[0-9]+([.][0-9]+)?$'
		    if [[ $timeSetting =~ $re ]] ; then
			echo "RunAllParSequence: found startTime=$timeSetting"
			sequenceStartTime=$timeSetting
			break
		    else
			echo "RunAllParSequence: ERROR invalid startTime setting: $timeSetting"
			return 1
		    fi
		fi
	    fi
	fi
    fi
done < "$tappingSettingsPath"		    

#- do all tasks before actual run
echo "RunAllParSequence: setting time to $sequenceStartTime"
perl $scriptPath/updateControlDict.pl $casePath/CFD/system/controlDict $casePath/CFD/system/controlDict $sequenceStartTime ' ' # set startTime before decomposing
currentTime=$sequenceStartTime
. $scriptPath/preRunAllPar.sh $casePath

if (( $sequenceStartTime == 0 )); then
    echo "RunAllParSequence: duplicating LIGGGHTS restart file"
    cp $casePath/DEM/post/restart/liggghts.restart $casePath/DEM/post/restart/liggghts.restartSequence # duplicate restart file to restartSequence
else
    echo "RunAllParSequence: starting from existing restartSequence file"
fi

#- read sequenceSettings file line by line
while IFS= read -r stepLine || [[ -n "$stepLine" ]]; do
    if [ ! -z "$stepLine" ]; then                    # skip empty lines
	if [ "${stepLine:0:1}" != "#" ]; then        # skip bash style comments
	    if [ "${stepLine:0:2}" != "//" ] ; then  # skip c style comments
		if [ "${stepLine:0:9}" != "startTime" ] ; then # skip startTime line
		    
	            stepWords=( $stepLine );
		    stepName=${stepWords[0]}
		    variable=${stepWords[1]}
		    operator=${stepWords[2]}
		    value=${stepWords[3]}
		    interval=${stepWords[4]}
		    nRunsMax=${stepWords[5]}

		    echo "runAllParTapping: running step $stepName until $variable $operator $value"

		    #- find and update boundary conditions
		    for timeDir in $casePath/CFD/processor0/*/ ; do         # loop over time directories
			timeName=$(basename $timeDir)                           # get time folder name

			# check if numeric
			re='^[0-9]+([.][0-9]+)?$'
			if [[ $timeName =~ $re ]] ; then

			    # find startTime folder
			    if (( $(echo "$timeName == $currentTime" |bc -l) )); then

				for rFile in $boundaryFilesPath/$stepName/* ; do                          # loop over replacement files
				    if [ "${rFile: -1}" != "~" ]; then                                    # skip temp files
					
					rFileName=$(basename $rFile)                                          # get file name
					echo "runAllParTapping: updating boundary conditions in processor*/$timeName/$rFileName"
					for procDir in $casePath/CFD/processor* ; do                          # loop over processor directories
				            oFile="$procDir/$timeName/$rFileName"                         # path to original file in CFD/processorX/time
				            perl $scriptPath/replaceBC.pl $oFile $rFile $oFile            # update boundary conditions
					done
				    fi
				done
			    fi
			fi
		    done
		    
		    #- loop up to the maximum number of runs
		    for (( iRun=1; iRun<=$nRunsMax; iRun++ )); do
			startTime=$currentTime
			endTime=$(echo $currentTime + $interval | bc)

			# update start and end times in controlDict
			perl $scriptPath/updateControlDict.pl $casePath/CFD/system/controlDict $casePath/CFD/system/controlDict $startTime $endTime
			
			#- run parallel CFD-DEM
			echo "runAllParTapping: Running CFD-DEM (nRuns = $iRun)"
			. $scriptPath/runCFDDEMPar.sh $casePath </dev/null
			
			#- update current time
			currentTime=$endTime

			#- check for maximum number of runs
			if (( $iRun == $nRunsMax )); then
			    echo "runAllParTapping: maximum number of runs ($nRunsMax) reached, proceeding to next step"
			    break
			fi

			#- find the variable
			variablePath="$casePath/CFD/postProcessing/$variable.dat"
			if [ ! -f "$variablePath" ]; then
			    echo "runAllParTapping: ERROR $variablePath does not exist"
			    return 1
			fi
			
			#- read variable
			lastLine=$( tail -n 1 $variablePath )
			lastWord=`echo ${lastLine##* }`
			echo "runAllParTapping: $variable = $lastWord"

			#- check condition
			condition="$lastWord $operator $value"
			goToNextStep=$(echo $condition | bc -l )

			#- do another run or skip to next step
			if (( $goToNextStep == 1 )); then
			    echo "runAllParTapping: $variable $operator $value returned TRUE, proceeding to next step"
			    break
			else
			    echo "runAllParTapping: $variable $operator $value returned FALSE, doing another run"
			fi
			
		    done
		fi
	    fi
	fi
    fi
done < "$tappingSettingsPath"

#- do all tasks after actual run
perl $scriptPath/updateControlDict.pl $casePath/CFD/system/controlDict $casePath/CFD/system/controlDict $sequenceStartTime ' ' # set startTime before reconstructing
. $scriptPath/postRunAllPar.sh $casePath



