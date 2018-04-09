#!/bin/sh
# Source run functions
. $WM_PROJECT_DIR/bin/tools/RunFunctions

# to be executed from one level above the source directory

if [ $# -eq 0 ]
  then
    sourceName="dataBase"
    targetName="dataBaseCoarse"
  else
    sourceName=$1
    targetName=$2
fi

cd $sourceName

for time in *
do
    if [ $time != "system" ] && [ $time != "constant" ];
    then
        cd ../$targetName
        echo "Found $time."
        sed -i "/^startTime/c\startTime \t$time;" ./system/controlDict
        grep 'startTime' ./system/controlDict
        mapFields ../$sourceName -sourceTime $time -consistent
        cd ../$sourceName
    fi
done

