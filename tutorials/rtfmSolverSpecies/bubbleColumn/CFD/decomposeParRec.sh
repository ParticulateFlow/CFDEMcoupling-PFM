#!/bin/sh
# Source run functions
. $WM_PROJECT_DIR/bin/tools/RunFunctions

if [ $# -eq 0 ]
  then
    dBname="dataBase"
  else
    dBname=$1
fi

rm log.decomposePar

runApplication decomposePar -force

cpdirs="system constant"
for f in $cpdirs
do
  cp -r $f $dBname/$f
done
cd $dBname
rm log.decomposePar
runApplication decomposePar -force -time 0:

rm -rf $cpdirs
rm -rf processor*/constant

for proc in processor*
do
    echo "Transferring decomposed recurrence fields of $dBname/$proc to $proc."
    mkdir ../$proc/$dBname
    mv $proc/* ../$proc/$dBname/
    rm -rf $proc
done
