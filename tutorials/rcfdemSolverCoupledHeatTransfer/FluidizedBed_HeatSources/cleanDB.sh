cd db$1
cd CFD
rm -r proc*
rm -r dynamicCode
rm -r postProcessing*
rm -r clockData
rm -r [0-9]*
rm -r -[0-9]*
rm -r dataBase*
rm log.*
rm *.log

cd ../DEM
cd post
rm *
cd ..
rm *.txt
rm liggghts.restartCFDEM*
cd ../..
