cd rCFD

cd CFD
rm -r proc*
rm -r dynamicCode
rm -r postProcessing*
rm -r clockData
rm -r [0-9]*
rm -r [0-9]*
rm log.*
rm *.log
rm recurrenceError
rm recurrenceMatrix
rm recurrencePath


cd ../DEM
cd post
rm *
cd ..
rm *.txt
rm liggghts.restartCFDEM*
cd ../..
