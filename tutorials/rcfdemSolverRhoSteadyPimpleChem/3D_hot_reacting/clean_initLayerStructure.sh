cd initLayerStructure/CFD
rm -r proc*
rm -r [0-9]*
rm -r postProcessing
rm -r dynamicCode
rm -r clockData
rm *.log
rm log.*

cd ../DEM
rm -r post*
rm -r monitor
rm init*
mkdir post
mkdir monitor
rm liggghts.restart*
rm restart.*
rm log.*
