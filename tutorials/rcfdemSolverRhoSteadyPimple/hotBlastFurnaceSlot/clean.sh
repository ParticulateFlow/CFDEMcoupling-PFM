rm runIndex
cd CFDDEM/CFD
rm -r proc*
rm -r 0
cd ../DEM
rm -r post*
mkdir post
rm liggghts.restart

cd ../../dataDrivenCFD/CFD
rm -r proc*
rm -r 0
rm ../DEM/temp_ave.txt
