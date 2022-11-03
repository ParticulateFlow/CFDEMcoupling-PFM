rm -r postProcessing
rm -r dynamicCode

for proc in processor*
do
    cd $proc
    rm -r [1-9]*
    rm -r 0.*
    cd ..
done
