#!/bin/sh
# Source run functions
. $WM_PROJECT_DIR/bin/tools/RunFunctions

echo "Clearing processor directories without removing decomposed database."

for proc in processor*
do
    echo "Clearing $proc directory."
    rm -rf $proc/0.*
    rm -rf $proc/[1-9]*
done
