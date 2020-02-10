import os
import shutil

with open('CFDDEM/CFD/postProcessing/cellObj1/0/volRegion.dat', 'r') as f:
    for line in f:
        if line.startswith('#'): 
            continue
        row = line.split()
    v = float(row[3].rstrip(")"))

scaleFac = -v/0.002

f = open("dataDrivenCFD/CFD/constant/tmp","w+")
with open("dataDrivenCFD/CFD/constant/couplingProperties") as old_file:
        for line in old_file:
            row = line.split()
            if (line != "\n"):
                if (row[0]=="scalingFactor"):
                    f.write("    scalingFactor " + str(scaleFac) + ";\n")
                else:
                    f.write(line)
            else:
                f.write(line)

shutil.move("dataDrivenCFD/CFD/constant/tmp", "dataDrivenCFD/CFD/constant/couplingProperties")
