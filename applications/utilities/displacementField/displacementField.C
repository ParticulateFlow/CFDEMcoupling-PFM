/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    displacementField

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "vectorList.H"
#include <sstream>
#include <string>
#include <sys/stat.h>
#include <unistd.h>
#include <set>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
void findPairs(labelList &, labelList &, labelPairList &);
void findPairsUnordered(labelList &, labelList &, labelPairList &);
void fillEmptyCells(fvMesh &, label, label, scalarList &, volVectorField &, volVectorField &, scalarList &, volVectorField &, volVectorField &, bool, scalar);
void nearestNeighborCells(fvMesh &, label, label, label, scalarList &, labelList &);
void normalizeFields(scalarList &, volVectorField &, volVectorField &);
void readDump(std::string, labelList &, scalarList &, vectorList &);
scalar weightFun(scalar);
label maxNumParticles = 1000000;
scalar minVol = 1e-12;
scalar Pi43 = 4.1888;
label posIndex = -1;
label posRad = -1;
label posX = -1;
label posY = -1;
label posZ = -1;

int main(int argc, char *argv[])
{
    argList::addOption
    (
        "totalProcs",
        "label",
        "total number of parallel processes, defaults to 1"
    );
    argList::addOption
    (
        "thisProc",
        "label",
        "number of current process, defaults to 0"
    );


    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    const label thisProc = args.optionLookupOrDefault("thisProc", 0);
    const label totalProcs = args.optionLookupOrDefault("totalProcs", 1);

    Info << "This is number " << thisProc << " of " << totalProcs << " processes." << endl;


    // user-defined input for each case
    IOdictionary displacementProperties
    (
        IOobject
        (
            "displacementProperties",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    label dumpIndexStart(readLabel(displacementProperties.lookup("dumpIndexStart")));
    label dumpIndexEnd(readLabel(displacementProperties.lookup("dumpIndexEnd")));
    label dumpIndexInputIncrement(readLabel(displacementProperties.lookup("dumpIndexInputIncrement")));
    label dumpIndexDisplacementIncrement(readLabel(displacementProperties.lookup("dumpIndexDisplacementIncrement")));
    label nNeighMin(readLabel(displacementProperties.lookup("nNeighMin")));
    label maxSearchLayers(displacementProperties.lookupOrDefault<label>("maxSearchLayers",0));
    posIndex = readLabel(displacementProperties.lookup("posIndex"));
    posRad = readLabel(displacementProperties.lookup("posRad"));
    posX = readLabel(displacementProperties.lookup("posX"));
    posY = readLabel(displacementProperties.lookup("posY"));
    posZ = readLabel(displacementProperties.lookup("posZ"));
    scalar timePerInputStep(readScalar(displacementProperties.lookup("timePerInputStep")));
    scalar timePerDisplacementStep(readScalar(displacementProperties.lookup("timePerDisplacementStep")));
    scalar startTime(readScalar(displacementProperties.lookup("startTime")));
    std::string filepath=string(displacementProperties.lookup("filepath"));
    std::string fileext=string(displacementProperties.lookupOrDefault<string>("fileextension",""));
    bool interpolate=bool(displacementProperties.lookupOrDefault<bool>("fillEmptyCells",true));
    bool averageMode=bool(displacementProperties.lookupOrDefault<bool>("averageMode",false));

    volVectorField defaultUs
    (
        IOobject
        (
            "defaultUDisp",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedVector("zero", dimensionSet(0,1,-1,0,0), vector::zero)
    );

    volVectorField defaultUsDirectedStdDev
    (
        IOobject
        (
            "defaultUDispDirectedStdDev",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedVector("zero", dimensionSet(0,1,-1,0,0), vector::zero)
    );

    scalar xmin=scalar(displacementProperties.lookupOrDefault<scalar>("xmin",-1e10));
    scalar xmax=scalar(displacementProperties.lookupOrDefault<scalar>("xmax",1e10));
    scalar ymin=scalar(displacementProperties.lookupOrDefault<scalar>("ymin",-1e10));
    scalar ymax=scalar(displacementProperties.lookupOrDefault<scalar>("ymax",1e10));
    scalar zmin=scalar(displacementProperties.lookupOrDefault<scalar>("zmin",-1e10));
    scalar zmax=scalar(displacementProperties.lookupOrDefault<scalar>("zmax",1e10));
    scalarList boundaries(6);
    boundaries[0]=xmin;
    boundaries[1]=xmax;
    boundaries[2]=ymin;
    boundaries[3]=ymax;
    boundaries[4]=zmin;
    boundaries[5]=zmax;

    vectorList probePoints=vectorList(displacementProperties.lookupOrDefault<vectorList>("probePoints",vectorList(0)));
    bool monitorProbes = false;
    if (probePoints.size()>0) monitorProbes = true;
    #include "OFstream.H"
    OFstream monitoringDataFile("monitoringData.txt");
    if (monitorProbes)
    {
        monitoringDataFile << "# monitoring data file" << endl;
        monitoringDataFile << "# format: time nPerCell[p1] UDisp[p1] UDispDirectedVariance[p1] nPerCell[p2] ... " << endl;
        for(label p=0;p<probePoints.size();p++)
        {
            vector pos = probePoints[p];
            monitoringDataFile << "# point[" << p << "] = " << pos << endl;
        }
    }


    label dumpIndex1 = dumpIndexStart + thisProc * dumpIndexInputIncrement;
    label dumpIndex2 = dumpIndex1 + dumpIndexDisplacementIncrement;

    volVectorField Us
    (
        IOobject
        (
            "UDisp",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedVector("zero", dimensionSet(0,1,-1,0,0), vector::zero)
    );

    volVectorField UsDirectedStdDev
    (
        IOobject
        (
            "UDispDirectedStdDev",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedVector("zero", dimensionSet(0,1,-1,0,0), vector::zero)
    );

    scalarList particleVolInCell(mesh.nCells(), 0.0);

    scalar currTime=startTime + thisProc * timePerInputStep;
    label timeIndex=thisProc;

    while(true)
    {
        runTime.setTime(currTime,timeIndex);
        // read dump files and check which particle indices are present in both
        labelList indices1, indices2;
        scalarList radii1, radii2;
        vectorList positions1, positions2;

        std::stringstream ss;
        ss << filepath << dumpIndex1 << fileext;
        std::string filename1 = ss.str();
        ss.str("");
        ss << filepath << dumpIndex2 << fileext;
        std::string filename2 = ss.str();

        if (access( filename1.c_str(), F_OK ) == -1 || access( filename2.c_str(), F_OK ) == -1 || dumpIndex2 > dumpIndexEnd)
        {
            if (averageMode)
            {
                normalizeFields(particleVolInCell, Us, UsDirectedStdDev);
                fillEmptyCells(mesh,nNeighMin,maxSearchLayers,particleVolInCell,Us,UsDirectedStdDev,boundaries,defaultUs,defaultUsDirectedStdDev,interpolate,timePerDisplacementStep);

                Us /= timePerDisplacementStep;
                UsDirectedStdDev /= timePerDisplacementStep;
                Us.write();
                UsDirectedStdDev.write();
            }
            break;
        }

        Info << "\nReading" << endl;
        Info << "\t" << filename1 << endl;
        Info << "\t" << filename2 << endl;
        Info << "corresponding to time = " << currTime << "." << endl;

        readDump(filename1, indices1, radii1, positions1);
        readDump(filename2, indices2, radii2, positions2);

        labelPairList pairs;
        findPairs(indices1,indices2,pairs);

        // average particle displacements and their variance
        Info << "Binning particle displacements on mesh." << endl;
        vector position, displacement;
        scalar radius, volume;
        label line1, line2;
        label cellI;

        if (!averageMode)
        {
            Us *= 0.0;
            UsDirectedStdDev *= 0.0;
            particleVolInCell.clear();
            particleVolInCell.setSize(mesh.nCells(), 0);
        }

        for (label partI = 0; partI < pairs.size(); partI++)
        {
            line1 = pairs[partI].first();
            line2 = pairs[partI].second();
            position = positions1[line1];
            cellI = mesh.findCell(position);
            if (cellI < 0)  continue;
            displacement = positions2[line2] - positions1[line1];
            radius = radii1[line1];
            volume = Pi43 * radius * radius * radius;
            particleVolInCell[cellI] += volume;
            Us[cellI] += displacement*volume;

            for (label comp=0;comp<3;comp++)
            {
                UsDirectedStdDev[cellI].component(comp) += displacement.component(comp)*displacement.component(comp)*volume;
            }
        }

        if (!averageMode)
        {
            normalizeFields(particleVolInCell, Us, UsDirectedStdDev);
            fillEmptyCells(mesh,nNeighMin,maxSearchLayers,particleVolInCell,Us,UsDirectedStdDev,boundaries,defaultUs,defaultUsDirectedStdDev,interpolate,timePerDisplacementStep);

            Us /= timePerDisplacementStep;
            UsDirectedStdDev /= timePerDisplacementStep;
            Us.write();
            UsDirectedStdDev.write();
        }

        if (averageMode && monitorProbes)
        {
            monitoringDataFile << currTime << " ";
            for(label p=0;p<probePoints.size();p++)
            {
                vector pos = probePoints[p];
                label cellP = mesh.findCell(pos);
                monitoringDataFile << " " << particleVolInCell[cellP] << " " << Us[cellP]/timePerDisplacementStep << " " << UsDirectedStdDev[cellP]/(timePerDisplacementStep*timePerDisplacementStep);
            }
            monitoringDataFile << endl;
        }

        dumpIndex1 += dumpIndexInputIncrement*totalProcs;
        dumpIndex2 += dumpIndexInputIncrement*totalProcs;
        currTime += timePerInputStep*totalProcs;
        timeIndex += totalProcs;
    }
    return 0;
}

void readDump(std::string filename, labelList &indices, scalarList &radii, vectorList &positions)
{
    #include <fstream>

    const label leadingLines = 9;
    label lineCounter = 0;
    label partIndex = 0;
    scalar r = 1.0, x = 0.0, y = 0.0, z = 0.0;

    indices.clear();
    radii.clear();
    positions.clear();

    indices.setSize(maxNumParticles);
    radii.setSize(maxNumParticles);
    positions.setSize(maxNumParticles);

    std::ifstream file(filename);
    std::string str;
    std::string word;
    label wordcounter;
    while (std::getline(file, str))
    {
        if (lineCounter >= leadingLines)
        {
            std::istringstream ss(str);
            wordcounter = 0;
            while (ss >> word)
            {
                if (wordcounter == posIndex)
                {
                    partIndex = stoi(word);
                }
                else if (wordcounter == posRad)
                {
                    r = stod(word);
                }
                else if (wordcounter == posX)
                {
                    x = stod(word);
                }
                else if (wordcounter == posY)
                {
                    y = stod(word);
                }
                else if (wordcounter == posZ)
                {
                    z = stod(word);
                }
                wordcounter++;
            }
//            sscanf(str.c_str(), "%d %lf %lf %lf", &partIndex, &x, &y, &z);
            indices[lineCounter-leadingLines] = partIndex;
            radii[lineCounter-leadingLines] = r;
            positions[lineCounter-leadingLines] = vector(x,y,z);
        }
        lineCounter++;
    }

    label readLines = lineCounter - leadingLines;
    indices.resize(readLines);
    radii.resize(readLines);
    positions.resize(readLines);
}

void findPairs(labelList &indices1, labelList &indices2, labelPairList &pairs)
{
    // remove all entries from first list if they are not present in second list
    // this assumes ordered entries

    pairs.clear();
    pairs.setSize(maxNumParticles);
    label pairCounter = 0;

    if (indices2.size() == 0) return;

    for (label i=0;i<indices1.size();i++)
    {
        label j1 = -1;
        label j2 = indices2.size();
        label jmid = 0;
        label index1 = indices1[i];
        while(true)
        {
            jmid = (j1+j2)/2;
            if (indices2[jmid] > index1) j2 = jmid;
            else if (indices2[jmid] < index1) j1 = jmid;
            else
            {
                pairs[pairCounter]=labelPair(i,jmid);
                pairCounter++;
                break;
            }
            if (j2-j1 == 1) break;
        }
    }
    pairs.resize(pairCounter);
    Info << "findPairs: " << pairs.size() << " pairs found." << endl;
}

void findPairsUnordered(labelList &indices1, labelList &indices2, labelPairList &pairs)
{
    // remove all entries from first list if they are not present in second list
    pairs.clear();
    pairs.setSize(maxNumParticles);
    label pairCounter = 0;

    for (label i=0;i<indices1.size();i++)
    {
        for (label j=0;j<indices2.size();j++)
        {
            if (indices1[i] == indices2[j])
            {
                pairs[pairCounter]=labelPair(i,j);
                pairCounter++;
                break;
            }
        }
    }
    pairs.resize(pairCounter);
    Info << "findPairs: " << pairs.size() << " pairs found." << endl;
}

void fillEmptyCells(fvMesh &mesh, label nNeighMin, label maxSearchLayers, scalarList &particleVolInCell, volVectorField &Us, volVectorField& UsDirectedStdDev,scalarList& boundaries, volVectorField &defaultUs, volVectorField &defaultUsDirectedStdDev, bool interpolate, scalar dt)
{
    labelList neighborsWithValues;
    scalar neighborSqrDistance;
    scalar weight;
    scalar weightSum;
    scalarList weights;

    Info << "Filling empty cells." << endl;
    forAll(mesh.C(), cellI)
    {
        if (particleVolInCell[cellI] > minVol) continue;

        vector position = mesh.C()[cellI];
        label outsideBox = 0;
        if (position.x() < boundaries[0] || position.x() > boundaries[1]) outsideBox++;
        if (position.y() < boundaries[2] || position.y() > boundaries[3]) outsideBox++;
        if (position.z() < boundaries[4] || position.z() > boundaries[5]) outsideBox++;

        if (outsideBox > 0 || !interpolate)
        {
            Us[cellI] = defaultUs[cellI]*dt;
            UsDirectedStdDev[cellI] = defaultUsDirectedStdDev[cellI]*dt;
            continue;
        }

        nearestNeighborCells(mesh, cellI, nNeighMin, maxSearchLayers, particleVolInCell, neighborsWithValues);
        weightSum = 0.0;
        weights.clear();
        for (label neighI=0; neighI<neighborsWithValues.size(); neighI++)
        {
            neighborSqrDistance = magSqr(mesh.C()[cellI] - mesh.C()[neighborsWithValues[neighI]]);
            weight = weightFun(neighborSqrDistance);
            weights.append(weight);
            weightSum += weight;
        }
        for (label neighI=0; neighI<neighborsWithValues.size(); neighI++)
        {
            weight = weights[neighI]/weightSum;
            Us[cellI] += weight*Us[neighborsWithValues[neighI]];
            UsDirectedStdDev[cellI] += weight*UsDirectedStdDev[neighborsWithValues[neighI]];
        }

        if (neighborsWithValues.size() == 0)
        {
            Us[cellI] = defaultUs[cellI]*dt;
            UsDirectedStdDev[cellI] = defaultUsDirectedStdDev[cellI]*dt;
        }

        // make sure no particles are placed outside of domain
        // TODO: correct following implementation (meshSearch) and test it
/*
        vector shiftedPosition = position + dt * Us[cellI];
        label cellJ = mesh.findCell(shiftedPosition);
        if (cellJ < 0)
        {
            label cellK = mesh.findNearestCellWalk(shiftedPosition,cellI);
            Us[cellI] = (mesh.C()[cellI] - mesh.C()[cellK]) / dt;
        }
*/
    }
}

void nearestNeighborCells(fvMesh &mesh, label refCell, label nNeighMin, label maxSearchLayers, scalarList &particleVolInCell, labelList &neighborsWithValues)
{
    label numSearchLayers = 0;
    std::set<label> neighbors;
    std::set<label> newNeighbors;
    std::set<label> recentNeighbors;

    neighbors.insert(refCell);
    recentNeighbors.insert(refCell);

    neighborsWithValues.clear();

    while(neighborsWithValues.size() < nNeighMin)
    {
        for (std::set<label>::iterator it=recentNeighbors.begin(); it!=recentNeighbors.end(); ++it)
        {
            labelList adjacent = mesh.cellCells()[*it];
            label adj;
            for (label j=0; j<adjacent.size(); j++)
            {
                adj = adjacent[j];
                std::set<label>::iterator it2 = neighbors.find(adj);
                if (it2 == neighbors.end())
                {
                    newNeighbors.insert(adj);
                    neighbors.insert(adj);
                    if (particleVolInCell[adj] > minVol) neighborsWithValues.append(adj);
                }
            }
        }

        numSearchLayers++;
        if (numSearchLayers > maxSearchLayers && maxSearchLayers > 0) return;

        if (newNeighbors.size() == 0) return;
        recentNeighbors.clear();
        recentNeighbors = newNeighbors;
        newNeighbors.clear();
    }
}

void normalizeFields(scalarList& particleVolInCell, volVectorField& Us, volVectorField & UsDirectedStdDev)
{
    for (label cellJ = 0; cellJ<particleVolInCell.size(); cellJ++)
    {
        if (particleVolInCell[cellJ] > minVol)
        {
            Us[cellJ] /= particleVolInCell[cellJ];
            UsDirectedStdDev[cellJ] /= particleVolInCell[cellJ];
            for (label comp=0;comp<3;comp++)
            {
                UsDirectedStdDev[cellJ].component(comp) -= Us[cellJ].component(comp)*Us[cellJ].component(comp);
                if (UsDirectedStdDev[cellJ].component(comp) > 0) UsDirectedStdDev[cellJ].component(comp) = Foam::sqrt(UsDirectedStdDev[cellJ].component(comp));
            }
        }
    }
}

scalar weightFun(scalar distSqr)
{
    // inverse distance weighting, order 2
    return 1.0/distSqr;
}
// ************************************************************************* //
