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
void interpolateCellValues(fvMesh &, label , labelList &, volVectorField &, volVectorField &);
void nearestNeighborCells(fvMesh &, label, label, labelList &, labelList &);
void readDump(std::string, labelList &, vectorList &);
scalar weightFun(scalar);

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

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
    label dumpIndexIncrement(readLabel(displacementProperties.lookup("dumpIndexIncrement")));
    label nNeighMin(readLabel(displacementProperties.lookup("nNeighMin")));
    scalar timePerStep(readScalar(displacementProperties.lookup("timePerStep")));
    scalar startTime(readScalar(displacementProperties.lookup("startTime")));
    std::string filepath=string(displacementProperties.lookup("filepath"));
    std::string fileext=string(displacementProperties.lookupOrDefault<string>("fileextension",""));
    bool fillEmptyCells=bool(displacementProperties.lookupOrDefault<bool>("fillEmptyCells",true));

    label dumpIndex1 = dumpIndexStart;
    label dumpIndex2 = dumpIndex1 + dumpIndexIncrement;

    volVectorField Us
    (
        IOobject
        (
            "Us",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedVector("zero", dimensionSet(0,1,-1,0,0), vector::zero)
    );

    volVectorField UsDirectedVariance
    (
        IOobject
        (
            "UsDirectedVariance",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedVector("zero", dimensionSet(0,1,-1,0,0), vector::zero)
    );

    scalar currTime=startTime;
    label timeIndex=0;

    while(true)
    {
        runTime.setTime(currTime,timeIndex);
        // read dump files and check which particle indices are present in both
        labelList indices1, indices2;
        vectorList positions1, positions2;

        std::stringstream ss;
        ss << filepath << dumpIndex1 << fileext;
        std::string filename1 = ss.str();
        ss.str("");
        ss << filepath << dumpIndex2 << fileext;
        std::string filename2 = ss.str();

        if (access( filename1.c_str(), F_OK ) == -1 || access( filename2.c_str(), F_OK ) == -1) break;

        Info << "\nReading" << endl;
        Info << "\t" << filename1 << endl;
        Info << "\t" << filename2 << endl;
        Info << "corresponding to time = " << currTime << "." << endl;

        readDump(filename1, indices1, positions1);
        readDump(filename2, indices2, positions2);

        labelPairList pairs;
        findPairs(indices1,indices2,pairs);

        // average particle displacements and their variance
        Info << "Binning particle displacements on mesh." << endl;
        labelList particlesInCell(mesh.nCells(), 0);
        vector position, displacement;
        label line1, line2;
        label cellI;

        Us *= 0.0;
        UsDirectedVariance *= 0.0;
        for (label partI = 0; partI < pairs.size(); partI++)
        {
            line1 = pairs[partI].first();
            line2 = pairs[partI].second();
            position = positions1[line1];
            displacement = positions2[line2] - positions1[line1];
            cellI = mesh.findCell(position);
            if (cellI < 0)  continue;
            particlesInCell[cellI] += 1;
            Us[cellI] += displacement;

            for (label comp=0;comp<3;comp++)
            {
                UsDirectedVariance[cellI].component(comp) += displacement.component(comp)*displacement.component(comp);
            }
        }


        for (label cellJ = 0; cellJ<particlesInCell.size(); cellJ++)
        {
            if (particlesInCell[cellJ] > 0)
            {
                Us[cellJ] /= particlesInCell[cellJ];
                UsDirectedVariance[cellJ] /= particlesInCell[cellJ];
                for (label comp=0;comp<3;comp++)
                {
                    UsDirectedVariance[cellJ].component(comp) -= Us[cellJ].component(comp)*Us[cellJ].component(comp);
                    if (UsDirectedVariance[cellJ].component(comp) > 0) UsDirectedVariance[cellJ].component(comp) = Foam::sqrt(UsDirectedVariance[cellJ].component(comp));
                }
            }
        }


        // interpolate values for empty cells
        if (fillEmptyCells)
        {
            Info << "Interpolating empty cells." << endl;
            interpolateCellValues(mesh,nNeighMin,particlesInCell,Us,UsDirectedVariance);
        }

        dumpIndex1 = dumpIndex2;
        dumpIndex2 = dumpIndex1 + dumpIndexIncrement;

        Us /= timePerStep;
        UsDirectedVariance /= timePerStep;
        Us.write();
        UsDirectedVariance.write();
        currTime += timePerStep;
        timeIndex++;
    }
    return 0;
}

void readDump(std::string filename, labelList &indices, vectorList &positions)
{
    #include <fstream>

    const label leadingLines = 9;
    label lineCounter = 0;
    label partIndex;
    scalar x, y, z;

    indices.clear();
    positions.clear();

    std::ifstream file(filename);
    std::string str; 
    while (std::getline(file, str))
    {
        if (lineCounter >= leadingLines)
        {
            sscanf(str.c_str(), "%d %lf %lf %lf", &partIndex, &x, &y, &z);
            indices.append(partIndex);
            positions.append(vector(x,y,z));
        }
        lineCounter++;
    }
}

void findPairs(labelList &indices1, labelList &indices2, labelPairList &pairs)
{
    // remove all entries from first list if they are not present in second list
    // this assumes ordered entries

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
                pairs.append(labelPair(i,jmid));
                break;
            }
            if (j2-j1 == 1) break;
        }
    }
    Info << "findPairs: " << pairs.size() << " pairs found." << endl;
}

void findPairsUnordered(labelList &indices1, labelList &indices2, labelPairList &pairs)
{
    // remove all entries from first list if they are not present in second list

    for (label i=0;i<indices1.size();i++)
    {
        for (label j=0;j<indices2.size();j++)
        {
            if (indices1[i] == indices2[j])
            {
                pairs.append(labelPair(i,j));
                break;
            }
        }
    }
    Info << "findPairs: " << pairs.size() << " pairs found." << endl;
}

void interpolateCellValues(fvMesh &mesh, label nNeighMin, labelList &particlesInCell, volVectorField &Us, volVectorField& UsDirectedVariance)
{
    labelList neighborsWithValues;
    scalar neighborSqrDistance;
    scalar weight;
    scalar weightSum;
    scalarList weights;
    forAll(mesh.C(), cellI)
    {
        if (particlesInCell[cellI] > 0) continue;
        nearestNeighborCells(mesh, cellI, nNeighMin, particlesInCell, neighborsWithValues);
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
            UsDirectedVariance[cellI] += weight*UsDirectedVariance[neighborsWithValues[neighI]];
        }
    }
}

void nearestNeighborCells(fvMesh &mesh, label refCell, label nNeighMin, labelList &particlesInCell, labelList &neighborsWithValues)
{
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
                    if (particlesInCell[adj] > 0) neighborsWithValues.append(adj);
                }
            }
        }

        if (newNeighbors.size() == 0) return;
        recentNeighbors.clear();
        recentNeighbors = newNeighbors;
        newNeighbors.clear();
    }
}

scalar weightFun(scalar distSqr)
{
    // inverse distance weighting, order 2
    return 1.0/distSqr;
}
// ************************************************************************* //
