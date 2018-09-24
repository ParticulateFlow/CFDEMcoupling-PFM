/*---------------------------------------------------------------------------*\
    CFDEMcoupling - Open Source CFD-DEM coupling

    CFDEMcoupling is part of the CFDEMproject
    www.cfdem.com
                                Christoph Goniva, christoph.goniva@cfdem.com
                                Copyright 2009-2012 JKU Linz
                                Copyright 2012-     DCS Computing GmbH, Linz
-------------------------------------------------------------------------------
License
    This file is part of CFDEMcoupling.

    CFDEMcoupling is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 3 of the License, or (at your
    option) any later version.

    CFDEMcoupling is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with CFDEMcoupling; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Description
    This code is designed to realize coupled CFD-DEM simulations using LIGGGHTS
    and OpenFOAM(R). Note: this code is not part of OpenFOAM(R) (see DISCLAIMER).
\*---------------------------------------------------------------------------*/

#include "error.H"

#include "particleProbe.H"
#include "addToRunTimeSelectionTable.H"
#include <mpi.h>
#include "IOmanip.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(particleProbe, 0);

addToRunTimeSelectionTable
(
    probeModel,
    particleProbe,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
particleProbe::particleProbe
(
    const dictionary& dict,
    cfdemCloud& sm,
    const word& typeName,
    const char* logFileName
)
:
    probeModel(dict,sm,typeName,logFileName),
    propsDict_(dict.subDict(typeName + "Props")),
    name_(typeName),
    particleCloud_(sm),
    verbose_(propsDict_.found("verbose")),
    verboseToFile_(propsDict_.found("verboseToFile")),
    writePrecision_(propsDict_.lookupOrDefault<int>("writePrecision", 3)),
    dirName_("particleProbes"),
    rank_(-1),
    sPtr(NULL),
    printEvery_(propsDict_.lookupOrDefault<int>("printEvery", 1)),
    sampleAll_(propsDict_.found("sampleAll")),
    probeDebug_(propsDict_.found("probeDebug")),
    includePosition_(propsDict_.found("includePosition")),
    particleIDsToSample_(propsDict_.lookup("particleIDsToSample")),
    itemCounter_(0),
    currItemId_(0),
    printCounter_(0),
    printNow_(false)
{
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

particleProbe::~particleProbe()
{
    forAll(sPtrList_, i)
        delete sPtrList_[i];
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void particleProbe::setOutputFile(const word& logFileName)
{
    if (itemCounter_ > 0 && verboseToFile_)
    {
        bool foundFile = false;
        forAll(itemsToSample_, i)
        {
            if (itemsToSample_[i] == logFileName)
            {
                probeIndex_ = i;
                foundFile = true;
            }
        }

        if(!foundFile)
             FatalError << "particleProbe::setOutputFile for logFileName " << logFileName << " : " << "File not found" << abort(FatalError);
        currItemId_ = probeIndex_ + 1;
        setCounter();
    }
}


void particleProbe::initialize(const word& modelName, const word& logFileName)
{
    //update the list of items to be sampled
    ++itemCounter_;
    itemsToSample_.append(logFileName);

    // init environment
    //propsDict_ = particleCloud_.couplingProperties().subDict(typeName + "Props");
    name_ = modelName;

    if (verboseToFile_)
    {
        Info << "Will sample these particle IDs: " << particleIDsToSample_ << " every " << printEvery_ << endl;

        //initialize the output files
        MPI_Comm_rank(MPI_COMM_WORLD, &rank_);

        //open a separate file for each processor
        char* filecurrent_ = new char[logFileName.length() + 1 + 4 + 1]; //reserve 4 chars for processor name
        sprintf(filecurrent_,"%s.%d", logFileName.c_str(), rank_);

        Info << "particleProbe for model " << name_ << " will write to file " << filecurrent_ << endl;

        //generate the file streams
        fileName probeSubDir = dirName_;
        if (particleCloud_.mesh().name() != polyMesh::defaultRegion)
        {
            probeSubDir = probeSubDir/particleCloud_.mesh().name();
        }
        probeSubDir = probeSubDir/particleCloud_.mesh().time().timeName();

        fileName probeDir_;
        if (Pstream::parRun())
        {
            // Put in undecomposed case
            // (Note: gives problems for distributed data running)
            probeDir_ = particleCloud_.mesh().time().path()/".."/probeSubDir;
        }
        else
        {
            probeDir_ = particleCloud_.mesh().time().path()/probeSubDir;
        }

        //manage files and OFstreams
        mkDir(probeDir_);
        sPtr = new OFstream(probeDir_/filecurrent_);
        sPtrList_.append(sPtr);

        delete [] filecurrent_;
        //Clear the containers for the fields to be probed
        scalarFields_.clear();
        vectorFields_.clear();
    }
}


void particleProbe::writeHeader() const
{
    if (verboseToFile_)
    {
        *sPtr << "#processor: " << rank_ << endl;
        *sPtr << "#index   time   " << "   ";
        *sPtr << "||  vectorData:  " << "   ";

        forAll(vectorFields_, iter)
        {
            if (!probeDebug_ && iter > 0) break;
            *sPtr << vectorFields_(iter) << "   ";
        }

        if (probeDebug_)
        {
            *sPtr << "||   scalarData:  "  << "   ";
            forAll(scalarFields_, iter)
            {
                *sPtr << scalarFields_(iter)  << "   ";
            }
        }

        if (includePosition_) *sPtr << " ||  position" << endl;
        else *sPtr << endl;
    }
}


void particleProbe::writeProbe(int index, Field<scalar> sValues, Field<vector> vValues)
{
    if (printNow_ && verboseToFile_ && checkIDForPrint(index))
    {
        sPtr = sPtrList_[probeIndex_]; //set the pointer to the output file from list

        //index and time
        *sPtr << setprecision(IOstream::defaultPrecision()+7);
        *sPtr << index  << tab << particleCloud_.mesh().time().value()  << "   ";
        *sPtr << "||   ";

        //vectorFields
        *sPtr << setprecision(writePrecision_);
        forAll(vValues, iter)
        {
            // if(!probeDebug_ && iter>0) break;
            *sPtr << vValues[iter][0] << "   ";
            *sPtr << vValues[iter][1] << "   ";
            *sPtr << vValues[iter][2] << "   ";
        }

        //scalarFields
        if(probeDebug_)
        {
            *sPtr << "||   ";
            forAll(sValues, iter)
            {
                *sPtr << sValues[iter] << "   ";
            }
        }

        if(includePosition_)
        {
            *sPtr << "||   ";
            *sPtr << particleCloud_.position(index)[0] << "   "
                  << particleCloud_.position(index)[1] << "   "
                  << particleCloud_.position(index)[2]
                  << endl;
        }
        else
        {
            *sPtr << endl;
        }
    }
}


bool particleProbe::checkIDForPrint(int index) const
{
    if(sampleAll_)
    {
        return true;
    }
    else
    {
        forAll(particleIDsToSample_, iSample)
        {
            if (index == particleIDsToSample_[iSample]) return true;
        }
    }
    return false;
}

void particleProbe::setCounter()
{
    //reset or increment counter for printing to file
    //Do only if called by first item in the list of items!
    if (currItemId_ == 1)
    {
        ++printCounter_;

        if (printCounter_ >= printEvery_)
        {
            printCounter_ = 0;
            printNow_ = true;
        }
        else printNow_ = false;
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
