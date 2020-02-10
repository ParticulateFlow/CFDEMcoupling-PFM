/*---------------------------------------------------------------------------*\
    CFDEMcoupling academic - Open Source CFD-DEM coupling

    Contributing authors:
    Thomas Lichtenegger, Gerhard Holzinger
    Copyright (C) 2015- Johannes Kepler University, Linz
-------------------------------------------------------------------------------
License
    This file is part of CFDEMcoupling academic.

    CFDEMcoupling academic is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    CFDEMcoupling academic is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with CFDEMcoupling academic.  If not, see <http://www.gnu.org/licenses/>.
\*---------------------------------------------------------------------------*/

#include "error.H"
#include "Random.H"
#include "gerhardsRecModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(gerhardsRecModel, 0);

addToRunTimeSelectionTable
(
    recModel,
    gerhardsRecModel,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
gerhardsRecModel::gerhardsRecModel
(
    const dictionary& dict,
    recBase& base
)
:
    recModel(dict,base),
    propsDict_(dict.subDict(typeName + "Props")),
    dataBaseName_(propsDict_.lookupOrDefault<word>("dataBase", word("dataBase"))),
    recTime(fileName(dataBaseName_), "", "../system", "../constant", false),
    timeDirs(recTime.times()),
    skipZero_(propsDict_.lookupOrDefault<Switch>("skipZero", Switch(false))),
    checkSkipZero_(checkSkipZero()),
    numRecFields_(skipZero_ ? label(timeDirs.size())-1 : label(timeDirs.size())),
    verboseVerbose_(propsDict_.lookupOrDefault<Switch>("verboseVerbose", Switch(false))),
    verboseVacate_(propsDict_.lookupOrDefault<Switch>("verboseVacate", Switch(false))),
    recurrenceMatrix_(numRecFields_, scalar(-1.0)),
    timeIndexList_(numRecFields_-1),
    revTimeIndexList_(numRecFields_-1),
    timeValueList_(numRecFields_-1),
    contTimeIndex(0),
    lowerSeqLim_(max(1, label(numRecFields_/20))),
    upperSeqLim_(label(numRecFields_/5)),
    numDataBaseFields_(propsDict_.lookupOrDefault<label>("numDataBaseFields",numRecFields_)),
    storageIndex_(numRecFields_),
    storageUsageList_(numRecFields_),
    volScalarFieldList_(volScalarFieldNames_.size()),
    volVectorFieldList_(volVectorFieldNames_.size()),
    surfaceScalarFieldList_(surfaceScalarFieldNames_.size()),
    nrOfReadsFromDisk_(0)
{
    /*
        Sanity checks

        "Sanity is just another weird form of Madness"
                                                --- Darkwell
    */
    if (numDataBaseFields_ < 2)
    {
        FatalError << "Number of fields in dataBase specified in "
            << propsDict_.name() << " is smaller than 2!"
            << abort(FatalError);
    }

    if (numDataBaseFields_ > numRecFields_)
    {
        FatalError << "Number of fields in dataBase specified in "
            << propsDict_.name() << " is larger than number of snapshots!"
            << abort(FatalError);
    }

    if (dataBaseName_ == "")
    {
        FatalError << "Empty dataBase path provided"
            << abort(FatalError);
    }

    if (numRecFields_ == 0)
    {
        FatalError << "No snapshots in dataBase! The dataBase is empty."
            << abort(FatalError);
    }

    if (upperSeqLim_ == 0)
    {
        FatalError << "Bad sequence limits! The dataBase is most probably too small."
            << abort(FatalError);
    }

    // initialize data structures
    forAll(storageUsageList_, i)
    {
        storageUsageList_[i] = storageUsageList_.size()+1;
    }

    if (verbose_)
    {
        // be informative on properties of the "recTime" Time-object
        Info << "recTime.rootPath() " << recTime.rootPath() << endl;
        Info << "recTime.caseName() " << recTime.caseName() << endl;
        Info << "recTime.path() " << recTime.path() << endl;
        Info << "recTime.timePath() " << recTime.timePath() << endl;
        Info << "recTime.timeName() " << recTime.timeName() << endl;
        Info << "timeDirs " << timeDirs << endl;
        Info << "ignoring 0 directory: " << skipZero_ << endl;
        Info << "number of snapshots: " << numRecFields_ << endl;
        Info << "number of dataBase slots: " << numDataBaseFields_ << endl;
    }
    readTimeSeries();

    recTimeStep_ = checkTimeStep();
    totRecSteps_ = 1 + static_cast<label>( (endTime_-startTime_) / recTimeStep_ );

    for(int i=0; i<volScalarFieldNames_.size(); i++)
    {
        volScalarFieldList_[i].setSize(numDataBaseFields_);
    }

    for(int i=0; i<volVectorFieldNames_.size(); i++)
    {
        volVectorFieldList_[i].setSize(numDataBaseFields_);
    }

    for(int i=0; i<surfaceScalarFieldNames_.size(); i++)
    {
        surfaceScalarFieldList_[i].setSize(numDataBaseFields_);
    }

    //   setRecFields();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

gerhardsRecModel::~gerhardsRecModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
scalar gerhardsRecModel::checkTimeStep()
{
   // check time step of provided data
    scalar dtCur(0.0);
    scalar dtOld(0.0);

    if (verbose_)
    {
        Info << "timeValueList : " << timeValueList_ << endl;
    }

    forAll(timeValueList_, i)
    {
        // skip zero
        if (skipZero_ and timeDirs[i].value() == 0)
        {
            if (verbose_)
            {
                Info << " ... skipping 0 in checkTimeStep()" << endl;
            }
            continue;
        }

        // compute time step
        if (timeDirs[i].value() == timeDirs.last().value())
        {
            if (verbose_)
            {
                Info << ".. leaving loop at " << timeDirs[i] << endl;
            }
            // leave loop
            break;
        }

        if (verbose_)
        {
            Info << "timeDirs.fcIndex(i)].value(),  timeDirs[i].value() : "
                    << timeDirs[timeDirs.fcIndex(i)].value() << "   " << timeDirs[i].value()
                    << endl;
        }

        // the documentation is in the code ;-)
        //    fcIndex() - return forward circular index, i.e. the next index
        dtCur = timeDirs[timeDirs.fcIndex(i)].value() - timeDirs[i].value();

        if (dtOld < SMALL)
        {
            dtOld = dtCur;
        }

        if (abs(dtOld - dtCur) > SMALL)
        {
            Info << "dtCur, dtOld = " << dtCur << "   " << dtOld << endl;
            FatalError << "    in setting up data" << nl
                            << "    non-constant time-step of provided simulation data"
                            << abort(FatalError);
        }
    }

    // set deltaT
#if OPENFOAM_VERSION_MAJOR < 6
    recTime.setDeltaT(dtCur, false);
#else
    recTime.setDeltaT(dtCur);
#endif

    if (verbose_)
    {
        Info << "Setting deltaRecT to " << dtCur << endl;
        Info << "Actual recTime.deltaT = " << recTime.deltaTValue() << endl;
        Info << "Actual runTime.deltaT = " << timeStep_ << endl;
    }

    return dtCur;
}

void gerhardsRecModel::readFieldSeries()
{
    Info << "Checking all " << numRecFields() << " fields\n" << endl;

    label size = timeDirs.size();
    label counter = 0;
    label percentage = 0;

    for (instantList::iterator it=timeDirs.begin(); it != timeDirs.end(); ++it)
    {
        if(counter >= 0.1 * percentage * size)
        {
            Info << "\t" << 10 * percentage << " \% done" << endl;
            percentage++;
        }
        counter++;

        // set time
        recTime.setTime(*it, it->value());

        // skip zero
        if (skipZero_ and recTime.timeName() == "0")
        {
            if (verbose_)
            {
                Info << " ... skipping 0 in readFieldSeries()" << endl;
            }

            continue;
        }

        // skip constant
        if (recTime.timeName() == "constant")
        {
            continue;
        }

        if (verbose_)
        {
            Info << "Checking fields at t = " << recTime.timeName() << endl;
        }

        for (int i=0; i<volScalarFieldNames_.size(); i++)
        {
            IOobject header
            (
                volScalarFieldNames_[i],
                recTime.timePath(),
                base_.mesh(),
                IOobject::MUST_READ
            );

            // Check if volScalarFieldNames_[i] is a valid field
            /*
                COMPAT: replace the method typeHeaderOk<volScalarField>(true) with headerOk()
                for OpenFOAM versions prior to OpenFOAM-5.0
                Do this the other way around for OpenFOAM-5.0 and potentially later versions
            */
#if OPENFOAM_VERSION_MAJOR < 5
            if (! header.headerOk())
#else
            if (! header.typeHeaderOk<volScalarField>(true))
#endif
            {
                FatalError
                << "Field " << volScalarFieldNames_[i] << " not found"
                << exit(FatalError);
            }
        }

        for (int i=0; i<volVectorFieldNames_.size(); i++)
        {
            IOobject header
            (
                volVectorFieldNames_[i],
                recTime.timePath(),
                base_.mesh(),
                IOobject::MUST_READ
            );

            // Check if volVectorFieldNames_[i] is a valid field
            /*
                COMPAT: replace the method typeHeaderOk<volScalarField>(true) with headerOk()
                for OpenFOAM versions prior to OpenFOAM-5.0
                Do this the other way around for OpenFOAM-5.0 and potentially later versions
            */
#if OPENFOAM_VERSION_MAJOR < 5
            if (! header.headerOk())
#else
            if (! header.typeHeaderOk<volVectorField>(true))
#endif
            {
                FatalError
                << "Field " << volVectorFieldNames_[i] << " not found"
                << exit(FatalError);
            }
        }

        for (int i=0; i<surfaceScalarFieldNames_.size(); i++)
        {
            IOobject header
            (
                surfaceScalarFieldNames_[i],
                recTime.timePath(),
                base_.mesh(),
                IOobject::MUST_READ
            );

            // Check if surfaceScalarFieldNames_[i] is a valid field
            /*
                COMPAT: replace the method typeHeaderOk<volScalarField>(true) with headerOk()
                for OpenFOAM versions prior to OpenFOAM-5.0
                Do this the other way around for OpenFOAM-5.0 and potentially later versions
            */
#if OPENFOAM_VERSION_MAJOR < 5
            if (! header.headerOk())
#else
            if (! header.typeHeaderOk<surfaceScalarField>(true))
#endif
            {
                FatalError
                << "Field " << surfaceScalarFieldNames_[i] << " not found"
                << exit(FatalError);
            }
        }

        // reset storageIndex_ to -1
        storageIndex_[timeIndexList_(recTime.timeName())] = -1;
    }

    Info << "Checking fields done" << nl << endl;


    Info << "Reading fields for " << numDataBaseFields() << " dataBase slots\n" << endl;

    size = numDataBaseFields();
    counter = 0;
    percentage = 0;


    label fieldCounter(0);
    for (instantList::iterator it=timeDirs.begin(); it != timeDirs.end(); ++it)
    {
        if(counter >= 0.1 * percentage * size)
        {
            Info << "\t" << 10 * percentage << " \% done" << endl;
            percentage++;
        }
        counter++;

        // set time
        recTime.setTime(*it, it->value());

        // skip constant
        if (recTime.timeName() == "constant")
        {
            continue;
        }

        // skip zero
        if (skipZero_ and recTime.timeName() == "0")
        {
            if (verbose_)
            {
                Info << " ... skipping 0 in readTimeSeries()" << endl;
            }

            continue;
        }

        if (verbose_)
        {
            Info << "Reading fields at t = " << recTime.timeName() << endl;
        }

        nrOfReadsFromDisk_++;

        for (int i=0; i<volScalarFieldNames_.size(); i++)
        {
            volScalarFieldList_[i].set
            (
                fieldCounter,
                new volScalarField
                (
                    IOobject
                    (
                        volScalarFieldNames_[i],
                        recTime.timePath(),
                        base_.mesh(),
                        IOobject::MUST_READ,
                        IOobject::NO_WRITE
                    ),
                    base_.mesh()
                )
            );
        }

        for (int i=0; i<volVectorFieldNames_.size(); i++)
        {
            volVectorFieldList_[i].set
            (
                fieldCounter,
                new volVectorField
                (
                    IOobject
                    (
                        volVectorFieldNames_[i],
                        recTime.timePath(),
                        base_.mesh(),
                        IOobject::MUST_READ,
                        IOobject::NO_WRITE
                    ),
                    base_.mesh()
                )
            );
        }

        for (int i=0; i<surfaceScalarFieldNames_.size(); i++)
        {
            surfaceScalarFieldList_[i].set
            (
                fieldCounter,
                new surfaceScalarField
                (
                    IOobject
                    (
                        surfaceScalarFieldNames_[i],
                        recTime.timePath(),
                        base_.mesh(),
                        IOobject::MUST_READ,
                        IOobject::NO_WRITE
                    ),
                    base_.mesh()
                )
            );
        }

        // set storageIndex_
        //  storageIndex_[index] translates the index of the snapshot
        //  to an appropriate index of the dataBase
        storageIndex_[timeIndexList_(recTime.timeName())] = fieldCounter;


        fieldCounter++;

        if (fieldCounter >= numDataBaseFields_)
        {
            break;
        }
    }

    Info << "Reading fields done" << endl;
}


void gerhardsRecModel::readTimeSeries()
{
    bool firsttime = true;
    // fill the data structure for the time indices
    for (instantList::iterator it=timeDirs.begin(); it != timeDirs.end(); ++it)
    {
        // set run-time
        recTime.setTime(*it, it->value());


        // skip constant
        if (recTime.timeName() == "constant")
        {
            continue;
        }

        // skip zero
        if (skipZero_ and recTime.timeName() == "0")
        {
            if (verbose_)
            {
                Info << " ... skipping 0 in readTimeSeries()" << endl;
            }

            continue;
        }

        if (firsttime)
        {
            firsttime = false;
            recStartTime_ = recTime.value();
        }
        recEndTime_ = recTime.value();

        // insert the time name into the hash-table with a continuous second index
        timeIndexList_.insert(recTime.timeName(), contTimeIndex);
        revTimeIndexList_.insert(contTimeIndex, recTime.timeName());


        if (verbose_)
        {
            Info << "current time " << recTime.timeName() << endl;
            Info << "insert " << recTime.timeName() << " , " << contTimeIndex << endl;
        }


        // insert the time value
        timeValueList_.insert(contTimeIndex, recTime.timeOutputValue());

        // increment continuousIndex
        contTimeIndex++;

        if (verbose_)
        {
            Info << "contTimeIndex " << contTimeIndex << endl;
        }
    }

    if (verbose_)
    {
        Info << endl;
        Info << "Found " << label(timeDirs.size()) << " time folders" << endl;
        Info << "Found " << label(timeIndexList_.size()) << " time steps" << endl;
        Info << "database start time = " << recStartTime_ << endl;
        Info << "database end time = " << recEndTime_ << endl;
    }
}



void gerhardsRecModel::fetchStateForDataBase(label index)
{
    if (verboseVerbose_)
    {
        Info << " ... calling fetchStateForDataBase(" << index << ")" << endl;
        Info << "     storageUsageList_ : " << storageUsageList_ << endl;
        Info << "     storageIndex_ : " << storageIndex_ << endl;
    }

    // decide which state to remove from the data base
    //  find maximum element in storageUsageList_, as it is the index of the least used state
    label oldIndex(0);
    label maxVal(0);

    // find maximum value in storageUsageList_
    forAll(storageUsageList_, i)
    {
        if (storageUsageList_[i] > maxVal and storageIndex_[i] > -1)
        {
            oldIndex = i;
            maxVal = storageUsageList_[i];
        }
    }
    if (verboseVerbose_ or verboseVacate_)
    {
        Info << " --> vacate dataBase  : " << oldIndex << endl;
    }

    // read state for index from disk and put it into the data base
    readNewSnapshot(index, oldIndex);


    // reset storageIndex
    //  the slot at oldIndex in the dataBase is replaced by the
    //  snapshot for index
    storageIndex_[index] = storageIndex_[oldIndex];
    storageIndex_[oldIndex] = -1;

    // update storageIndex
    updateStorageUsageList(index);

    if (verboseVerbose_)
    {
        Info << " ... finished fetchStateForDataBase(" << index << ")" << endl;
        Info << "     storageUsageList_ : " << storageUsageList_ << endl;
        Info << "     storageIndex_ : " << storageIndex_ << endl;
    }
}


void gerhardsRecModel::updateStorageUsageList(label index)
{
    if (verboseVerbose_)
    {
        Info << " ... calling updateStorageUsageList(" << index << ")" << endl;
        Info << "     storageUsageList_ : " << storageUsageList_ << endl;
        Info << "   storageIndex_ : " << storageIndex_ << endl;
    }

    label oldValue(storageUsageList_[index]);

    /*
        increment (make less recent) all elements more recent (smaller)
        than the element at slot index
    */
    forAll(storageUsageList_, i)
    {
        if (storageUsageList_[i] < oldValue)
        {
            storageUsageList_[i]++;
        }
    }

    // storageUsageList_[index] is the most recently used element in dataBase
    storageUsageList_[index] = 1;

    if (verboseVerbose_)
    {
        Info << " ... finished updateStorageUsageList(" << index << ")" << endl;
        Info << "     storageUsageList_ : " << storageUsageList_ << endl;
        Info << "     storageIndex_ : " << storageIndex_ << nl << endl;
    }
}


void gerhardsRecModel::readNewSnapshot(label index, label oldIndex)
{
    if (verboseVerbose_)
    {
        Info << "      ... calling readNewSnapshot(" << index << ", " << oldIndex << ")" << endl;
        Info << "          index = " << index << "; oldIndex = " << oldIndex << endl;
        Info << "          storageUsageList_ : " << storageUsageList_ << endl;
        Info << "          storageIndex_ : " << storageIndex_ << nl << endl;
    }

    /*
        Replace all fields in the slot "oldIndex" of the dataBase with fields
        for the snapshot at "index"
    */

    // vacate slot
    for (int i=0; i<volScalarFieldNames_.size(); i++)
    {
        volScalarFieldList_[i][storageIndex_[oldIndex]].clear();
    }
    for (int i=0; i<volVectorFieldNames_.size(); i++)
    {
        volVectorFieldList_[i][storageIndex_[oldIndex]].clear();
    }
    for (int i=0; i<surfaceScalarFieldNames_.size(); i++)
    {
        surfaceScalarFieldList_[i][storageIndex_[oldIndex]].clear();
    }


    // build path to read from
    word readPath(dataBaseName_+"/"+timeDirs[index].name());

    if (dataBaseName_ == ".")
    {
        readPath = timeDirs[index].name();
    }

    nrOfReadsFromDisk_++;

    // repopulate slot
    for (int i=0; i<volScalarFieldNames_.size(); i++)
    {
        if (verboseVerbose_)
        {
            Info << " .... reading " << volScalarFieldNames_[i] << " at time "
                << timeDirs[index] << endl;
        }

        volScalarFieldList_[i].set
        (
            storageIndex_[oldIndex],
            new volScalarField
            (
                IOobject
                (
                    volScalarFieldNames_[i],
                    readPath,
                    base_.mesh(),
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                ),
                base_.mesh()
            )
        );
    }

    for (int i=0; i<volVectorFieldNames_.size(); i++)
    {
        volVectorFieldList_[i].set
        (
            storageIndex_[oldIndex],
            new volVectorField
            (
                IOobject
                (
                    volVectorFieldNames_[i],
                    readPath,
                    base_.mesh(),
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                ),
                base_.mesh()
            )
        );
    }

    for (int i=0; i<surfaceScalarFieldNames_.size(); i++)
    {
        surfaceScalarFieldList_[i].set
        (
            storageIndex_[oldIndex],
            new surfaceScalarField
            (
                IOobject
                (
                    surfaceScalarFieldNames_[i],
                    readPath,
                    base_.mesh(),
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                ),
                base_.mesh()
            )
        );
    }

}




void gerhardsRecModel::exportVolScalarField(word fieldname, volScalarField& field)
{
    field = exportVolScalarField(fieldname, virtualTimeIndex);
}


void gerhardsRecModel::exportVolVectorField(word fieldname, volVectorField& field)
{
    field = exportVolVectorField(fieldname, virtualTimeIndex);
}


void gerhardsRecModel::exportSurfaceScalarField(word fieldname, surfaceScalarField& field)
{
    field = exportSurfaceScalarField(fieldname, virtualTimeIndex);
}




const volScalarField& gerhardsRecModel::exportVolScalarField(word fieldname, label index)
{
    const label fieldI = getVolScalarFieldIndex(fieldname, index);

    if (verboseVerbose_)
    {
        Info << nl << " ... calling exportVolScalarField(" << fieldname << ", " << index << ")" << endl;
    }

    if (numDataBaseFields_ == numRecFields_)
    {
        return volScalarFieldList_[fieldI][index];
    }

    // check whether field with index is in the dataBase
    //  translate index to storageIndex
    if (storageIndex_[index] < 0)
    {
        fetchStateForDataBase(index);
    }
    else
    {
        updateStorageUsageList(index);
    }

    return volScalarFieldList_[fieldI][storageIndex_[index]];
}



const volVectorField& gerhardsRecModel::exportVolVectorField(word fieldname, label index)
{
    const label fieldI = getVolVectorFieldIndex(fieldname, index);

    if (verboseVerbose_)
    {
        Info << nl << " ... calling exportVolVectorField(" << fieldname << ", " << index << ")" << endl;
    }

    if (numDataBaseFields_ == numRecFields_)
    {
        return volVectorFieldList_[fieldI][index];
    }


    // check whether field with index is in the dataBase
    //  translate intex to storageIndex
    if (storageIndex_[index] < 0)
    {
        fetchStateForDataBase(index);
    }
    else
    {
        updateStorageUsageList(index);
    }

    return volVectorFieldList_[fieldI][storageIndex_[index]];
}

const surfaceScalarField& gerhardsRecModel::exportSurfaceScalarField(word fieldname, label index)
{
    const label fieldI = getSurfaceScalarFieldIndex(fieldname, index);

    if (numDataBaseFields_ == numRecFields_)
    {
        return surfaceScalarFieldList_[fieldI][index];
    }

    // check whether field with index is in the dataBase
    //  translate intex to storageIndex
    if (storageIndex_[index] < 0)
    {
        fetchStateForDataBase(index);
    }
    else
    {
        updateStorageUsageList(index);
    }

    return surfaceScalarFieldList_[fieldI][storageIndex_[index]];
}


PtrList<volScalarField>& gerhardsRecModel::exportVolScalarFieldList(word fieldname)
{
    const label fieldI = getVolScalarFieldIndex(fieldname);

    return volScalarFieldList_[fieldI];
}

PtrList<volVectorField>& gerhardsRecModel::exportVolVectorFieldList(word fieldname)
{
    const label fieldI = getVolVectorFieldIndex(fieldname);

    return volVectorFieldList_[fieldI];
}

PtrList<surfaceScalarField>& gerhardsRecModel::exportSurfaceScalarFieldList(word fieldname)
{
    const label fieldI = getSurfaceScalarFieldIndex(fieldname);

    return surfaceScalarFieldList_[fieldI];
}



SymmetricSquareMatrix<scalar>& gerhardsRecModel::recurrenceMatrix()
{
    return recurrenceMatrix_;
}


const HashTable<label,word>& gerhardsRecModel::timeIndexList() const
{
    return timeIndexList_;
}


label gerhardsRecModel::lowerSeqLim() const
{
    return lowerSeqLim_;
}

label gerhardsRecModel::upperSeqLim() const
{
    return upperSeqLim_;
}


label gerhardsRecModel::numRecFields() const
{
    return numRecFields_;
}

label gerhardsRecModel::numDataBaseFields() const
{
    return numDataBaseFields_;
}

void gerhardsRecModel::updateRecFields()
{
    virtualTimeIndex=virtualTimeIndexNext;
    virtualTimeIndexNext++;
    if (virtualTimeIndexNext>sequenceEnd)
    {
        virtualTimeIndexListPos_++;
        sequenceStart=virtualTimeIndexList_[virtualTimeIndexListPos_].first();
        sequenceEnd=virtualTimeIndexList_[virtualTimeIndexListPos_].second();
        virtualTimeIndexNext=sequenceStart;
    }

    if (verbose_)
    {
        Info << "\nUpdating virtual time index to " << virtualTimeIndex << ".\n" << endl;
    }
}


void gerhardsRecModel::writeRecMatrix() const
{
    if (writeRecMat_)
    {
        OFstream matrixFile(recMatName_);
        matrixFile << recurrenceMatrix_;

        Info << nl << "Nr. of reads from disk : " << nrOfReadsFromDisk_
            << "; compared to a theoretical minimum of "
            << numRecFields() << " reads." << endl;
    }
}




Switch gerhardsRecModel::checkSkipZero()
{
    if (skipZero_)
    {
        bool foundZero(false);

        forAll(timeDirs, i)
        {
            if (timeDirs[i].value() == 0)
            {
                foundZero = true;
            }
        }

        if (! foundZero)
        {
            skipZero_ = false;
        }
    }

    return skipZero_;
}



// tmp<surfaceScalarField> gerhardsRecModel::exportAveragedSurfaceScalarField(word fieldname, scalar threshold, label index)
// {
//     label timeIndex;
//     if (index < 0)
//     {
//         timeIndex = virtualTimeIndex;
//     }
//     else
//     {
//         timeIndex = index;
//     }
//     const label fieldI = getSurfaceScalarFieldIndex(fieldname, timeIndex);
//
//     tmp<surfaceScalarField> tAveragedSurfaceScalarField(surfaceScalarFieldList_[fieldI][timeIndex]);
//
//     label counter = 1;
//     scalar recErr;
//     label delay = 10;
//     label lastMin = -1000;
//
//     for(int runningTimeIndex = 1; runningTimeIndex < numRecFields_-1 ; runningTimeIndex++)
//     {
//         recErr = recurrenceMatrix_[timeIndex][runningTimeIndex];
//         if(recErr > threshold) continue;
//         if(recErr > recurrenceMatrix_[timeIndex][runningTimeIndex-1]) continue;
//         if(recErr > recurrenceMatrix_[timeIndex][runningTimeIndex+1]) continue;
//         if(abs(runningTimeIndex - timeIndex) < delay) continue;
//         if(abs(runningTimeIndex - lastMin) < delay) continue;
//
//         lastMin = runningTimeIndex;
//         counter++;
//         tAveragedSurfaceScalarField += surfaceScalarFieldList_[fieldI][runningTimeIndex];
//     }
//
//     tAveragedSurfaceScalarField /= counter;
//     return tAveragedSurfaceScalarField;
// }

void gerhardsRecModel::exportAveragedVolVectorField(volVectorField& smoothfield, word fieldname, scalar threshold, label index) const
{
    label timeIndex;
    if (index < 0)
    {
        timeIndex = virtualTimeIndex;
    }
    else
    {
        timeIndex = index;
    }
    const label fieldI = getVolVectorFieldIndex(fieldname, timeIndex);

    smoothfield = volVectorFieldList_[fieldI][timeIndex];

    label counter = 1;
    scalar recErr;
    label delay = 1;
    label lastMin = -1000;

    for(int runningTimeIndex = 0; runningTimeIndex < numRecFields_ ; runningTimeIndex++)
    {
        recErr = recurrenceMatrix_[timeIndex][runningTimeIndex];
        if(recErr > threshold) continue;
   //     if(recErr > recurrenceMatrix_[timeIndex][runningTimeIndex-1]) continue;
   //     if(recErr > recurrenceMatrix_[timeIndex][runningTimeIndex+1]) continue;
        if(abs(runningTimeIndex - timeIndex) < delay) continue;
        if(abs(runningTimeIndex - lastMin) < delay) continue;

        lastMin = runningTimeIndex;
        counter++;
        smoothfield += volVectorFieldList_[fieldI][runningTimeIndex];
    }
    Info << "time index = " << index << ", counter = " << counter << endl;
    smoothfield /= counter;
}



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
