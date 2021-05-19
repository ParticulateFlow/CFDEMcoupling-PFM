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
#include "standardRecModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(standardRecModel, 0);

addToRunTimeSelectionTable
(
    recModel,
    standardRecModel,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
standardRecModel::standardRecModel
(
    const dictionary& dict,
    recBase& base
)
:
    recModel(dict,base),
    propsDict_(dict.subDict(typeName + "Props")),
    dataBaseNames_(propsDict_.lookupOrDefault<wordList>("dataBases", wordList(1,"dataBase"))),
    numDataBases_(dataBaseNames_.size()),
    timeDirs_(),
    skipZero_(propsDict_.lookupOrDefault<Switch>("skipZero", Switch(false))),
    numRecFields_(),
    cumulativeNumRecFields_(),
    totNumRecFields_(0),
    storeAveragedFields_(propsDict_.lookupOrDefault<bool>("storeAveragedFields",false)),
    checkTimeStep_(propsDict_.lookupOrDefault<bool>("checkTimeStep",true)),
    recurrenceMatrix_(1,scalar(-1.0)),
    timeIndexList_(),
    timeValueList_(),
    contTimeIndex(0),
    recTimeDilationFactor_(propsDict_.lookupOrDefault<scalar>("timeDilationFactor", 1.0)),
    volScalarFieldList_(volScalarFieldNames_.size()),
    volVectorFieldList_(volVectorFieldNames_.size()),
    surfaceScalarFieldList_(surfaceScalarFieldNames_.size()),
    aveVolScalarFieldList_(volScalarFieldNames_.size()),
    aveVolVectorFieldList_(volVectorFieldNames_.size()),
    aveSurfaceScalarFieldList_(surfaceScalarFieldNames_.size())
{
    for(label i=0;i<numDataBases_;i++)
    {
        Foam::Time recTime(fileName(dataBaseNames_[i]), "", "../system", "../constant", false);

        instantList timeDirs(recTime.times());
        if (timeDirs.size() == 0)
        {
            FatalError <<"database " << dataBaseNames_[i] << " does not exist or is empty\n" << abort(FatalError);
        }
        else if (timeDirs.size() < 3)
        {
            FatalError <<"database " << dataBaseNames_[i] << " contains only one or two entries. "<<
            "To describe a steady state, create copies of the steady fields so that the database contains "<<
            "entries with three identical sets of fields.\n" << abort(FatalError);
        } 
        timeDirs_.append(timeDirs);

        numRecFields_.append(label(timeDirs_[i].size()));
        if(skipZero_) numRecFields_[i]--;

        label sum=0;
        for(label j=0;j<=i;j++)
        {
            sum += numRecFields_[j];
        }
        cumulativeNumRecFields_.append(sum);

        totNumRecFields_ += numRecFields_[i];

        if (verbose_)
        {
            // be informative on properties of the "recTime" Time-object
            Info << "information on database " << i << endl;
            Info << "recTime.rootPath() " << recTime.rootPath() << endl;
            Info << "recTime.caseName() " << recTime.caseName() << endl;
            Info << "recTime.path() " << recTime.path() << endl;
            Info << "recTime.timePath() " << recTime.timePath() << endl;
            Info << "recTime.timeName() " << recTime.timeName() << endl;
            Info << "timeDirs " << timeDirs_[i] << endl;
            Info << "ignoring 0 directory: " << skipZero_ << endl;
        }
    }
    Info << "total number of recurrence fields = " << totNumRecFields_ << endl;
    recurrenceMatrix_.setSize(totNumRecFields_,totNumRecFields_);

    lowerSeqLim_ = max(1, label(totNumRecFields_/20));
    upperSeqLim_ = max(1, label(totNumRecFields_/5));

    readTimeSeries();

    // check if time steps in databases are consistent
    // if no initial number of time steps has been specified, create path for full runtime immediately
    if(checkTimeStep_)
    {
        recTimeStep_ = checkTimeStep();
    }

    if(totRecSteps_ < 0)
    {
        totRecSteps_ = 1 + static_cast<label>( (endTime_-startTime_) / recTimeStep_ );
    }

    for(int i=0; i<volScalarFieldNames_.size(); i++)
    {
        volScalarFieldList_[i].setSize(totNumRecFields_);
        if(storeAveragedFields_) aveVolScalarFieldList_[i].setSize(numDataBases_);
    }

    for(int i=0; i<volVectorFieldNames_.size(); i++)
    {
        volVectorFieldList_[i].setSize(totNumRecFields_);
        if(storeAveragedFields_) aveVolVectorFieldList_[i].setSize(numDataBases_);
    }

    for(int i=0; i<surfaceScalarFieldNames_.size(); i++)
    {
        surfaceScalarFieldList_[i].setSize(totNumRecFields_);
        if(storeAveragedFields_) aveSurfaceScalarFieldList_[i].setSize(numDataBases_);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

standardRecModel::~standardRecModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
scalar standardRecModel::checkTimeStep()
{
    // check time step of provided data
    scalar dtCur(1.e10);
    scalar dtOld(1.e10);
    bool first = true;
    scalar frac = 0.0;
    scalar tolerance = 1e-5;
    for(int i=0;i<numDataBases_;i++)
    {
        for(int j=0;j<numRecFields_[i]-1;j++)
        {
            dtOld = dtCur;
            dtCur = timeDirs_[i][j+1].value() - timeDirs_[i][j].value();
            frac = 1 - dtOld/dtCur;
            if(!first && fabs(frac) > tolerance)
            {
                FatalError <<"detected different time steps in database(s)\n" << abort(FatalError);
            }
        }
    }

    if (verbose_)
    {
        Info << "Setting deltaRecT to " << dtCur << endl;
        Info << "Actual runTime.deltaT = " << timeStep_ << endl;
    }

    return dtCur*recTimeDilationFactor_;
}


void standardRecModel::init()
{
    recModel::init();

    const objectRegistry& objReg = base_.mesh().thisDb();

    for(int j=0; j<volScalarFieldNames_.size(); j++)
    {
        objReg.checkIn(volScalarFieldList_[j][virtualTimeIndex]);
    }

    for(int j=0; j<volVectorFieldNames_.size(); j++)
    {
        objReg.checkIn(volVectorFieldList_[j][virtualTimeIndex]);
    }

    for(int j=0; j<surfaceScalarFieldNames_.size(); j++)
    {
        objReg.checkIn(surfaceScalarFieldList_[j][virtualTimeIndex]);
    }


    for(int i = 0; i < numDataBases_; i++)
    {
        if (virtualTimeIndex < cumulativeNumRecFields_[i])
        {
            currDataBase_ = i;
            break;
        }
    }

    currDataBaseNext_ = currDataBase_;
}


void standardRecModel::averageFieldSeries()
{
    for(int i=0;i<numDataBases_;i++)
    {
        Foam::Time recTime(fileName(dataBaseNames_[i]), "", "../system", "../constant", false);

        // perform averaging over all volScalarFields
        for(int j=0; j<volScalarFieldNames_.size(); j++)
        {
            dimensionedScalar dimZero("zero",volScalarFieldList_[j][0].dimensions(),0.0);
            aveVolScalarFieldList_[j].set
            (
                i,
                new volScalarField
                (
                    IOobject
                    (
                        "ave"+volScalarFieldNames_[j],
                        recTime.timePath(),
                        base_.mesh(),
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ),
                    base_.mesh(),
                    dimZero
                )
            );

            label k1 = 0;
            if(i>0) k1 = cumulativeNumRecFields_[i-1];
            label k2 = cumulativeNumRecFields_[i];
            for(label k = k1; k < k2; k++)
            {
                aveVolScalarFieldList_[j][i] += volScalarFieldList_[j][k];
            }
            aveVolScalarFieldList_[j][i] /= k2 - k1;
        }

        // perform averaging over all volVectorFields
        for(int j=0; j<volVectorFieldNames_.size(); j++)
        {
            dimensionedVector dimZero("zero",volVectorFieldList_[j][0].dimensions(),vector::zero);
            aveVolVectorFieldList_[j].set
            (
                i,
                new volVectorField
                (
                    IOobject
                    (
                        "ave"+volVectorFieldNames_[j],
                        recTime.timePath(),
                        base_.mesh(),
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ),
                    base_.mesh(),
                    dimZero
                )
            );

            label k1 = 0;
            if(i>0) k1 = cumulativeNumRecFields_[i-1];
            label k2 = cumulativeNumRecFields_[i];
            for(label k = k1; k < k2; k++)
            {
                aveVolVectorFieldList_[j][i] += volVectorFieldList_[j][k];
            }
            aveVolVectorFieldList_[j][i] /= k2 - k1;
        }

        // perform averaging over all surfaceScalarFields
        for(int j=0; j<surfaceScalarFieldNames_.size(); j++)
        {
            dimensionedScalar dimZero("zero",surfaceScalarFieldList_[j][0].dimensions(),0.0);
            aveSurfaceScalarFieldList_[j].set
            (
                i,
                new surfaceScalarField
                (
                    IOobject
                    (
                        "ave"+surfaceScalarFieldNames_[j],
                        recTime.timePath(),
                        base_.mesh(),
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ),
                    base_.mesh(),
                    dimZero
                )
            );

            label k1 = 0;
            if(i>0) k1 = cumulativeNumRecFields_[i-1];
            label k2 = cumulativeNumRecFields_[i];
            for(label k = k1; k < k2; k++)
            {
                aveSurfaceScalarFieldList_[j][i] += surfaceScalarFieldList_[j][k];
            }
            aveSurfaceScalarFieldList_[j][i] /= k2 - k1;
        }
    }
}


void standardRecModel::readFieldSeries()
{
    label fieldcounter = 0;
    for(int i=0;i<numDataBases_;i++)
    {
        Info << "Reading fields of database " << dataBaseNames_[i] << "\n" << endl;

        Foam::Time recTime(fileName(dataBaseNames_[i]), "", "../system", "../constant", false);

        label size = timeDirs_[i].size();
        label counter = 0;
        label percentage = 0;

        for (instantList::iterator it=timeDirs_[i].begin(); it != timeDirs_[i].end(); ++it)
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
                Info << "Reading at t = " << recTime.timeName() << endl;
            }

            for(int j=0; j<volScalarFieldNames_.size(); j++)
            {
                volScalarFieldList_[j].set
                (
                    fieldcounter,
                    new volScalarField
                    (
                        IOobject
                        (
                            volScalarFieldNames_[j],
                            recTime.timePath(),
                            base_.mesh(),
                            IOobject::MUST_READ,
                            IOobject::NO_WRITE
                        ),
                        base_.mesh()
                    )
                );
            }

            for(int j=0; j<volVectorFieldNames_.size(); j++)
            {
                volVectorFieldList_[j].set
                (
                    fieldcounter,
                    new volVectorField
                    (
                        IOobject
                        (
                            volVectorFieldNames_[j],
                            recTime.timePath(),
                            base_.mesh(),
                            IOobject::MUST_READ,
                            IOobject::NO_WRITE
                        ),
                        base_.mesh()
                    )
                );
            }

            for(int j=0; j<surfaceScalarFieldNames_.size(); j++)
            {
                surfaceScalarFieldList_[j].set
                (
                    fieldcounter,
                    new surfaceScalarField
                    (
                        IOobject
                        (
                            surfaceScalarFieldNames_[j],
                            recTime.timePath(),
                            base_.mesh(),
                            IOobject::MUST_READ,
                            IOobject::NO_WRITE
                        ),
                        base_.mesh()
                    )
                );
            }
            fieldcounter++;
        }
        Info << "Reading fields of database " << dataBaseNames_[i] <<" done" << endl;
    }

    if(storeAveragedFields_)
    {
        Info << "Calculating field averages." << endl;
        averageFieldSeries();
        Info << "Calculating field averages done." << endl;
    }
}


void standardRecModel::readTimeSeries()
{
    for(int i=0;i<numDataBases_;i++)
    {
        Foam::Time recTime(fileName(dataBaseNames_[i]), "", "../system", "../constant", false);
        bool firsttime = true;
        // fill the data structure for the time indices
        for (instantList::iterator it=timeDirs_[i].begin(); it != timeDirs_[i].end(); ++it)
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

        recEndTime_ = recStartTime_ + (recEndTime_ - recStartTime_) * recTimeDilationFactor_;

        if (verbose_)
        {
            Info << endl;
            Info << "Found " << label(timeDirs_[i].size()) << " time folders" << endl;
            Info << "Found " << label(timeIndexList_.size()) << " time steps" << endl;
            Info << "database start time = " << recStartTime_ << endl;
            Info << "database end time = " << recEndTime_ << endl;
        }
    }
}


void standardRecModel::exportVolScalarFieldAve(word fieldname, volScalarField& field, label db)
{
    if(!storeAveragedFields_)
    {
        FatalError <<"no averaged fields available, need to activate \"storeAveragedFields\"\n" << abort(FatalError);
    }
    if (db >= numDataBases_)
    {
        FatalError <<"can't find database with number " << db << abort(FatalError);
    }
    const label fieldI = getVolScalarFieldIndex(fieldname, 0);
    field = aveVolScalarFieldList_[fieldI][db];
}


void standardRecModel::exportVolVectorFieldAve(word fieldname, volVectorField& field, label db)
{
    if(!storeAveragedFields_)
    {
        FatalError <<"no averaged fields available, need to activate \"storeAveragedFields\"\n" << abort(FatalError);
    }
    if (db >= numDataBases_)
    {
        FatalError <<"can't find database with number " << db << abort(FatalError);
    }
    const label fieldI = getVolVectorFieldIndex(fieldname, 0);
    field = aveVolVectorFieldList_[fieldI][db];
}


void standardRecModel::exportSurfaceScalarFieldAve(word fieldname, surfaceScalarField& field, label db)
{
    if(!storeAveragedFields_)
    {
        FatalError <<"no averaged fields available, need to activate \"storeAveragedFields\"\n" << abort(FatalError);
    }
    if (db >= numDataBases_)
    {
        FatalError <<"can't find database with number " << db << abort(FatalError);
    }
    const label fieldI = getSurfaceScalarFieldIndex(fieldname, 0);
    field = aveSurfaceScalarFieldList_[fieldI][db];
}

tmp<volScalarField> standardRecModel::exportVolScalarFieldAve(word fieldname, label db)
{
    if(!storeAveragedFields_)
    {
        FatalError <<"no averaged fields available, need to activate \"storeAveragedFields\"\n" << abort(FatalError);
    }
    if (db >= numDataBases_)
    {
        FatalError <<"can't find database with number " << db << abort(FatalError);
    }
    const label fieldI = getVolScalarFieldIndex(fieldname, 0);
    return aveVolScalarFieldList_[fieldI][db];
}

tmp<volVectorField> standardRecModel::exportVolVectorFieldAve(word fieldname, label db)
{
    if(!storeAveragedFields_)
    {
        FatalError <<"no averaged fields available, need to activate \"storeAveragedFields\"\n" << abort(FatalError);
    }
    if (db >= numDataBases_)
    {
        FatalError <<"can't find database with number " << db << abort(FatalError);
    }
    const label fieldI = getVolVectorFieldIndex(fieldname, 0);
    return aveVolVectorFieldList_[fieldI][db];
}

tmp<surfaceScalarField> standardRecModel::exportSurfaceScalarFieldAve(word fieldname, label db)
{
    if(!storeAveragedFields_)
    {
        FatalError <<"no averaged fields available, need to activate \"storeAveragedFields\"\n" << abort(FatalError);
    }
    if (db >= numDataBases_)
    {
        FatalError <<"can't find database with number " << db << abort(FatalError);
    }
    const label fieldI = getSurfaceScalarFieldIndex(fieldname, 0);
    return aveSurfaceScalarFieldList_[fieldI][db];
}

void standardRecModel::exportVolScalarField(word fieldname, volScalarField& field)
{
    field = exportVolScalarField(fieldname, virtualTimeIndex);
}


void standardRecModel::exportVolVectorField(word fieldname, volVectorField& field)
{
    field = exportVolVectorField(fieldname, virtualTimeIndex);
}


void standardRecModel::exportSurfaceScalarField(word fieldname, surfaceScalarField& field)
{
    field = exportSurfaceScalarField(fieldname, virtualTimeIndex);
}


const volScalarField& standardRecModel::exportVolScalarField(word fieldname, label index)
{
    const label fieldI = getVolScalarFieldIndex(fieldname, index);

    return volScalarFieldList_[fieldI][index];
}

const volVectorField& standardRecModel::exportVolVectorField(word fieldname, label index)
{
    const label fieldI = getVolVectorFieldIndex(fieldname, index);

    return volVectorFieldList_[fieldI][index];
}

const surfaceScalarField& standardRecModel::exportSurfaceScalarField(word fieldname, label index)
{
    const label fieldI = getSurfaceScalarFieldIndex(fieldname, index);

    return surfaceScalarFieldList_[fieldI][index];
}

PtrList<volScalarField>& standardRecModel::exportVolScalarFieldList(word fieldname)
{
    const label fieldI = getVolScalarFieldIndex(fieldname);

    return volScalarFieldList_[fieldI];
}

PtrList<volVectorField>& standardRecModel::exportVolVectorFieldList(word fieldname)
{
    const label fieldI = getVolVectorFieldIndex(fieldname);

    return volVectorFieldList_[fieldI];
}

PtrList<surfaceScalarField>& standardRecModel::exportSurfaceScalarFieldList(word fieldname)
{
    const label fieldI = getSurfaceScalarFieldIndex(fieldname);

    return surfaceScalarFieldList_[fieldI];
}


SymmetricSquareMatrix<scalar>& standardRecModel::recurrenceMatrix()
{
    return recurrenceMatrix_;
}


const HashTable<label,word>& standardRecModel::timeIndexList() const
{
    return timeIndexList_;
}


label standardRecModel::lowerSeqLim() const
{
    return lowerSeqLim_;
}

label standardRecModel::upperSeqLim() const
{
    return upperSeqLim_;
}

label standardRecModel::numIntervals() const
{
    return numDataBases_;
}

label standardRecModel::numRecFields() const
{
    return totNumRecFields_;
}

label standardRecModel::numRecFields(label i) const
{
    return numRecFields_[i];
}


label standardRecModel::numDataBaseFields() const
{
    return totNumRecFields_;
}


void standardRecModel::updateRecFields()
{
    // make fields of upcoming virtualTimeIndex available in object registry
    const objectRegistry& objReg = base_.mesh().thisDb();

    for(int j=0; j<volScalarFieldNames_.size(); j++)
    {
        objReg.checkOut(volScalarFieldList_[j][virtualTimeIndex]);
        objReg.checkIn(volScalarFieldList_[j][virtualTimeIndexNext]);
    }

    for(int j=0; j<volVectorFieldNames_.size(); j++)
    {
        objReg.checkOut(volVectorFieldList_[j][virtualTimeIndex]);
        objReg.checkIn(volVectorFieldList_[j][virtualTimeIndexNext]);
    }

    for(int j=0; j<surfaceScalarFieldNames_.size(); j++)
    {
        objReg.checkOut(surfaceScalarFieldList_[j][virtualTimeIndex]);
        objReg.checkIn(surfaceScalarFieldList_[j][virtualTimeIndexNext]);
    }

    virtualTimeIndex = virtualTimeIndexNext;
    virtualTimeIndexNext++;

    currDataBase_ = currDataBaseNext_;

    if (virtualTimeIndexNext > sequenceEnd)
    {
        virtualTimeIndexListPos_++; // this is problematic with noPath

        sequenceStart = virtualTimeIndexList_[virtualTimeIndexListPos_].first();
        sequenceEnd = virtualTimeIndexList_[virtualTimeIndexListPos_].second();

        virtualTimeIndexNext = sequenceStart;

        if (verbose_)
        {
            Info << " new sequence (start/end) : " << sequenceStart << " / " << sequenceEnd << endl;
        }

        // check in which DB the new interval lies
        for(int i = 0; i < numDataBases_; i++)
        {
            if (virtualTimeIndexNext < cumulativeNumRecFields_[i])
            {
                currDataBaseNext_ = i;
                break;
            }
        }
        lastJumpTime_.set("lastJumpTime",base_.mesh().time().timeOutputValue());
    }

    if (verbose_)
    {
        Info << "\nUpdating virtual time index to " << virtualTimeIndex << ".\n" << endl;
    }
}

label standardRecModel::currentTimeIndex() const
{
    return virtualTimeIndex;
}

void standardRecModel::writeRecMatrix() const
{
    if (writeRecMat_)
    {
        OFstream matrixFile(recMatName_);
        matrixFile << recurrenceMatrix_;
    }
}



// tmp<surfaceScalarField> standardRecModel::exportAveragedSurfaceScalarField(word fieldname, scalar threshold, label index)
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

void standardRecModel::exportAveragedVolVectorField(volVectorField& smoothfield, word fieldname, scalar threshold, label index) const
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

    for(int runningTimeIndex = 0; runningTimeIndex < totNumRecFields_ ; runningTimeIndex++)
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
