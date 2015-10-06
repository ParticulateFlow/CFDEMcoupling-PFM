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

#ifndef recModel_H
#define recModel_H


#include "cfdemCloud.H"
#include "fvCFD.H"
#include "HashTable.H"
#include "labelPair.H"
#include <mpi.h>

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class recModel Declaration
\*---------------------------------------------------------------------------*/

class recModel
{
protected:

    // Protected data
    const dictionary& dict_;
    cfdemCloud& particleCloud_;
    const IOdictionary& controlDict_;
    bool verbose_;
    
    SymmetricSquareMatrix<scalar> recurrenceMatrix;
    Foam::Time recTime;
    instantList timeDirs;
    
    // create a data structure for the time indices
    //	constant will not be contained
    //	runTimeIndex -> continuousIndex
    HashTable<label,word> timeIndexList;
    label contTimeIndex;
    
    // create a data structure for the time values
    //	constant will not be contained
    //	continuousIndex -> time.value()
    HashTable<label,scalar> timeValueList;
    
    label sequenceStart;
    label sequenceEnd;
    const label lowerSeqLim;
    const label upperSeqLim;
   
    label virtualStartIndex;
    
    scalar startTime_;
    scalar endTime_;
    scalar timeStep_;
    scalar recTimeStep_;
    
    label numRecFields;
    label totRecSteps;
    
    label virtualTimeIndex;
    // list of pairs (start value, end value)
    labelPairList virtualTimeIndexList;
    label virtualTimeIndexListPos;
    
    virtual void computeRecMatrix() = 0;
    virtual void computeRecPath() = 0;
    virtual void readFieldSeries() = 0;
    
private:
  
    void readTimeSeries();
    scalar checkTimeStep();

public:

    //- Runtime type information
    TypeName("recModel");

    // Declare runtime constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            recModel,
            dictionary,
            (
                const dictionary& dict,
                cfdemCloud& sm
            ),
            (dict,sm)
        );


    // Constructors

        //- Construct from components
        recModel
        (
            const dictionary& dict,
            cfdemCloud& sm
        );


    // Destructor

        virtual ~recModel();


    // Selector

        static autoPtr<recModel> New
        (
            const dictionary& dict,
            cfdemCloud& sm
        );


    // Member Functions
    
    virtual void initRecFields() = 0; 
    virtual void updateRecFields() = 0;
    virtual void writeRecFields() = 0;
    scalar recTimeStep()
    {
        return recTimeStep_;
    }
    void writeRecMatrix();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
