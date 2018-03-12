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
#include "recModel.H"
#include <unistd.h>


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(recModel, 0);

defineRunTimeSelectionTable(recModel, dictionary);


// * * * * * * * * * * * * * private Member Functions  * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
recModel::recModel
(
    const dictionary& dict,
    recBase& base
)
:
    base_(base),
    recProperties_(dict),
    controlDict_
    (
        IOobject
        (
            "controlDict",
            base.mesh().time().system(),
            base.mesh(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    verbose_(dict.lookupOrDefault<Switch>("verbose", false)),
    volScalarFieldNames_(recProperties_.lookup("volScalarFields")),
    volVectorFieldNames_(recProperties_.lookup("volVectorFields")),
    surfaceScalarFieldNames_(recProperties_.lookup("surfaceScalarFields")),
    startTime_(readScalar(controlDict_.lookup("startTime"))),
    endTime_(readScalar(controlDict_.lookup("endTime"))),
    timeStep_(readScalar(controlDict_.lookup("deltaT"))),
    recStartTime_(-1.0),
    recEndTime_(-1.0),
    totRecSteps_(recProperties_.lookupOrDefault<label>("initialRecSteps",-1)),
    sequenceStart(0),
    sequenceEnd(0),
    virtualStartIndex(0),
    virtualTimeIndex(0),
    virtualTimeIndexNext(1),
    virtualTimeIndexList_(0),
    virtualTimeIndexListPos(0)
{
    recTimeStep_ = -1.0;
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

recModel::~recModel()
{}

// * * * * * * * * * * * * * public Member Functions  * * * * * * * * * * * * //

label recModel::totRecSteps() const
{
    return totRecSteps_;
}
// * * * * * * * * * * * * * private Member Functions  * * * * * * * * * * * * //


void recModel::init()
{
    sequenceStart = virtualTimeIndexList_[0].first();
    sequenceEnd = virtualTimeIndexList_[0].second();
    virtualTimeIndex=sequenceStart;
    virtualTimeIndexNext=virtualTimeIndex+1; 
}

labelPairList& recModel::virtualTimeIndexList()
{
    return virtualTimeIndexList_;
}

void recModel::writeRecPath() const
{
    OFstream listFile("recurrencePath");
    listFile << virtualTimeIndexList_;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
