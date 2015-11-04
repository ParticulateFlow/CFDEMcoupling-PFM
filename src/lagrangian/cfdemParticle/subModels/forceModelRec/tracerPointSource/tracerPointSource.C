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

#include "tracerPointSource.H"
#include "addToRunTimeSelectionTable.H"
#include "dataExchangeModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(tracerPointSource, 0);

addToRunTimeSelectionTable
(
    forceModelRec,
    tracerPointSource,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
tracerPointSource::tracerPointSource
(
    const dictionary& dict,
    cfdemCloudRec& sm
)
:
    forceModelRec(dict,sm),
    propsDict_(dict.subDict(typeName + "Props")),
    tracerConcentrationField_
    (   IOobject
        (
            "tracerConcentration",
            sm.mesh().time().timeName(),
            sm.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        sm.mesh(),
        dimensionedScalar("zero", dimensionSet(0,0,0,0,0), 0.0)
    ),
    tracerConcentrationPart_(NULL),
    sourcePosition_(propsDict_.lookup("sourcePosition")),
    sourceStrength_(readScalar(propsDict_.lookup("sourceStrength"))),
    beginTime_(propsDict_.lookupOrDefault<scalar>("beginTime",0.0)),
    endTime_(propsDict_.lookupOrDefault<scalar>("endTime",GREAT)),
    numReceiverPart_(propsDict_.lookupOrDefault<label>("numReceiverPart",1)),
    itID_(NULL),
    itDist_(NULL)
{
    if (numReceiverPart_<1)
    {
         FatalError
            << "tracerPointSource::tracerPointSource() : "
            << endl
            << "    trying to set invalid number of receiver particles. "<< endl << endl
            << abort(FatalError);  
    }
    
    // distribute source strength by number of receiving particles
    sourceStrength_/=numReceiverPart_;
    
    receiverIDs_.reserve(numReceiverPart_);
    receiverDist_.reserve(numReceiverPart_);
    allocateMyArrays();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

tracerPointSource::~tracerPointSource()
{
    delete tracerConcentrationPart_;
}

// * * * * * * * * * * * * * * * private Member Functions  * * * * * * * * * * * * * //
void tracerPointSource::allocateMyArrays() const
{
    // get memory for 2d arrays
    double initVal=0.0;
    particleCloud_.dataExchangeM().allocateArray(tracerConcentrationPart_,initVal,1);
}
// * * * * * * * * * * * * * * * public Member Functions  * * * * * * * * * * * * * //

void tracerPointSource::setForce() const
{
    scalar currTime=particleCloud_.mesh().time().value();
    if (currTime < beginTime_ || currTime > endTime_)
        return;
    
    // realloc the arrays
    allocateMyArrays();
    
    // reset Scalar field
    tracerConcentrationField_.internalField() = 0.0;

    // get DEM data
    particleCloud_.dataExchangeM().getData("tracerConcentration","scalar-atom",tracerConcentrationPart_);
    
    label cellI=0;
    scalar distance=0.0;
    vector position;
    receiverIDs_.assign(numReceiverPart_, -1);
    receiverDist_.assign(numReceiverPart_, GREAT);
    
    for(int index = 0;index < particleCloud_.numberOfParticles(); ++index)
    {
        cellI = particleCloud_.cellIDs()[index][0];
        if(cellI >= 0)
        {
            position = particleCloud_.position(index);
	    distance = mag(position-sourcePosition_);
            if (distance < receiverDist_[numReceiverPart_-1])
	    {
	        int pos;
		for(pos=numReceiverPart_-2; pos>=0; pos--)
		    if( distance > receiverDist_[pos] )
		        break;
		receiverIDs_.pop_back();
		receiverDist_.pop_back();
                pos++;
                itID_ = receiverIDs_.begin()+pos;
                itDist_ = receiverDist_.begin()+pos;
		receiverIDs_.insert(itID_,index);
		receiverDist_.insert(itDist_,distance);
	    }           
        }
    }
    
    // check if enough particles were found
    if(receiverIDs_[numReceiverPart_-1] < 0)
    {
         FatalError
            << "tracerPointSource::setForce() : "
            << endl
            << "    not enough receiver particles found. "<< endl << endl
            << abort(FatalError);  
    }
    
    
    for(int i=0;i<numReceiverPart_;i++)
    {
        int index=receiverIDs_[i];
        tracerConcentrationPart_[index][0] += sourceStrength_;
    }

    particleCloud_.averagingM().setScalarSum
    (
        tracerConcentrationField_,
        tracerConcentrationPart_,
        particleCloud_.particleWeights(),
        NULL
    );
    // give DEM data
    particleCloud_.dataExchangeM().giveData("tracerConcentration","scalar-atom", tracerConcentrationPart_);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
