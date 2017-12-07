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
#include <mpi.h>

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
    volScalarFieldList_(volScalarFieldNames_.size()),
    volVectorFieldList_(volVectorFieldNames_.size()),
    surfaceScalarFieldList_(surfaceScalarFieldNames_.size())
{
    for(int i=0; i<volScalarFieldNames_.size(); i++)
        volScalarFieldList_[i].setSize(numRecFields_);
    
    for(int i=0; i<volVectorFieldNames_.size(); i++)
        volVectorFieldList_[i].setSize(numRecFields_);
    
    for(int i=0; i<surfaceScalarFieldNames_.size(); i++)
        surfaceScalarFieldList_[i].setSize(numRecFields_);
       

    
 //   setRecFields();
    

}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

standardRecModel::~standardRecModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void standardRecModel::updateRecFields()
{
  virtualTimeIndex=virtualTimeIndexNext;
  virtualTimeIndexNext++;
  if (virtualTimeIndexNext>sequenceEnd)
  {
    virtualTimeIndexListPos++;
    sequenceStart=virtualTimeIndexList_[virtualTimeIndexListPos].first();
    sequenceEnd=virtualTimeIndexList_[virtualTimeIndexListPos].second();
    virtualTimeIndexNext=sequenceStart;
  }
  
  if (verbose_)
      Info << "\nUpdating virtual time index to " << virtualTimeIndex << ".\n" << endl;  
}




const volScalarField& standardRecModel::exportVolScalarField(word fieldname, label index) const
{
    const label fieldI = getVolScalarFieldIndex(fieldname, index);
    
    return volScalarFieldList_[fieldI][index];
}

const volVectorField& standardRecModel::exportVolVectorField(word fieldname, label index) const
{
    const label fieldI = getVolVectorFieldIndex(fieldname, index);
    
    return volVectorFieldList_[fieldI][index];
}

const surfaceScalarField& standardRecModel::exportSurfaceScalarField(word fieldname, label index) const
{
    const label fieldI = getSurfaceScalarFieldIndex(fieldname, index);
    
    return surfaceScalarFieldList_[fieldI][index];
}

void standardRecModel::exportVolScalarField(word fieldname, volScalarField& field) const
{
    field = exportVolScalarField(fieldname, virtualTimeIndex);
}

void standardRecModel::exportVolVectorField(word fieldname, volVectorField& field) const
{
    field = exportVolVectorField(fieldname, virtualTimeIndex);
}

void standardRecModel::exportSurfaceScalarField(word fieldname, surfaceScalarField& field) const
{
    field = exportSurfaceScalarField(fieldname, virtualTimeIndex);
}

// tmp<surfaceScalarField> standardRecModel::exportAveragedSurfaceScalarField(word fieldname, scalar threshold) const
// {
//     const label fieldI = getSurfaceScalarFieldIndex(fieldname, virtualTimeIndex);
//     
//     tmp<surfaceScalarField> tAveragedSurfaceScalarField(surfaceScalarFieldList_[fieldI][virtualTimeIndex]);
//     
//     label counter = 1;
//     scalar recErr;
//     label delay = 10;
//     label lastMin = -1000;
//     
//     for(int runningTimeIndex = 1; runningTimeIndex < numRecFields_-1 ; runningTimeIndex++)
//     {
//         recErr = recurrenceMatrix_[virtualTimeIndex][runningTimeIndex];
//         if(recErr > threshold) continue;
//         if(recErr > recurrenceMatrix_[virtualTimeIndex][runningTimeIndex-1]) continue;
//         if(recErr > recurrenceMatrix_[virtualTimeIndex][runningTimeIndex+1]) continue;
//         if(abs(runningTimeIndex - virtualTimeIndex) < delay) continue;
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


void standardRecModel::readFieldSeries()
{
    Info << "Reading fields\n" << endl;
    
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
        
        // skip constant
        if (recTime.timeName() == "constant")
        {
        	continue;
        }
        
        if (verbose_)
    	{
        	Info << "Reading at t = " << recTime.timeName() << endl;
        }
        
        for(int i=0; i<volScalarFieldNames_.size(); i++)
            volScalarFieldList_[i].set
            (
	        timeIndexList_(recTime.timeName()),
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
	    
	for(int i=0; i<volVectorFieldNames_.size(); i++)
            volVectorFieldList_[i].set
            (
	        timeIndexList_(recTime.timeName()),
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
	    
	for(int i=0; i<surfaceScalarFieldNames_.size(); i++)
            surfaceScalarFieldList_[i].set
            (
	        timeIndexList_(recTime.timeName()),
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
    Info << "Reading fields done" << endl;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
