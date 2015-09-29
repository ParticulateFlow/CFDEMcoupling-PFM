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
#include "IOModel.H"
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
    cfdemCloud& sm
)
:
    recModel(dict,sm)
{
    readFieldSeries();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

standardRecModel::~standardRecModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void standardRecModel::readFieldSeries()
{
     Info << "Reading fields" << endl;
    
    for (instantList::iterator it=timeDirs.begin(); it != timeDirs.end(); ++it)
    {
        // set time
        recTime.setTime(*it, it->value());
        
        // skip constant
        if (recTime.timeName() == "constant")
        {
        	continue;
        }
        
        if (verbose)
    	{
        	Info << "Reading at t = " << recTime.timeName() << endl;
        }
        
        alpha1pl.set
        (
            timeIndexList(recTime.timeName()),
            new volScalarField
            (
                IOobject
                (
                    "voidfraction",
                    recTime.timePath(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                ),
                mesh
            )
        );
        
        IOobject headerAlpha2
        (
            "particlefraction",
            recTime.timePath(),
            mesh,
            IOobject::NO_READ
        );
        
        if (headerAlpha2.headerOk())
        {
	        alpha2pl.set
    	    (
    	        timeIndexList(recTime.timeName()),
    	        new volScalarField
    	        (
    	            IOobject
    	            (
    	                "particlefraction",
                    	recTime.timePath(),
    	                mesh,
    	                IOobject::MUST_READ,
    	                IOobject::NO_WRITE
    	            ),
    	            mesh
    	        )
    	    );
    	}
    	else
    	{
    		alpha2pl.set
    	    (
    	        timeIndexList(recTime.timeName()),
    	        new volScalarField
    	        (
    	            IOobject
    	            (
    	                "particlefraction",
                    	recTime.timePath(),
    	                mesh,
    	                IOobject::NO_READ,
    	                IOobject::NO_WRITE
    	            ),
    	            scalar(1.0) - alpha1pl[timeIndexList(recTime.timeName())]
    	        )
    	    );
    	}
        
        U1pl.set
        (
        	timeIndexList(recTime.timeName()),
            new volVectorField
            (
                IOobject
                (
                    "U",
                    recTime.timePath(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                ),
                mesh
            )
        );
        U2pl.set
        (
        	timeIndexList(recTime.timeName()),
            new volVectorField
            (
                IOobject
                (
                    "Us",
                    recTime.timePath(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                ),
                mesh
            )
        );
        
    } 
}

void standardRecModel::computeRecMatrix()
{
    Info<< "\nComputing recurrence matrix\n" << endl;
    recurrenceMatrix(timeIndexList.size(), timeIndexList.size(), scalar(0));
    scalar maxElemVal(0.0);
    
    // compute recurrence matrix elements
    forAll(timeIndexList, ti)
    {
    	forAll(timeIndexList, tj)
    	{
    		// main diagonal
    		if (ti == tj)
    		{
    			recurrenceMatrix[ti][tj] = 1;
    			continue;
    		}
    		
    		// skip one half of the matrix
    		if (ti > tj)
    		{
    			continue;
    		}
    		
    		// compute elements
    		recurrenceMatrix[ti][tj]
    			= sumSqr(alpha1pl[ti].internalField() - alpha1pl[tj].internalField())
    			+ sum(magSqr(U1pl[ti].internalField() - U1pl[tj].internalField()));
    		
    		recurrenceMatrix[tj][ti] = recurrenceMatrix[ti][tj];
    		
    		if (maxElemVal < recurrenceMatrix[ti][tj])
    		{
    			maxElemVal = recurrenceMatrix[ti][tj];
    		}
    	}
    }
    
    // normalize matrix elements
    forAll(timeIndexList, ti)
    {
    	forAll(timeIndexList, tj)
    	{
    		recurrenceMatrix[ti][tj] /= maxElemVal;
    	}
    } 
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
