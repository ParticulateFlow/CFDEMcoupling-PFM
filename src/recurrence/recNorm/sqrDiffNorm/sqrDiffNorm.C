/*---------------------------------------------------------------------------*\
    CFDEMcoupling academic - Open Source CFD-DEM coupling
    
    Contributing authors:
    Thomas Lichtenegger
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
#include "sqrDiffNorm.H"
#include "recModel.H"
#include "addToRunTimeSelectionTable.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(sqrDiffNorm, 0);

addToRunTimeSelectionTable
(
    recNorm,
    sqrDiffNorm,
    dictionary
);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
sqrDiffNorm::sqrDiffNorm
(
    const dictionary& dict,
    recBase& base
)
:
    recNorm(dict, base),
    propsDict_(dict.subDict(typeName + "Props")),
    fieldType_(propsDict_.lookup("fieldType")),
    fieldName_(propsDict_.lookup("fieldName"))
{}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

sqrDiffNorm::~sqrDiffNorm()
{}

// * * * * * * * * * * * * * protected Member Functions  * * * * * * * * * * * * //

void sqrDiffNorm::computeRecMatrix()
{
    Info<< "\nsqrDiffNorm: computing recurrence matrix\n" << endl;
    
    const HashTable<label,word>& timeIndexList( base_.recM().timeIndexList() );
    SymmetricSquareMatrix<scalar>& recurrenceMatrix( base_.recM().recurrenceMatrix() );
    
    scalar normIJ(0.0);
    scalar maxNormIJ(0.0);

    // perform un-normalized calculation for lower half of the recurrence matrix
    forAll(timeIndexList, ti)
    {
    	forAll(timeIndexList, tj)
    	{
	        if(verbose_)
		    Info<<"\n Doing calculation for element " << ti << " " << tj << "\n" << endl;
		
    		// skip main diagonal and upper half
    		if (ti >= tj)
    		{
    			recurrenceMatrix[ti][tj] = 0;
    			continue;
    		}
    		
    		// compute elements
    		if (fieldType_ == "volScalarField")
		    normIJ = normVSF(ti,tj);
		else if (fieldType_ == "volVectorField")
		    normIJ = normVVF(ti,tj);
		else if (fieldType_ == "surfaceScalarField")
		    normIJ = normSSF(ti,tj);
		else
		    FatalError<<"sqrDiffNorm: Unknown field type " << fieldType_ << abort(FatalError);
		
    		recurrenceMatrix[ti][tj] = normIJ;
		
		if (normIJ > maxNormIJ)
		    maxNormIJ = normIJ;
    	}
    }
    
    // normalize matrix and copy lower to upper half
    forAll(timeIndexList, ti)
    {
    	forAll(timeIndexList, tj)
    	{
	    if (ti >= tj)
	        continue;
	    
	    recurrenceMatrix[ti][tj] /= maxNormIJ;
	    recurrenceMatrix[tj][ti] = recurrenceMatrix[ti][tj];
	}
    }
    
    Info<< "\nComputing recurrence matrix done\n" << endl;
}

scalar sqrDiffNorm::normVSF(label ti, label tj)
{
    const volScalarField& t1( base_.recM().exportVolScalarField(fieldName_,ti) );
    const volScalarField& t2( base_.recM().exportVolScalarField(fieldName_,tj) );
    dimensionedScalar tNorm( fvc::domainIntegrate( sqr( t1 - t2 ) ) );
    return tNorm.value();
}

scalar sqrDiffNorm::normVVF(label ti, label tj)
{
    const volVectorField& t1( base_.recM().exportVolVectorField(fieldName_,ti) );
    const volVectorField& t2( base_.recM().exportVolVectorField(fieldName_,tj) );
    dimensionedScalar tNorm( fvc::domainIntegrate( magSqr( t1 - t2 ) ) );
    return tNorm.value();
}

scalar sqrDiffNorm::normSSF(label ti, label tj)
{
    const surfaceScalarField& t1( base_.recM().exportSurfaceScalarField(fieldName_,ti) );
    const surfaceScalarField& t2( base_.recM().exportSurfaceScalarField(fieldName_,tj) );
    volVectorField t12 (fvc::reconstruct( t1-t2 ) );
    dimensionedScalar tNorm( fvc::domainIntegrate( magSqr( t12 ) ) );
    return tNorm.value();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
