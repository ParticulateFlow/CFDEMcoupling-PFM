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
#include "standardRecStatAnalysis.H"
#include "recModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(standardRecStatAnalysis, 0);

addToRunTimeSelectionTable
(
    recStatAnalysis,
    standardRecStatAnalysis,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
standardRecStatAnalysis::standardRecStatAnalysis
(
    const dictionary& dict,
    recBase& base
)
:
    recStatAnalysis(dict,base),
    propsDict_(dict.subDict(typeName + "Props")),
    aveSimSpecialRefLength_(propsDict_.lookupOrDefault<label>("aveSimSpecialRefLength",0)),
    ignoreStepsInit_(propsDict_.lookupOrDefault<label>("ignoreStepsInit",0)),
    ignoreStepsDiag_(propsDict_.lookupOrDefault<label>("ignoreStepsDiag",0)),
    aveSimilarity_(),
    mostDistinct_()
{

}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

standardRecStatAnalysis::~standardRecStatAnalysis()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
void standardRecStatAnalysis::init()
{
    label size = base_.recM().numRecFields() - ignoreStepsInit_;
    aveSimilarity_.resize(size,0.0);
    correlationSum_.resize(size,0.0);
    mostDistinct_.resize(size,0.0);
}

void standardRecStatAnalysis::statistics()
{
    aveSimGeneral();
    if(aveSimSpecialRefLength_ > ignoreStepsInit_)
        aveSimSpecial();
    mostDistinct();
}

void standardRecStatAnalysis::aveSimGeneral()
{
    Info << "Calculating general average similarity.\n" << endl;
    label Nmax = base_.recM().numRecFields();
    scalar r = 0.0;
    scalar Rij = 0.0;
    scalar aveS = 0.0;

    for(int n=ignoreStepsInit_; n<Nmax; n++)
    {
        aveS = 0.0;
        for(int i = ignoreStepsInit_; i<=n; i++)
	{
	    r = 0.0;
	    for(int j = ignoreStepsInit_; j <= i - ignoreStepsDiag_; j++)
	    {
             Rij = 1 - base_.recM().recurrenceMatrix()[i][j];
	        if( Rij > r)
		    r = Rij;
	    }
	    aveS += r;
	}
	aveS /= (n - ignoreStepsInit_ + 1);
        aveSimilarity_[n - ignoreStepsInit_] = aveS;
    }
    
    OFstream aveSimilarityFile("aveSimilarityGeneral");
    aveSimilarityFile << aveSimilarity_;
}

void standardRecStatAnalysis::aveSimSpecial()
{
    Info << "Calculating special average similarity.\n" << endl;
    label Nmax = base_.recM().numRecFields();
    scalar r = 0.0;
    scalar Rij = 0.0;
    scalar aveS = 0.0;
    aveSimilarity_=0.0;

    for(int n=aveSimSpecialRefLength_; n<Nmax; n++)
    {
        aveS = 0.0;
        for(int i = ignoreStepsInit_; i<=aveSimSpecialRefLength_; i++)
        {
            r = 0.0;
            for(int j = i + ignoreStepsDiag_; j <= n; j++)
            {
             Rij = 1 - base_.recM().recurrenceMatrix()[i][j];
                if( Rij > r)
                    r = Rij;
            }
            aveS += r;
        }
        aveS /= (aveSimSpecialRefLength_ - ignoreStepsInit_ + 1);
        aveSimilarity_[n - ignoreStepsInit_] = aveS;
    }
    
    OFstream aveSimilarityFile("aveSimilaritySpecial");
    aveSimilarityFile << aveSimilarity_;

}

void standardRecStatAnalysis::mostDistinct()
{
    Info << "Calculating most distinct states.\n" << endl;
    label Nmax = base_.recM().numRecFields();
    scalar r = 0.0;
    scalar Rij = 0.0;

    for(int n=ignoreStepsInit_; n<Nmax; n++)
    {
        r = 1.0;
        for(int i = ignoreStepsInit_; i<=n; i++)
	{
	    for(int j = ignoreStepsInit_; j <= i - ignoreStepsDiag_; j++)
	    {
             Rij = 1 -  base_.recM().recurrenceMatrix()[i][j];
	        if( Rij < r)
		    r = Rij;
	    }
	}
        mostDistinct_[n - ignoreStepsInit_] = r;
    }
    
    OFstream mostDistinctFile("mostDistinct");
    mostDistinctFile << mostDistinct_;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
