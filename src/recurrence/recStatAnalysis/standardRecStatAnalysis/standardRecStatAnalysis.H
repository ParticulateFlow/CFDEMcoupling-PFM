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

#ifndef standardRecStatAnalysis_H
#define standardRecStatAnalysis_H

#include "recStatAnalysis.H"
#include "scalarList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class standardRecStatAnalysis Declaration
\*---------------------------------------------------------------------------*/

class standardRecStatAnalysis
:
    public recStatAnalysis
{
protected:
  
    

public:

    //- Runtime type information
    TypeName("standardRecStatAnalysis");


    // Constructors

        //- Construct from components
        standardRecStatAnalysis
        (
            const dictionary& dict,
            recBase& base
        );

    // Destructor

        ~standardRecStatAnalysis();

    void init();
	
    void statistics();
    
	
private:
  
    dictionary propsDict_;

    label aveSimSpecialRefLength_;
    
    label ignoreStepsInit_;
    
    label ignoreStepsDiag_;
    
    scalarList aveSimilarity_;
    
    scalarList correlationSum_;
    
    scalarList mostDistinct_;
    
    void aveSimGeneral();

    void aveSimSpecial();
    
    void mostDistinct();

    
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
