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

#ifndef autocorrelation_H
#define autocorrelation_H

#include "recStatAnalysis.H"
#include "scalarList.H"
#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class autocorrelation Declaration
\*---------------------------------------------------------------------------*/

class autocorrelation
:
    public recStatAnalysis
{
protected:


public:

    //- Runtime type information
    TypeName("autocorrelation");


    // Constructors

        //- Construct from components
        autocorrelation
        (
            const dictionary& dict,
            recBase& base
        );

    // Destructor

        ~autocorrelation();

    void init();

    void statistics();


private:

    dictionary propsDict_;

    word fieldname_;

    word fieldtype_;

    label delaySteps_;

    label refCell_;

    vector refPoint_;

    bool normalize_;

    PtrList<volScalarField> scalarFieldList_;

    PtrList<volVectorField> vectorFieldList_;

    volScalarField autoCorrField_;

    bool suppressMatrixAndPath_;

    void autocorr();

    scalar autocorrSummand(label, label, label, label);

    volScalarField meanProd();

    bool suppressMatrixAndPath() {return suppressMatrixAndPath_;}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
