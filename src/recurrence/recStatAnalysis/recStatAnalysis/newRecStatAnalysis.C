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
#include "recStatAnalysis.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

autoPtr<recStatAnalysis> recStatAnalysis::New
(
    const dictionary& dict,
    recBase& base
)
{
    word recStatAnalysisType
    (
        dict.lookupOrDefault<word>("recStatAnalysis","off")
    );

    Info << "Selecting recStatAnalysis "
         << recStatAnalysisType << endl;


    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(recStatAnalysisType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalError
            << "recStatAnalysis::New(const dictionary&, const spray&) : "
            << endl
            << "    unknown recStatAnalysisType type "
            << recStatAnalysisType
            << ", constructor not in hash table" << endl << endl
            << "    Valid recStatAnalysis types are :"
            << endl;
        Info << dictionaryConstructorTablePtr_->toc()
            << abort(FatalError);
    }

    return autoPtr<recStatAnalysis>(cstrIter()(dict,base));
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
