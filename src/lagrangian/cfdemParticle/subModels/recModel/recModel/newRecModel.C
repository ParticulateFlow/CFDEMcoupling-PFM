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
#include "standardSearch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

autoPtr<recModel> recModel::New
(
    const dictionary& dict,
    cfdemCloud& sm
)
{
    word recModelType
    (
        dict.lookup("recModel")
    );

    Info<< "Selecting recModel "
         << recModelType << endl;


    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(recModelType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalError
            << "recModel::New(const dictionary&, const spray&) : "
            << endl
            << "    unknown recModelType type "
            << recModelType
            << ", constructor not in hash table" << endl << endl
            << "    Valid recModel types are :"
            << endl;
        Info<< dictionaryConstructorTablePtr_->toc()
            << abort(FatalError);
    }

    return autoPtr<recModel>(cstrIter()(dict,sm));
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
