/*---------------------------------------------------------------------------*\
License

    This is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This code is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with this code.  If not, see <http://www.gnu.org/licenses/>.

    Copyright (C) 2015- Thomas Lichtenegger, JKU Linz, Austria

\*---------------------------------------------------------------------------*/

#include "error.H"

#include "chemistryModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

autoPtr<chemistryModel> chemistryModel::New
(
    const dictionary& dict,
    cfdemCloudEnergy& sm,
    word chemistryModelType
)
{
    Info<< "Selecting chemistryModel "
         << chemistryModelType << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(chemistryModelType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalError
            << "chemistryModel::New(const dictionary&, cfdemCloudEnergy&) : "
            << endl
            << "    unknown chemistryModelType type "
            << chemistryModelType
            << ", constructor not in hash table" << endl << endl
            << "    Valid chemistryModel types are :"
            << endl;
        Info<< dictionaryConstructorTablePtr_->toc()
            << abort(FatalError);
    }

    return autoPtr<chemistryModel>(cstrIter()(dict,sm));
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
