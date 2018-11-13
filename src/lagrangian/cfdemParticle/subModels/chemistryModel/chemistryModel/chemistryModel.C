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
                        M. Efe Kinaci, JKU Linz, Austria

\*---------------------------------------------------------------------------*/

#include "error.H"
#include "chemistryModel.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(chemistryModel, 0);

defineRunTimeSelectionTable(chemistryModel, dictionary);

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<volScalarField> chemistryModel::Smi(label i) const
{
    FatalError<<"the solver calls for Smi()\n"
              <<"please activate 'speciesModel' in 'chemistryModels'"
              <<abort(FatalError);

    tmp<volScalarField> tsource;
    return tsource;
}

tmp<volScalarField> chemistryModel::Sm() const
{
    FatalError<<"the solver calls for Sm()\n"
              <<"please activate 'speciesModel' in 'chemistryModels'"
              <<abort(FatalError);

    tmp<volScalarField> tsource;
    return tsource;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
chemistryModel::chemistryModel
(
    const dictionary& dict,
    cfdemCloudEnergy& sm
)
:
    dict_(dict),
    particleCloud_(sm)
{}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

chemistryModel::~chemistryModel()
{}

// * * * * * * * * * * * * * * * * Member Fct  * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
