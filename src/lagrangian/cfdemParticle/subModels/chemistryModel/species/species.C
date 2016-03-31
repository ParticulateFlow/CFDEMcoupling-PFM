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
#include "species.H"
#include "addToRunTimeSelectionTable.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(species, 0);

addToRunTimeSelectionTable(chemistryModel, species, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
species::species
(
    const dictionary& dict,
    cfdemCloudEnergy& sm
)
:
    chemistryModel(dict,sm),
    propsDict_(dict.subDict(typeName + "Props")),
    interpolation_(propsDict_.lookupOrDefault<bool>("interpolation",false))
{
     allocateMyArrays();

    
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

species::~species()
{

}

// * * * * * * * * * * * * * * * private Member Functions  * * * * * * * * * * * * * //
void species::allocateMyArrays() const
{
  // could be implemented similarly as forcemodel LaEuScalarTemp
}

// * * * * * * * * * * * * * * * * Member Fct  * * * * * * * * * * * * * * * //

void species::execute()
{
  
  // get Y_i, T, rho at particle positions, fill arrays with them and push to LIGGGHTS
  
  // pull changeOfSpeciesMass_, transform onto fields changeOfSpeciesMassFields_, add them up on changeOfGasMassField_
  
}

//tmp<Foam::fvScalarMatrix> species::Smi(const label i) const
//{
//    return tmp<fvScalarMatrix>(new fvScalarMatrix(changeOfSpeciesMassFields_[i], dimMass/dimTime)); 
//}

//tmp<Foam::fvScalarMatrix> species::Sm() const
//{
//    return tmp<fvScalarMatrix>(new fvScalarMatrix(changeOfGasMassField_, dimMass/dimTime)); 
//}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
