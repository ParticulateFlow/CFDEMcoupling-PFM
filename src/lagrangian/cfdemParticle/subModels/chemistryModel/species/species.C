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

#include "dataExchangeModel.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(species, 0);

addToRunTimeSelectionTable
(
        chemistryModel,
        species,
        dictionary
);

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
    interpolation_(propsDict_.lookupOrDefault<bool>("interpolation",false)),
    //speciesNames_(propsDict_.lookup("speciesNames"));
    //Y_(sm.mesh().lookup<volScalarField> (speciesNames_)),


    //concentrations_(),
    //changeOfSpeciesMass_(),

    //changeOfSpeciesMassFields_(sm.mesh().lookup<volScalarField> (changeOfSpeciesMassFields)),
    //changeOfGasMassFields_(sm.mesh().lookup<volScalarField> (changeOfGasMassField_)),

    tempFieldName_(propsDict_.lookupOrDefault<word>("tempFieldName","T")),
    tempField_(sm.mesh().lookupObject<volScalarField> (tempFieldName_)),
    partTempName_(propsDict_.lookup("partTempName")),
    partTemp_(NULL),

    densityFieldName_(propsDict_.lookupOrDefault<word>("densityFieldName","rho")),
    rho_(sm.mesh().lookupObject<volScalarField> (densityFieldName_)),
    partRhoName_(propsDict_.lookup("partRhoName")),
    partRho_(NULL),


    /*  voidfractionFieldName_(propsDict_.lookup("voidfractionFieldName")),
    voidfraction_(sm.mesh().lookupObject<volScalarField> (voidfractionFieldName_)),
    velFieldName_(propsDict_.lookup("velFieldName")),
    U_(sm.mesh().lookup<volVectorField> (velFieldName_)),*/

{
     allocateMyArrays();

}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

species::~species()
{
    delete partTemp_;
    delete partRho_;
}

// * * * * * * * * * * * * * * * private Member Functions  * * * * * * * * * * * * * //
void species::allocateMyArrays() const
{
  // could be implemented similarly as forcemodel LaEuScalarTemp

    // get memory for 2d arrays
    double initVal=0.0;
    particleCloud_.dataExchangeM().allocateArray(partRho_,initVal,1);
    particleCloud_.dataExchangeM().allocateArray(partTemp_,initVal,1);
    //particleCloud_.dataExchangeM().allocateArray(Y_,initVal,1);

}

// * * * * * * * * * * * * * * * * Member Fct  * * * * * * * * * * * * * * * //

void species::execute()
{
  // realloc the arrays
    allocateMyArrays();

  // get DEM data
    particleCloud_.dataExchangeM().getData(partTempName_,'scalar-atom',partTemp_);
    //particleCloud_.dataExchangeM().getData(partRhoName_,'scalar-atom',partRho_);

  // get Y_i, T, rho at particle positions, fill arrays with them and push to LIGGGHTS

    scalar Tfluid(0);
    label  cellI=0;
    scalar rho(0);
    scalar Yfluid(0);

    interpolationCellPoint<scalar> TInterpolator_(tempField_);
    interpolationCellPoint<scalar> rhoInterpolator_(rho_);
    itnerpolationCellPoint<scalar> YInterpolator_(Y_);

    for (int index=0; index<particleCloud_.numberOfParticles(); index++)
    {
        cellI=particleCloud_.cellIDs()[index][0];
        if (cellI >=0)
        {
            if(interpolation_)
          //if(chemistryM(0).interpolation())
            {
                vector position = particleCloud_.position(index);
                Tfluid = TInterpolator_.interpolate(poistion,cellI);
                rho=rhoInterpolator_.interpolate(position,cellI);
                Yfluid=YInterpolator_.interpolate(position,cellI);
            }
            else
            {
                Tfluid = tempField_[cellI];
                rho=rho_[cellI];
                Yfluid=Y_[cellI];
            }

            // calculate T, rho, Y_i

            //fill arrays
            partTemp_[index][0]=;
            partRho_[index][0]=;


        }
    }




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
