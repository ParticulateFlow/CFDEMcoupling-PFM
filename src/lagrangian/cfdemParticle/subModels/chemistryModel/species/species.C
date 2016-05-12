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
#include "IFstream.H"



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
    mesh_(sm.mesh()),
    // define a file name in the coupling properties that contains the species
    specDict_
    (
        IFstream
        (
            fileName(propsDict_.lookup("ChemistryFile")).expand()
        )()
    ),
    // create a list from the Species table in the specified species dictionary
    speciesNames_(specDict_.lookup("species")),
    Y_(speciesNames_.size()),
    concentrations_(speciesNames_.size(),NULL),
    changeOfSpeciesMass_(speciesNames_.size(),NULL),
    changeOfSpeciesMassFields_(speciesNames_.size()),
    changeOfGasMassField_
    (
        IOobject
        (
            "changeOfGasMassField_",
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
         ),
        mesh_,
        dimensionedScalar("zero",dimensionSet(0,0,0,0,0,0,0),0.0)
    ),
    tempFieldName_(propsDict_.lookupOrDefault<word>("tempFieldName","T")),
    tempField_(sm.mesh().lookupObject<volScalarField> (tempFieldName_)),
    partTempName_(propsDict_.lookup("partTempName")),
    partTemp_(NULL),
    densityFieldName_(propsDict_.lookupOrDefault<word>("densityFieldName","rho")),
    rho_(sm.mesh().lookupObject<volScalarField> (densityFieldName_)),
    partRhoName_(propsDict_.lookup("partRhoName")),
    partRho_(NULL)
  // voidfraction and velocity fields can be included by wish
  /*  voidfractionFieldName_(propsDict_.lookup("voidfractionFieldName")),
      voidfraction_(sm.mesh().lookupObject<volScalarField> (voidfractionFieldName_)),
      velFieldName_(propsDict_.lookup("velFieldName")),
      U_(sm.mesh().lookup<volVectorField> (velFieldName_)),*/

{ 
    Info << " Read species list from: " << specDict_.name() << endl;
    Info << " Reading species list: " << speciesNames_ << endl;

    for (int i=0; i<speciesNames_.size(); i++)
    {
        // Defining the Species volume scalar fields
        Info << " Looking up species fields " << speciesNames_[i] << endl;
        volScalarField& Y = const_cast<volScalarField&>
                (sm.mesh().lookupObject<volScalarField>(speciesNames_[i]));
        Y_.set(i, &Y);
        // Create new volScalarFields for the changed values of the species mass fields
        changeOfSpeciesMassFields_.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                "New"+Y_[i].name(),
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
                ),
                mesh_,
                dimensionedScalar("0",mesh_.lookupObject<volScalarField>(speciesNames_[i]).dimensions(), 0)
            )
         );
    }
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

    for (int i=0; i<speciesNames_.size(); i++)
    {
        particleCloud_.dataExchangeM().allocateArray(concentrations_[i],initVal,1);
        particleCloud_.dataExchangeM().allocateArray(changeOfSpeciesMass_[i],initVal,1);
    }
}

// * * * * * * * * * * * * * * * * Member Fct  * * * * * * * * * * * * * * * //

void species::execute()
{
  // realloc the arrays
    allocateMyArrays();

  // get Y_i, T, rho at particle positions, fill arrays with them and push to LIGGGHTS

    label  cellI=0;
    scalar Tfluid(0);
    scalar rhofluid(0);
    List<scalar> Yfluid_;
    Yfluid_.setSize(speciesNames_.size());
    List<scalar> changedField_;
    changedField_.setSize(speciesNames_.size());

    // defining interpolators for T, rho
    interpolationCellPoint <scalar> TInterpolator_(tempField_);
    interpolationCellPoint <scalar> rhoInterpolator_(rho_);

    for (int index=0; index<particleCloud_.numberOfParticles(); index++)
    {
        cellI=particleCloud_.cellIDs()[index][0];
        if (cellI >=0)
        {
            if(interpolation_)
            {
                vector position     =   particleCloud_.position(index);
                Tfluid              =   TInterpolator_.interpolate(position,cellI);
                rhofluid            =   rhoInterpolator_.interpolate(position,cellI);
            }
            else
            {
                Tfluid = tempField_[cellI];
                rhofluid=rho_[cellI];
                for (int i=0; i<speciesNames_.size();i++)
            }

            //fill arrays
            partTemp_[index][0]=Tfluid;
            partRho_[index][0]=rhofluid;
            for (int i=0; i<speciesNames_.size();i++)
            {
	        Yfluid_[i] = Y_[i][cellI];
                concentrations_[i][index][0]=Yfluid_[i];
            }
        }
    }

    // give DEM data
    particleCloud_.dataExchangeM().giveData(partTempName_, "scalar-atom", partTemp_);
    particleCloud_.dataExchangeM().giveData(partRhoName_, "scalar-atom", partRho_);
    for (int i=0; i<speciesNames_.size();i++)
    {
        particleCloud_.dataExchangeM().giveData(Y_[i].name(),"scalar-atom",concentrations_[i]);
    };

  // pull changeOfSpeciesMass_, transform onto fields changeOfSpeciesMassFields_, add them up on changeOfGasMassField_
  changeOfGasMassField_.internalField() = 0.0;
  changeOfGasMassField_.boundaryField() = 0.0;
  for (int i=0; i<speciesNames_.size();i++)
  {
      particleCloud_.dataExchangeM().getData(Y_[i].name(),"scalar-atom",changeOfSpeciesMass_[i]);

      changeOfSpeciesMassFields_[i].internalField() = 0.0;
      changeOfSpeciesMassFields_[i].boundaryField() = 0.0;
      particleCloud_.averagingM().setScalarSum
      (
        changeOfSpeciesMassFields_[i],
        changeOfSpeciesMass_[i],
        particleCloud_.particleWeights(),
        NULL
      );
      // take care for implementation in LIGGGHTS: species produced from particles defined positive
      changeOfSpeciesMassFields_[i].internalField() /= changeOfSpeciesMassFields_[i].mesh().V();
      changeOfSpeciesMassFields_[i].correctBoundaryConditions();   
      changeOfGasMassField_ += changeOfSpeciesMassFields_[i];
      Info << "total conversion of species " << speciesNames_[i] << " = " << gSum(changeOfSpeciesMassFields_[i]*1.0*changeOfSpeciesMassFields_[i].mesh().V()) << endl; 
  }
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
