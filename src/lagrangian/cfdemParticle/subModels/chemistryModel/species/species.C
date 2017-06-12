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
                        M.Efe Kinaci, JKU Linz, Austria

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
    verbose_(propsDict_.lookupOrDefault<bool>("verbose",false)),
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
    mod_spec_names_(speciesNames_.size()),
    Y_(speciesNames_.size()),                           //volumeScalarFields created in the ts folders
    concentrations_(speciesNames_.size(),NULL),         //the value of species concentration for every species
    changeOfSpeciesMass_(speciesNames_.size(),NULL),    //the values that are received from DEM with the name of Modified_+species name
    changeOfSpeciesMassFields_(speciesNames_.size()),   //the scalar fields generated with the values from Modified_+species names
    changeOfGasMassField_                               //the total change of Gas Mass field (when the Modified species
    (
        IOobject
        (
            "changeOfGasMassField",
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
         ),
        mesh_,
        dimensionedScalar("zero",dimMass/(dimVol*dimTime),0.0)
    ),
    tempFieldName_(propsDict_.lookup("tempFieldName")),
    tempField_(sm.mesh().lookupObject<volScalarField> (tempFieldName_)),
    partTempName_(propsDict_.lookup("partTempName")),
    partTemp_(NULL),
    densityFieldName_(propsDict_.lookup("densityFieldName")),
    rho_(sm.mesh().lookupObject<volScalarField> (densityFieldName_)),
    partRhoName_(propsDict_.lookup("partRhoName")),
    partRho_(NULL),
    voidfractionFieldName_(propsDict_.lookup("voidfractionFieldName")),
    voidfraction_(sm.mesh().lookupObject<volScalarField>(voidfractionFieldName_)),
    // total mole field
    totalMoleFieldName_(propsDict_.lookup("totalMoleFieldName")),
    N_(sm.mesh().lookupObject<volScalarField>(totalMoleFieldName_)),
    partMoleName_(propsDict_.lookup("partMoleName")),
    partN_(NULL),
    pressureFieldName_(propsDict_.lookup("pressureFieldName")),
    P_(sm.mesh().lookupObject<volScalarField>(pressureFieldName_)),
    partPName_(propsDict_.lookup("partPName")),
    partP_(NULL)
{
    Info << " Read species list from: " << specDict_.name() << endl;
    Info << " Reading species list: " << speciesNames_ << endl;

    for (int i=0; i<speciesNames_.size(); i++)
    {
        // Defining the Species volume scalar fields
        Info << " Looking up species fields \n " << speciesNames_[i] << endl;
        volScalarField& Y = const_cast<volScalarField&>
                (sm.mesh().lookupObject<volScalarField>(speciesNames_[i]));
        Y_.set(i, &Y);

         Info << "The concentration fields (Y_i): \n" << Y_[i].name() << endl;
        // define the modified species names
        mod_spec_names_[i] = "Modified_" + speciesNames_[i];

        // Check if mod species are correct
        Info << "Modified species names are: \n" << mod_spec_names_[i] << endl;

        // Create new volScalarFields for the changed values of the species mass fields
        changeOfSpeciesMassFields_.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                "ModSpeciesMassField_"+Y_[i].name(),
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
                ),
                mesh_,
                dimensionedScalar("0",dimMass/(dimVol*dimTime), 0)
            )
         );

        particleCloud_.checkCG(false);
    }

    allocateMyArrays();
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

species::~species()
{
    int nP_ = particleCloud_.numberOfParticles();

    particleCloud_.dataExchangeM().destroy(partTemp_,nP_);
    particleCloud_.dataExchangeM().destroy(partRho_,nP_);
    particleCloud_.dataExchangeM().destroy(partN_,nP_);
    particleCloud_.dataExchangeM().destroy(partP_,nP_);

    for (int i=0; i<speciesNames_.size(); i++)
    {
        particleCloud_.dataExchangeM().destroy(concentrations_[i],nP_);
        particleCloud_.dataExchangeM().destroy(changeOfSpeciesMass_[i],nP_);
    }
}

// * * * * * * * * * * * * * * * private Member Functions  * * * * * * * * * * * * * //
void species::allocateMyArrays() const
{
    double initVal=0.0;
    if (particleCloud_.dataExchangeM().maxNumberOfParticles() > 0)
    {
        // get memory for 2d arrays
        particleCloud_.dataExchangeM().allocateArray(partRho_,initVal,1,"nparticles");
        particleCloud_.dataExchangeM().allocateArray(partTemp_,initVal,1,"nparticles");
        particleCloud_.dataExchangeM().allocateArray(partN_,initVal,1,"nparticles");
        particleCloud_.dataExchangeM().allocateArray(partP_,initVal,1,"nparticles");

        for (int i=0; i<speciesNames_.size(); i++)
        {
            particleCloud_.dataExchangeM().allocateArray(concentrations_[i],initVal,1,"nparticles");
            particleCloud_.dataExchangeM().allocateArray(changeOfSpeciesMass_[i],initVal,1,"nparticles");

        }
    }
}


void species::reAllocMyArrays() const
{
    if (particleCloud_.numberOfParticlesChanged())
    {
        double initVal=0.0;
        particleCloud_.dataExchangeM().allocateArray(partRho_,initVal,1);
        particleCloud_.dataExchangeM().allocateArray(partTemp_,initVal,1);
        particleCloud_.dataExchangeM().allocateArray(partN_,initVal,1);
        particleCloud_.dataExchangeM().allocateArray(partP_,initVal,1);

        for (int i=0; i<speciesNames_.size(); i++)
        {
            particleCloud_.dataExchangeM().allocateArray(concentrations_[i],initVal,1);
            particleCloud_.dataExchangeM().allocateArray(changeOfSpeciesMass_[i],initVal,1);
        }
    }
}

// * * * * * * * * * * * * * * * * Member Fct  * * * * * * * * * * * * * * * //

void species::execute()
{
    // realloc the arrays
    reAllocMyArrays();

  // get Y_i, T, rho at particle positions, fill arrays with them and push to LIGGGHTS

    label  cellI=0;
    scalar Tfluid(0);
    scalar rhofluid(0);
    List<scalar> Yfluid_;
    scalar voidfraction(1);
    Yfluid_.setSize(speciesNames_.size());
    scalar Nfluid(0);
    scalar Pfluid(0);


    // defining interpolators for T, rho, voidfraction, N
    interpolationCellPoint <scalar> TInterpolator_(tempField_);
    interpolationCellPoint <scalar> rhoInterpolator_(rho_);
    interpolationCellPoint <scalar> voidfractionInterpolator_(voidfraction_);
    interpolationCellPoint <scalar> NInterpolator_(N_);
    interpolationCellPoint <scalar> PInterpolator_(P_);


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
                voidfraction        =   voidfractionInterpolator_.interpolate(position,cellI);
                Nfluid              =   NInterpolator_.interpolate(position,cellI);
                Pfluid              =   PInterpolator_.interpolate(position,cellI);
            }
            else
            {
                Tfluid          =   tempField_[cellI];
                rhofluid        =   rho_[cellI];
                voidfraction    =   voidfraction_[cellI];
                Nfluid          =   N_[cellI];
                Pfluid          =   P_[cellI];

                for (int i = 0; i<speciesNames_.size();i++)
                {
                    Yfluid_[i] = Y_[i][cellI];
                }
            }

            //fill arrays
            partTemp_[index][0] =   Tfluid;
            partRho_[index][0]  =   rhofluid*voidfraction;
            partN_[index][0]    =   Nfluid;
            partP_[index][0]    =   Pfluid;

            for (int i=0; i<speciesNames_.size();i++)
            {
                concentrations_[i][index][0]=Yfluid_[i];
                //concentrations_[i][index][0]=1*particleCloud_.mesh().time().value();
            }
        }

        if(verbose_ && index >=0 && index < 2)
        {
            for(int i =0; i<speciesNames_.size();i++)
            {
                Info << "Y_i = " << Y_[i].name() << endl;
                Info << "concentrations = " << concentrations_[i][index][0] << endl;
                Info << "partRho_[index][0] = " << partRho_[index][0] << endl;
                Info << "rhofluid =" << rhofluid << endl;
                Info << "Yfluid = " << Yfluid_[i] << endl;
                Info << "partTemp_[index][0] = " << partTemp_[index][0] << endl;
                Info << "Tfluid = " << Tfluid << endl  ;
                Info << "voidfraction =" << voidfraction << endl;
                Info << "N_" << N_ << endl;
            }
        }
    }

        // give DEM data
        particleCloud_.dataExchangeM().giveData(partTempName_, "scalar-atom", partTemp_);
        particleCloud_.dataExchangeM().giveData(partRhoName_, "scalar-atom", partRho_);
        particleCloud_.dataExchangeM().giveData(partMoleName_, "scalar-atom", partN_);
        particleCloud_.dataExchangeM().giveData(partPName_, "scalar-atom", partP_);

        for (int i=0; i<speciesNames_.size();i++)
        {
            particleCloud_.dataExchangeM().giveData(speciesNames_[i],"scalar-atom",concentrations_[i]);
        };

        Info << "give data done" << endl;

        // pull changeOfSpeciesMass_, transform onto fields changeOfSpeciesMassFields_, add them up on changeOfGasMassField_
        changeOfGasMassField_.primitiveFieldRef() = 0.0;
        changeOfGasMassField_.boundaryFieldRef() = 0.0;
        for (int i=0; i<speciesNames_.size();i++)
        {
            particleCloud_.dataExchangeM().getData(mod_spec_names_[i],"scalar-atom",changeOfSpeciesMass_[i]);
            changeOfSpeciesMassFields_[i].primitiveFieldRef() = 0.0;
            changeOfSpeciesMassFields_[i].boundaryFieldRef() = 0.0;
            particleCloud_.averagingM().setScalarSum
            (
                changeOfSpeciesMassFields_[i],
                changeOfSpeciesMass_[i],
                particleCloud_.particleWeights(),
                NULL
            );

            // take care for implementation in LIGGGHTS: species produced from particles defined positive
            changeOfSpeciesMassFields_[i].primitiveFieldRef() /= changeOfSpeciesMassFields_[i].mesh().V();
            changeOfSpeciesMassFields_[i].correctBoundaryConditions();
            changeOfGasMassField_ += changeOfSpeciesMassFields_[i];
            Info << "total conversion of species" << speciesNames_[i] << " = " << gSum(changeOfSpeciesMassFields_[i]*1.0*changeOfSpeciesMassFields_[i].mesh().V()) << endl;
        }
        Info << "get data done" << endl;
}

tmp<volScalarField> species::Smi (const label i) const
{
    return tmp<volScalarField> (changeOfSpeciesMassFields_[i]);
}

tmp<volScalarField> species::Sm () const
{
    return tmp<volScalarField> (changeOfGasMassField_);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
