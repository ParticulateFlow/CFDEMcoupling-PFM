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
    X_(speciesNames_.size()),                           //volumeScalarFields of molarFractions
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
    tempFieldName_(propsDict_.lookupOrDefault<word>("tempFieldName","T")),
    tempField_(sm.mesh().lookupObject<volScalarField> (tempFieldName_)),
    partTempName_(propsDict_.lookupOrDefault<word>("partTempName","partTemp")),
    densityFieldName_(propsDict_.lookupOrDefault<word>("densityFieldName","rho")),
    rho_(sm.mesh().lookupObject<volScalarField> (densityFieldName_)),
    partRhoName_(propsDict_.lookupOrDefault<word>("partRhoName","partRho")),
    voidfractionFieldName_(propsDict_.lookupOrDefault<word>("voidfractionFieldName","voidfraction")),
    voidfraction_(sm.mesh().lookupObject<volScalarField>(voidfractionFieldName_)),
    // total mole field
    molarConcFieldName_(propsDict_.lookupOrDefault<word>("totalMoleFieldName","molarConc")),
    molarConc_(sm.mesh().lookupObject<volScalarField>(molarConcFieldName_)),
    partMolarConcName_(propsDict_.lookupOrDefault<word>("partMoleName","partMolarConc")),
    loopCounter_(-1),
    Nevery_(propsDict_.lookupOrDefault<label>("Nevery",1)),
    massSourceCurr_(0.0),
    massSourceTot_(0.0),
    initialized_(false)
{
    particleCloud_.checkCG(false);
    particleCloud_.registerParticleProperty<double**>(partTempName_);
    particleCloud_.registerParticleProperty<double**>(partRhoName_);
    particleCloud_.registerParticleProperty<double**>(partMolarConcName_);

    for (int i=0; i<speciesNames_.size(); i++)
    {
        particleCloud_.registerParticleProperty<double**>("X_"+speciesNames_[i]);
        particleCloud_.registerParticleProperty<double**>("Modified_"+speciesNames_[i]);
    }
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

species::~species()
{
}

// * * * * * * * * * * * * * * * private Member Functions  * * * * * * * * * * * * * //

void species::reAllocMyArrays() const
{
    double initVal=0.0;
    double**& partRho_ = particleCloud_.getParticlePropertyRef<double**>(partRhoName_);
    double**& partTemp_ = particleCloud_.getParticlePropertyRef<double**>(partTempName_);
    double**& partMolarConc_ = particleCloud_.getParticlePropertyRef<double**>(partMolarConcName_);

    particleCloud_.dataExchangeM().allocateArray(partRho_,initVal,1);
    particleCloud_.dataExchangeM().allocateArray(partTemp_,initVal,1);
    particleCloud_.dataExchangeM().allocateArray(partMolarConc_,initVal,1);

    for (int i=0; i<speciesNames_.size(); i++)
    {
        double**& molarFractions_ = particleCloud_.getParticlePropertyRef<double**>("X_"+speciesNames_[i]);
        double**& changeOfSpeciesMass_ = particleCloud_.getParticlePropertyRef<double**>("Modified_"+speciesNames_[i]);
        particleCloud_.dataExchangeM().allocateArray(molarFractions_,initVal,1);
        particleCloud_.dataExchangeM().allocateArray(changeOfSpeciesMass_,initVal,1);
    }
}

void species::init()
{
    if(verbose_)
    {
        Info << " Read species list from: " << specDict_.name() << endl;
        Info << " Reading species list: " << speciesNames_ << endl;
    }

    for (int i=0; i<speciesNames_.size(); i++)
    {
        // Define the Species volume scalar fields
        volScalarField& X = const_cast<volScalarField&>
                (mesh_.lookupObject<volScalarField>("X_"+speciesNames_[i]));
        X_.set(i, &X);

        // define the modified species names
        mod_spec_names_[i] = "Modified_" + speciesNames_[i];

        if(verbose_)
        {
            Info << "The molar fraction fields (X_i): \n" << X_[i].name() << endl;
            // Check if mod species are correct
            Info << "Modified species names are: \n" << mod_spec_names_[i] << endl;
        }

        // Create new volScalarFields for the changed values of the species mass fields -- gas species source term
        changeOfSpeciesMassFields_.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                "TotalChangeOfMassField_"+speciesNames_[i],
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
                ),
                mesh_,
                dimensionedScalar("zero",dimMass/(dimVol*dimTime), 0.0)
            )
        );
    }
    initialized_ = true;
}

// * * * * * * * * * * * * * * * * Member Fct  * * * * * * * * * * * * * * * //

void species::execute()
{
    if(!initialized_)
    {
        init();
    }

    loopCounter_++;
    if (loopCounter_ % Nevery_ != 0)
    {
        return;
    }
    // realloc the arrays
    reAllocMyArrays();

    // get X_i, T, rho at particle positions
    label  cellI = 0;
    scalar Tfluid(0);
    scalar rhofluid(0);
    scalar voidfraction(1);
    scalar molarConcfluid(0);

    // defining interpolators for T, rho, voidfraction, molarConc
    interpolationCellPoint <scalar> TInterpolator_(tempField_);
    interpolationCellPoint <scalar> rhoInterpolator_(rho_);
    interpolationCellPoint <scalar> voidfractionInterpolator_(voidfraction_);
    interpolationCellPoint <scalar> molarConcInterpolator_(molarConc_);

    double**& partRho_ = particleCloud_.getParticlePropertyRef<double**>(partRhoName_);
    double**& partTemp_ = particleCloud_.getParticlePropertyRef<double**>(partTempName_);
    double**& partMolarConc_ = particleCloud_.getParticlePropertyRef<double**>(partMolarConcName_);

    for (int index=0; index<particleCloud_.numberOfParticles(); ++index)
    {
        cellI=particleCloud_.cellIDs()[index][0];
        if (cellI >= 0)
        {
            if(interpolation_)
            {
                vector position = particleCloud_.position(index);
                Tfluid          = TInterpolator_.interpolate(position,cellI);
                rhofluid        = rhoInterpolator_.interpolate(position,cellI);
                voidfraction    = voidfractionInterpolator_.interpolate(position,cellI);
                molarConcfluid  = molarConcInterpolator_.interpolate(position,cellI);
            }
            else
            {
                Tfluid          = tempField_[cellI];
                rhofluid        = rho_[cellI];
                voidfraction    = voidfraction_[cellI];
                molarConcfluid  = molarConc_[cellI];
            }

            partTemp_[index][0] = Tfluid;
            partRho_[index][0]  = rhofluid;
            partMolarConc_[index][0] = molarConcfluid;

            for (int i=0; i<speciesNames_.size();i++)
            {
                // attention for indices when not communicating all species
                double**& molarFractions_ = particleCloud_.getParticlePropertyRef<double**>("X_"+speciesNames_[i]);
                molarFractions_[index][0] = X_[i][cellI];
            }
        }
    }


    if(verbose_)
    {
        for(int i =0; i<speciesNames_.size(); i++)
        {
            double**& molarFractions_ = particleCloud_.getParticlePropertyRef<double**>("X_"+speciesNames_[i]);
            Info << "X_i = " << X_[i].name() << endl;
            Info << "molarFractions_= " << molarFractions_[0][0] << endl;
            Info << "partRho_[index][0] = " << partRho_[0][0] << endl;
            Info << "rhofluid = " << rhofluid << endl;
            Info << "partTemp_[index][0] = " << partTemp_[0][0] << endl;
            Info << "Tfluid = " << Tfluid << endl;
            Info << "voidfraction = " << voidfraction << endl;
        }
    }

    // give DEM data
    {
        particleCloud_.dataExchangeM().giveData(partTempName_, "scalar-atom", partTemp_);
        particleCloud_.dataExchangeM().giveData(partRhoName_,  "scalar-atom", partRho_);
        particleCloud_.dataExchangeM().giveData(partMolarConcName_, "scalar-atom", partMolarConc_);

        for (int i=0; i<speciesNames_.size();i++)
        {
            double**& molarFractions_ = particleCloud_.getParticlePropertyRef<double**>("X_"+speciesNames_[i]);
            particleCloud_.dataExchangeM().giveData("X_"+speciesNames_[i],"scalar-atom",molarFractions_);
        }

        if (verbose_) Info << "give data done" << endl;
    }

    // pull changeOfSpeciesMass_, transform onto fields changeOfSpeciesMassFields_, add them up on changeOfGasMassField_
    {
        scalar timestep = mesh_.time().deltaTValue();
        changeOfGasMassField_.primitiveFieldRef() = 0.0;
        changeOfGasMassField_.boundaryFieldRef() = 0.0;
        for (int i=0; i<speciesNames_.size();i++)
        {
            double**& changeOfSpeciesMass_ = particleCloud_.getParticlePropertyRef<double**>("Modified_"+speciesNames_[i]);
            changeOfSpeciesMassFields_[i].primitiveFieldRef() = 0.0;
            changeOfSpeciesMassFields_[i].boundaryFieldRef() = 0.0;

            particleCloud_.dataExchangeM().getData(mod_spec_names_[i],"scalar-atom",changeOfSpeciesMass_,particleCloud_.dataExchangeM().couplingInterval());

            if (verbose_) Info << "changeOfSpeciesMass received from DEM = " << changeOfSpeciesMass_[0][0] << endl;

            particleCloud_.averagingM().setScalarSumCentre
            (
                changeOfSpeciesMassFields_[i],
                changeOfSpeciesMass_,
                particleCloud_.particleWeights(),
                NULL
            );

            // take care for implementation in LIGGGHTS: species produced from particles defined positive
            // changeOf...Fields need to be mass per volume per timestep
            changeOfSpeciesMassFields_[i].primitiveFieldRef() /= (changeOfSpeciesMassFields_[i].mesh().V() * Nevery_ * timestep);
            changeOfSpeciesMassFields_[i].correctBoundaryConditions();
            changeOfGasMassField_ += changeOfSpeciesMassFields_[i];

            if (verbose_)
            {
                Info << "total conversion of species" << speciesNames_[i] << " = "
                     << gSum(changeOfSpeciesMassFields_[i]*1.0*changeOfSpeciesMassFields_[i].mesh().V() * Nevery_ * timestep) << endl;
            }
        }
        massSourceCurr_ = gSum(changeOfGasMassField_*1.0*changeOfGasMassField_.mesh().V() * Nevery_ * timestep);
        massSourceTot_ += massSourceCurr_;

        if (verbose_)
        {
            Info << "total conversion of mass:\n\tcurrent source = " << massSourceCurr_ << "\n\ttotal source = " << massSourceTot_ << "\n" << endl;
            Info << "get data done" << endl;
        }
    }
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
