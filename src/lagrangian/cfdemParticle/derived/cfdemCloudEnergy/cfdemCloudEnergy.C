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

#include "cfdemCloudEnergy.H"
#include "energyModel.H"
#include "massTransferModel.H"
#include "thermCondModel.H"
#include "diffCoeffModel.H"
#include "chemistryModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
cfdemCloudEnergy::cfdemCloudEnergy
(
    const fvMesh& mesh
)
:
    cfdemCloud(mesh),
    energyModels_(couplingProperties_.lookup("energyModels")),
    massTransferModels_(couplingProperties_.lookupOrDefault<wordList>("massTransferModels",wordList::null())),
    implicitEnergyModel_(false),
    implicitMassTransferModel_(false),
    chemistryModels_(couplingProperties_.lookup("chemistryModels")),
    energyModel_(nrEnergyModels()),
    massTransferModel_(nrMassTransferModels()),
    thermCondModel_
    (
        thermCondModel::New
        (
            couplingProperties_,
            *this
        )
    ),
    diffCoeffModel_
    (
        diffCoeffModel::New
        (
            couplingProperties_,
            *this
        )
    ),
    chemistryModel_(nrChemistryModels())
{
    forAll(energyModels_, modeli)
    {
        energyModel_.set
        (
            modeli,
            energyModel::New
            (
                couplingProperties_,
                *this,
                energyModels_[modeli]
            )
        );
    }

    forAll(massTransferModels_, modeli)
    {
        massTransferModel_.set
        (
            modeli,
            massTransferModel::New
            (
                couplingProperties_,
                *this,
                massTransferModels_[modeli]
            )
        );
    }

    forAll(chemistryModels_, modeli)
    {
        chemistryModel_.set
        (
            modeli,
            chemistryModel::New
            (
                couplingProperties_,
                *this,
                chemistryModels_[modeli]
            )
        );
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

cfdemCloudEnergy::~cfdemCloudEnergy()
{

}

// * * * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * * //

void cfdemCloudEnergy::calcEnergyContributions()
{
    forAll(energyModel_, modeli)
    {
        energyModel_[modeli].calcEnergyContribution();
    }
}

void cfdemCloudEnergy::calcMassContributions()
{
    forAll(massTransferModel_, modeli)
    {
        massTransferModel_[modeli].calcMassContribution();
    }
}

void cfdemCloudEnergy::speciesExecute()
{
    forAll(chemistryModel_, modeli)
    {
        chemistryModel_[modeli].execute();
    }
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
label cfdemCloudEnergy::nrEnergyModels() const
{
    return energyModels_.size();
}

label cfdemCloudEnergy::nrMassTransferModels() const
{
    return massTransferModels_.size();
}

int cfdemCloudEnergy::nrChemistryModels()
{
    return chemistryModels_.size();
}

bool cfdemCloudEnergy::implicitEnergyModel() const
{
    return implicitEnergyModel_;
}

bool cfdemCloudEnergy::implicitMassTransferModel() const
{
    return implicitMassTransferModel_;
}

const energyModel& cfdemCloudEnergy::energyM(int i)
{
    return energyModel_[i];
}

const massTransferModel& cfdemCloudEnergy::massTransferM(int i)
{
    return massTransferModel_[i];
}


const chemistryModel& cfdemCloudEnergy::chemistryM(int i)
{
    return chemistryModel_[i];
}

const thermCondModel& cfdemCloudEnergy::thermCondM()
{
    return thermCondModel_;
}

const diffCoeffModel& cfdemCloudEnergy::diffCoeffM()
{
    return diffCoeffModel_;
}

void cfdemCloudEnergy::energyContributions(volScalarField& Qsource)
{
    Qsource.primitiveFieldRef()=0.0;
    Qsource.boundaryFieldRef()=0.0;
    forAll(energyModel_, modeli)
    {
        energyM(modeli).addEnergyContribution(Qsource);
    }
}

void cfdemCloudEnergy::energyCoefficients(volScalarField& Qcoeff)
{
    Qcoeff.primitiveFieldRef()=0.0;
    Qcoeff.boundaryFieldRef()=0.0;
    forAll(energyModel_, modeli)
    {
        energyM(modeli).addEnergyCoefficient(Qcoeff);
    }
}

void cfdemCloudEnergy::massContributions(volScalarField& Sm)
{
    Sm.primitiveFieldRef()=0.0;
    Sm.boundaryFieldRef()=0.0;
    forAll(massTransferModel_, modeli)
    {
        massTransferM(modeli).addMassContribution(Sm);
    }
}

void cfdemCloudEnergy::massCoefficients(volScalarField& Smi)
{
    Smi.primitiveFieldRef()=0.0;
    Smi.boundaryFieldRef()=0.0;
    forAll(massTransferModel_, modeli)
    {
        massTransferM(modeli).addMassCoefficient(Smi);
    }
}

bool cfdemCloudEnergy::evolve
(
    volScalarField& alpha,
    volVectorField& Us,
    volVectorField& U
)
{
    if (cfdemCloud::evolve(alpha, Us, U))
    {
        // calc energy contributions including thermal conductivity
        // position 26 was already defined as Flow in clockModels and RhoPimpleChem solver.
        clockM().start(27,"calcEnergyContributions");
        if(verbose_) Info << "- calcEnergyContributions" << endl;
        calcEnergyContributions();
        thermCondModel_().calcThermCond();
        if(verbose_) Info << "calcEnergyContributions done." << endl;
        clockM().stop("calcEnergyContributions");

        // execute chemical model species
        clockM().start(28,"speciesExecute");
        if(verbose_) Info << "- speciesExecute()" << endl;
        speciesExecute();
        if(verbose_) Info << "speciesExecute done" << endl;
        clockM().stop("speciesExecute");

        // calculate mass contributions
        clockM().start(33,"calcMassContributions");
        if(verbose_) Info << "- calcMassContributions" << endl;
        calcMassContributions();
        if(verbose_) Info << "calcMassContributions done." << endl;
        clockM().stop("calcMassContributions");

        return true;
    }
    return false;
}

void cfdemCloudEnergy::postFlow()
{
    cfdemCloud::postFlow();
    forAll(energyModel_, modeli)
    {
        energyModel_[modeli].postFlow();
    }
    forAll(massTransferModel_, modeli)
    {
        massTransferModel_[modeli].postFlow();
    }
}

void cfdemCloudEnergy::solve()
{
    forAll(energyModel_, modeli)
    {
        energyModel_[modeli].solve();
    }
    forAll(massTransferModel_, modeli)
    {
        massTransferModel_[modeli].solve();
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
