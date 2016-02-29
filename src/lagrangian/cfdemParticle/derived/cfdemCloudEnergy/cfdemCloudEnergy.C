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
#include "thermCondModel.H"

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
    thermCondModel_
    (
        thermCondModel::New
        (
            couplingProperties_,
            *this
        )
    )
{
    energyModel_ = new autoPtr<energyModel>[nrEnergyModels()];
    for (int i=0;i<nrEnergyModels();i++)
    {
        energyModel_[i] = energyModel::New
        (
            couplingProperties_,
            *this,
            energyModels_[i]
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
    for (int i=0;i<nrEnergyModels();i++)
        energyModel_[i]().calcEnergyContribution();
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
int cfdemCloudEnergy::nrEnergyModels()
{
    return energyModels_.size();
}

const energyModel& cfdemCloudEnergy::energyM(int i)
{
    return energyModel_[i];
}

const thermCondModel& cfdemCloudEnergy::thermCondM()
{
    return thermCondModel_;
}

void cfdemCloudEnergy::energyContributions(volScalarField& Qsource)
{
    Qsource=0.0;
    for (int i=0;i<nrEnergyModels();i++)
        energyM(i).addEnergyContribution(Qsource);
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
        // calc energy contributions
        clockM().start(26,"calcEnergyContributions");
        if(verbose_) Info << "- calcEnergyContributions" << endl;
        calcEnergyContributions();
        if(verbose_) Info << "calcEnergyContributions done." << endl;
        clockM().stop("calcEnergyContributions");
	return true;
    }
    return false;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
