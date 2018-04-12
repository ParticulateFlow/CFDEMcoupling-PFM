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
#include "heatTransferGunnPartField.H"
#include "addToRunTimeSelectionTable.H"
#include "thermCondModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(heatTransferGunnPartField, 0);

addToRunTimeSelectionTable(energyModel, heatTransferGunnPartField, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
heatTransferGunnPartField::heatTransferGunnPartField
(
    const dictionary& dict,
    cfdemCloudEnergy& sm
)
:
    heatTransferGunn(dict,sm),
    fvOptions(fv::options::New(sm.mesh())),
    thermCondModel_
    (
        thermCondModel::New
        (
            propsDict_,
            sm
        )
    )
{
    if(!implicit_)
    {
        FatalError << "heatTransferGunnPartField requires implicit heat transfer treatment." << abort(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

heatTransferGunnPartField::~heatTransferGunnPartField()
{
}

// * * * * * * * * * * * * * * * private Member Functions  * * * * * * * * * * * * * //
// * * * * * * * * * * * * * * * * Member Fct  * * * * * * * * * * * * * * * //
void heatTransferGunnPartField::calcEnergyContribution()
{
    heatTransferGunn::calcEnergyContribution();
    // if heat sources in particles present, pull them here
}

void heatTransferGunnPartField::addEnergyContribution(volScalarField& Qsource) const
{
    Qsource -= QPartFluidCoeff_*partTempField_;
}

void heatTransferGunnPartField::giveData()
{
    particleCloud_.dataExchangeM().giveData(partTempName_,"scalar-atom", partTemp_);
}

void heatTransferGunnPartField::postFlow()
{
    label cellI;
    scalar Tpart(0.0);
    interpolationCellPoint<scalar> partTInterpolator_(partTempField_);

    for(int index = 0;index < particleCloud_.numberOfParticles(); ++index)
    {
            cellI = particleCloud_.cellIDs()[index][0];
            if(cellI >= 0)
            {
                if(interpolation_)
                {
                    vector position = particleCloud_.position(index);
                    Tpart = partTInterpolator_.interpolate(position,cellI);
                }
                else
                {
                    Tpart = partTempField_[cellI];
                }
                partTemp_[index][0] = Tpart;
            }
    }

    giveData();
}

void heatTransferGunnPartField::solve()
{
//use gasT and QCoeff to solve particle temp eqn
    volScalarField Qsource = QPartFluidCoeff_*tempField_;
    volScalarField alphaP = 1.0 - voidfraction_;
    volScalarField thCond = thermCondModel_().thermCond();

    fvScalarMatrix partTEqn
    (
        Qsource
      - fvm::Sp(QPartFluidCoeff_, partTempField_)
      - fvm::laplacian(alphaP*thCond,partTempField_)
     ==
        fvOptions(rho_, partTempField_)
    );
  // if transient add time derivative - need particle density and specific heat fields
  // if sources activated add sources
  // if convection activated add convection

    partTEqn.relax();

    fvOptions.constrain(partTEqn);

    partTEqn.solve();

    fvOptions.correct(partTempField_);

    Info<< "partT max/min : " << max(partTempField_).value() << " " << min(partTempField_).value() << endl;

}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
