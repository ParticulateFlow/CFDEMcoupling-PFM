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
    partCpField_
    (
        IOobject
        (
            "partCp",
            sm.mesh().time().timeName(),
            sm.mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        sm.mesh(),
        dimensionedScalar("zero", dimensionSet(0,2,-2,-1,0,0,0), 0.0),
        "zeroGradient"
    ),
    partRhoField_(sm.mesh().lookupObject<volScalarField>("partRho")),
    typeCp_(propsDict_.lookupOrDefault<scalarList>("specificHeatCapacities",scalarList(1,-1.0))),
    partCp_(NULL),
    pTMax_(dimensionedScalar("pTMax",dimensionSet(0,0,0,1,0,0,0), -1.0)),
    pTMin_(dimensionedScalar("pTMin",dimensionSet(0,0,0,1,0,0,0), -1.0)),
    thermCondModel_
    (
        thermCondModel::New
        (
            propsDict_,
            sm
        )
    ),
    fvOptions(fv::options::New(sm.mesh()))
{
    if (!implicit_)
    {
        FatalError << "heatTransferGunnPartField requires implicit heat transfer treatment." << abort(FatalError);
    }

    if (typeCp_[0] < 0.0)
    {
        FatalError << "heatTransferGunnPartField: provide list of specific heat capacities." << abort(FatalError);
    }

    if (propsDict_.found("pTMax"))
    {
        pTMax_.value()=scalar(readScalar(propsDict_.lookup("pTMax")));
    }

    if (propsDict_.found("pTMin"))
    {
        pTMin_.value()=scalar(readScalar(propsDict_.lookup("pTMin")));
    }

    partTempField_.writeOpt() = IOobject::AUTO_WRITE;

    allocateMyArrays();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

heatTransferGunnPartField::~heatTransferGunnPartField()
{
    particleCloud_.dataExchangeM().destroy(partCp_,1);
}

// * * * * * * * * * * * * * * * private Member Functions  * * * * * * * * * * * * * //

void heatTransferGunnPartField::allocateMyArrays() const
{
    double initVal=0.0;
    particleCloud_.dataExchangeM().allocateArray(partCp_,initVal,1);
}
// * * * * * * * * * * * * * * * * Member Fct  * * * * * * * * * * * * * * * //
void heatTransferGunnPartField::calcEnergyContribution()
{
    allocateMyArrays();
    heatTransferGunn::calcEnergyContribution();

    // if heat sources in particles present, pull them here

    // loop over all particles to fill partCp_ based on type
    label cellI=0;
    label partType = 0;
    for(int index = 0;index < particleCloud_.numberOfParticles(); ++index)
    {
            cellI = particleCloud_.cellIDs()[index][0];
            if(cellI >= 0)
            {
                partType = particleCloud_.particleType(index);
                // LIGGGGHTS counts types 1, 2, ..., C++ array starts at 0
                partCp_[index][0] = typeCp_[partType - 1];
            }
    }

    partCpField_.primitiveFieldRef() = 0.0;
    particleCloud_.averagingM().resetWeightFields();
    particleCloud_.averagingM().setScalarAverage
    (
        partCpField_,
        partCp_,
        particleCloud_.particleWeights(),
        particleCloud_.averagingM().UsWeightField(),
        NULL
    );

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

    particleCloud_.dataExchangeM().allocateArray(partTemp_,0.0,1);

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
    // for some weird reason, the particle-fluid heat transfer fields were defined with a negative sign

    volScalarField alphaP = 1.0 - voidfraction_;
    volScalarField correctedQPartFluidCoeff(QPartFluidCoeff_);
    // no heattransfer in empty cells -- for numerical stability, add small amount
    forAll(correctedQPartFluidCoeff,cellI)
    {
        if (-correctedQPartFluidCoeff[cellI] < SMALL) correctedQPartFluidCoeff[cellI] = -SMALL;
    }

    volScalarField Qsource = correctedQPartFluidCoeff*tempField_;
    volScalarField partCpEff = alphaP*partRhoField_*partCpField_;
    volScalarField thCondEff = alphaP*thermCondModel_().thermCond();
//    thCondEff.correctBoundaryConditions();



    fvScalarMatrix partTEqn
    (
     //   fvm::ddt(partCpEff, partTempField_)
     // + Qsource
        Qsource
      - fvm::Sp(correctedQPartFluidCoeff, partTempField_)
      - fvm::laplacian(thCondEff,partTempField_)
     ==
        fvOptions(partCpEff, partTempField_)
    );
  // if transient add time derivative - need particle density and specific heat fields
  // if sources activated add sources
  // if convection activated add convection

    partTEqn.relax();

    fvOptions.constrain(partTEqn);

    partTEqn.solve();

    partTempField_.relax();

    fvOptions.correct(partTempField_);

    if (pTMax_.value() > 0.0) partTempField_ = min(partTempField_, pTMax_);
    if (pTMin_.value() > 0.0) partTempField_ = max(partTempField_, pTMin_);


    Info<< "partT max/min : " << max(partTempField_).value() << " " << min(partTempField_).value() << endl;

}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
