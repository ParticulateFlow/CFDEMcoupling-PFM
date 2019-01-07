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
#include "heatTransferGranConduction.H"
#include "addToRunTimeSelectionTable.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(heatTransferGranConduction, 0);

addToRunTimeSelectionTable(energyModel, heatTransferGranConduction, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
heatTransferGranConduction::heatTransferGranConduction
(
    const dictionary& dict,
    cfdemCloudEnergy& sm
)
:
    energyModel(dict,sm),
    propsDict_(dict.subDict(typeName + "Props")),
    multiTypes_(false),
    verbose_(propsDict_.lookupOrDefault<bool>("verbose",false)),
    QPartPartName_(propsDict_.lookupOrDefault<word>("QPartPartName","QPartPart")),
    QPartPart_
    (   IOobject
        (
            QPartPartName_,
            sm.mesh().time().timeName(),
            sm.mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        sm.mesh(),
        dimensionedScalar("zero", dimensionSet(1,-1,-3,0,0,0,0), 0.0)
    ),
    partEffThermCondField_
    (
        IOobject
        (
            "partEffThermCondField",
            sm.mesh().time().timeName(),
            sm.mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        sm.mesh(),
        dimensionedScalar("one", dimensionSet(1, 1, -3, -1,0,0,0), 1.0),
	"zeroGradient"
    ),
    partThermCondField_
    (
        IOobject
        (
            "partThermCondField",
            sm.mesh().time().timeName(),
            sm.mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        sm.mesh(),
        dimensionedScalar("one", dimensionSet(1, 1, -3, -1,0,0,0), 1.0),
	"zeroGradient"
    ),
    partTempField_(sm.mesh().lookupObject<volScalarField>("partTemp")),
    voidfractionFieldName_(propsDict_.lookupOrDefault<word>("voidfractionFieldName","voidfraction")),
    voidfraction_(sm.mesh().lookupObject<volScalarField> (voidfractionFieldName_)),
    partHeatFluxName_(propsDict_.lookupOrDefault<word>("partHeatFluxName","conductiveHeatFlux")),
    partHeatFlux_(NULL),
    typePartThermCond_(propsDict_.lookupOrDefault<scalarList>("thermalConductivities",scalarList(1,-1.0))),
    partThermCond_(NULL)
{
    allocateMyArrays();

    if (typePartThermCond_[0] < 0.0)
    {
        FatalError << "heatTransferGranConduction: provide list of thermal conductivities." << abort(FatalError);
    }

    if (typePartThermCond_.size() > 1)
    {
        multiTypes_ = true;
    }

    if (multiTypes_ = true && !particleCloud_.getParticleTypes())
    {
        FatalError << "heatTransferGranConduction needs data for more than one type, but types are not communicated." << abort(FatalError);
    }

    if (verbose_)
    {
        QPartPart_.writeOpt() = IOobject::AUTO_WRITE;
        QPartPart_.write();
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

heatTransferGranConduction::~heatTransferGranConduction()
{
    particleCloud_.dataExchangeM().destroy(partHeatFlux_,1);
    particleCloud_.dataExchangeM().destroy(partThermCond_,1);
}

// * * * * * * * * * * * * * * * private Member Functions  * * * * * * * * * * * * * //
void heatTransferGranConduction::allocateMyArrays() const
{
    // get memory for 2d arrays
    double initVal=0.0;
    particleCloud_.dataExchangeM().allocateArray(partHeatFlux_,initVal,1);
    particleCloud_.dataExchangeM().allocateArray(partThermCond_,initVal,1);
}
// * * * * * * * * * * * * * * * * Member Fct  * * * * * * * * * * * * * * * //

void heatTransferGranConduction::calcEnergyContribution()
{
   // realloc the arrays
    allocateMyArrays();

    calcPartEffThermCond();

    QPartPart_ = fvc::laplacian(partEffThermCondField_,partTempField_);

    label cellI=0;
    scalar partVolume(0);
    scalar QPartPart(0);
    scalar voidfraction(1);

    for(int index = 0;index < particleCloud_.numberOfParticles(); ++index)
    {
        cellI = particleCloud_.cellIDs()[index][0];
        if(cellI >= 0)
        {
            voidfraction = voidfraction_[cellI];
            if (voidfraction < 0.01)
                voidfraction = 0.01;

            partVolume = particleCloud_.particleVolume(index);
            QPartPart = QPartPart_[cellI];

            heatFlux(index, partVolume, voidfraction, QPartPart);
        }
    }
}

void heatTransferGranConduction::calcPartEffThermCond()
{
    calcPartThermCond();

    scalar volFrac = 0.0;

    forAll(partEffThermCondField_, cellI)
    {
        volFrac = 1.0 - voidfraction_[cellI];
        if (volFrac < 0.334)
        {
            partEffThermCondField_[cellI] = 0.0;
        }
        else
        {
            partEffThermCondField_[cellI] = 0.5 * (3*volFrac - 1) * partThermCondField_[cellI];
        }
    }
}

void heatTransferGranConduction::calcPartThermCond()
{
    label cellI=0;
    label partType = 1;
    for(int index = 0;index < particleCloud_.numberOfParticles(); ++index)
    {
            cellI = particleCloud_.cellIDs()[index][0];
            if (cellI >= 0)
            {
                if (multiTypes_)
                {
                    partType = particleCloud_.particleType(index);
                }
                // LIGGGGHTS counts types 1, 2, ..., C++ array starts at 0
                partThermCond_[index][0] = typePartThermCond_[partType - 1];
            }
    }

    partThermCondField_.primitiveFieldRef() = 0.0;
    particleCloud_.averagingM().resetWeightFields();
    particleCloud_.averagingM().setScalarAverage
    (
        partThermCondField_,
        partThermCond_,
        particleCloud_.particleWeights(),
        particleCloud_.averagingM().UsWeightField(),
        NULL
    );
}

void heatTransferGranConduction::heatFlux(label index, scalar vol, scalar voidfraction, scalar QPartPart)
{
        partHeatFlux_[index][0] = vol / (1.0 - voidfraction) * QPartPart;
}

void heatTransferGranConduction::giveData()
{
    Info << "total conductive particle-particle heat flux [W] (Eulerian) = " << gSum(QPartPart_*1.0*QPartPart_.mesh().V()) << endl;

    particleCloud_.dataExchangeM().giveData(partHeatFluxName_,"scalar-atom", partHeatFlux_);
}

void heatTransferGranConduction::postFlow()
{
    giveData();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

