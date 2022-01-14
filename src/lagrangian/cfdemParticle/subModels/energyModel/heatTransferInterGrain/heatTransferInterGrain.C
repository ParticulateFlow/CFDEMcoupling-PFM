
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
#include "heatTransferInterGrain.H"
#include "addToRunTimeSelectionTable.H"
#include "fvCFD.H"
#include "physicoChemicalConstants.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(heatTransferInterGrain, 0);

addToRunTimeSelectionTable(energyModel, heatTransferInterGrain, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
heatTransferInterGrain::heatTransferInterGrain
(
    const dictionary& dict,
    cfdemCloudEnergy& sm
)
:
    energyModel(dict,sm),
    propsDict_(dict.subDict(typeName + "Props")),
    multiTypes_(false),
    verbose_(propsDict_.lookupOrDefault<bool>("verbose",false)),
    implicit_(propsDict_.lookupOrDefault<bool>("implicit",false)),
    calcTotalHeatFlux_(propsDict_.lookupOrDefault<bool>("calcTotalHeatFlux",false)),
    radiativeHeatTransfer_(propsDict_.lookupOrDefault<bool>("radiativeHeatTransfer",false)),
    totalHeatFlux_(0.0),
    partTempName_(""),
    QPartPartName_(propsDict_.lookupOrDefault<word>("QPartPartName","QPartPart")),
    QPartPart_
    (
        IOobject
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
    partThermCapField_
    (
        IOobject
        (
            "partThermCapField",
            sm.mesh().time().timeName(),
            sm.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        sm.mesh(),
        dimensionedScalar("zero", dimensionSet(1,-1,-2,-1,0,0,0), 0.0),
        "zeroGradient"
    ),
    partThermRadField_
    (
        IOobject
        (
            "partThermRadField",
            sm.mesh().time().timeName(),
            sm.mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        sm.mesh(),
        dimensionedScalar("zero", dimensionSet(1, 1, -3, -1,0,0,0), 0.0),
        "zeroGradient"
    ),
    partTempField_(sm.mesh().lookupObject<volScalarField>("partTemp")),
    prescribedVoidfractionFieldName_(propsDict_.lookupOrDefault<word>("prescribedVoidfractionFieldName","voidfraction")),
    prescribedVoidfraction_(sm.mesh().lookupObject<volScalarField> (prescribedVoidfractionFieldName_)),
    voidfractionFieldName_(propsDict_.lookupOrDefault<word>("voidfractionFieldName","voidfraction")),
    voidfraction_(sm.mesh().lookupObject<volScalarField> (voidfractionFieldName_)),
    typeCG_(propsDict_.lookupOrDefault<scalarList>("coarseGrainingFactors",scalarList(1,1.0))),
    maxTypeCG_(typeCG_.size()),
    partHeatFluxName_(propsDict_.lookupOrDefault<word>("partHeatFluxName","conductiveHeatFlux")),
    typePartThermCond_(propsDict_.lookupOrDefault<scalarList>("thermalConductivities",scalarList(1,-1.0))),
    partThermCondRegName_(typeName + "partThermCond"),
    typePartThermCap_(propsDict_.lookupOrDefault<scalarList>("thermalCapacities",scalarList(1,-1.0))),
    partThermCapRegName_(typeName + "partThermCap"),
    partThermRadRegName_(typeName + "partThermRad"),
    typePartEmissivity_(propsDict_.lookupOrDefault<scalarList>("thermalEmissivities",scalarList(1,-1.0))),
    kMax_(propsDict_.lookupOrDefault<scalar>("kMax",-1.0))

{
    particleCloud_.registerParticleProperty<double**>(partHeatFluxName_,1);
    particleCloud_.registerParticleProperty<double**>(partThermCondRegName_,1);
    if (radiativeHeatTransfer_)
    {
        particleCloud_.registerParticleProperty<double**>(partThermRadRegName_,1);
        partTempName_ = word(propsDict_.lookup("partTempName"));
    }

    if (typePartThermCond_[0] < 0.0)
    {
        FatalError << "heatTransferInterGrain: provide list of thermal conductivities." << abort(FatalError);
    }

    if (typePartEmissivity_[0] < 0.0 && radiativeHeatTransfer_)
    {
        FatalError << "heatTransferGranRadiation: provide list of thermal emissivities." << abort(FatalError);
    }

    if (typePartEmissivity_[0] > 0.0 && !radiativeHeatTransfer_)
    {
        FatalError << "heatTransferGranRadiation: thermal emissivities provided but calculation deactivated." << abort(FatalError);
    }

    if (typePartThermCond_.size() > 1)
    {
        multiTypes_ = true;
    }

    if (multiTypes_ && !particleCloud_.getParticleTypes())
    {
        FatalError << "heatTransferInterGrain needs data for more than one type, but types are not communicated." << abort(FatalError);
    }

    if (verbose_)
    {
        QPartPart_.writeOpt() = IOobject::AUTO_WRITE;
        partEffThermCondField_.writeOpt() = IOobject::AUTO_WRITE;
        partThermCondField_.writeOpt() = IOobject::AUTO_WRITE;
        partThermCapField_.writeOpt() = IOobject::AUTO_WRITE;
        partThermRadField_.writeOpt() = IOobject::AUTO_WRITE;
        QPartPart_.write();
        partEffThermCondField_.write();
        partThermCondField_.write();
        partThermCapField_.write();
        partThermRadField_.write();
    }

    if (implicit_)
    {
        particleCloud_.registerParticleProperty<double**>(partThermCapRegName_,1);
        if (typePartThermCap_[0] < 0.0)
        {
            FatalError << "heatTransferInterGrain: provide list of thermal capacities." << abort(FatalError);
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

heatTransferInterGrain::~heatTransferInterGrain()
{
}

// * * * * * * * * * * * * * * * private Member Functions  * * * * * * * * * * * * * //
// * * * * * * * * * * * * * * * * Member Fct  * * * * * * * * * * * * * * * //

void heatTransferInterGrain::calcEnergyContribution()
{
    calcPartEffThermCond();

    if (kMax_ > 0.0)
    {
        dimensionedScalar kMax("kMax",dimensionSet(1,1,-3,-1,0,0,0),kMax_);
        partEffThermCondField_ = min(partEffThermCondField_,kMax);
    }

    if (!implicit_)
    {
        QPartPart_ = fvc::laplacian(partEffThermCondField_,partTempField_);
    }
    else
    {
        double**& partThermCap_ = particleCloud_.getParticlePropertyRef<double**>(partThermCapRegName_);
        label cellI = -1;
        label partType = -1;
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
                partThermCap_[index][0] = typePartThermCap_[partType - 1] * particleCloud_.particleDensity(index);
            }
        }

        partThermCapField_.primitiveFieldRef() = 0.0;
        particleCloud_.averagingM().resetWeightFields();
        particleCloud_.averagingM().setScalarAverage
        (
            partThermCapField_,
            partThermCap_,
            particleCloud_.particleWeights(),
            particleCloud_.averagingM().UsWeightField(),
            NULL
        );

        partThermCapField_ *= (1.0 - voidfraction_);
        dimensionedScalar CMin(dimensionedScalar("CMin",dimensionSet(1,-1,-2,-1,0,0,0),SMALL));
        partThermCapField_ = max(partThermCapField_, CMin);
        partThermCapField_.oldTime() = partThermCapField_;

        volScalarField partTempField(partTempField_);
        partTempField.oldTime() = partTempField;
        fvScalarMatrix TpEqn
        (
            fvm::ddt(partThermCapField_,partTempField)
            ==
            fvm::laplacian(partEffThermCondField_,partTempField)
        );

        TpEqn.solve();
        QPartPart_ = fvc::laplacian(partEffThermCondField_,partTempField);
    }

    label cellI=0;
    scalar partVolume(0);
    scalar QPartPart(0);
    scalar voidfraction(1);

    totalHeatFlux_ = 0.0;
    double**& partHeatFlux_ = particleCloud_.getParticlePropertyRef<double**>(partHeatFluxName_);

    for(int index = 0;index < particleCloud_.numberOfParticles(); ++index)
    {
        cellI = particleCloud_.cellIDs()[index][0];
        if(cellI >= 0)
        {
            voidfraction = voidfraction_[cellI];

            partVolume = particleCloud_.particleVolume(index);
            QPartPart = QPartPart_[cellI];

            heatFlux(index, partVolume, voidfraction, QPartPart);
            if (calcTotalHeatFlux_) totalHeatFlux_ += partHeatFlux_[index][0];
        }
    }
}

void heatTransferInterGrain::calcPartEffThermCond()
{
    calcPartThermCond();
    if (radiativeHeatTransfer_) calcPartThermRad();

    scalar volFrac = 0.0;

    if (!radiativeHeatTransfer_)
    {
        forAll(partEffThermCondField_, cellI)
        {
            volFrac = 1.0 - prescribedVoidfraction_[cellI];
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
    else
    {
        scalar thRad = 0.0;
        forAll(partEffThermCondField_, cellI)
        {
            volFrac = 1.0 - prescribedVoidfraction_[cellI];
            thRad = volFrac * partThermRadField_[cellI]; // volume fraction factor for rad. contribution according to Qian et al.
            if (volFrac < 0.334)
            {
                partEffThermCondField_[cellI] = thRad;
            }
            else
            {
                partEffThermCondField_[cellI] = 0.5 * (3*volFrac - 1) * partThermCondField_[cellI] + thRad;
            }
        }
    }
}

void heatTransferInterGrain::calcPartThermCond()
{
    label cellI=0;
    label partType = 1;
    double**& partThermCond_ = particleCloud_.getParticlePropertyRef<double**>(partThermCondRegName_);

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

void heatTransferInterGrain::calcPartThermRad()
{
    label cellI=0;
    label partType = 1;
    scalar cg = typeCG_[0];
    scalar ds = 0.0;
    scalar L = 0.0;
    scalar prefac = 0.0;
    scalar Tp = 0.0;
    scalar voidfraction = 0.0;

    double**& partTemp_ = particleCloud_.getParticlePropertyRef<double**>(partTempName_);
    double**& partThermCond_ = particleCloud_.getParticlePropertyRef<double**>(partThermCondRegName_);
    double**& partThermRad_ = particleCloud_.getParticlePropertyRef<double**>(partThermRadRegName_);

    for(int index = 0;index < particleCloud_.numberOfParticles(); ++index)
    {
            cellI = particleCloud_.cellIDs()[index][0];
            if (cellI >= 0)
            {
                if (multiTypes_)
                {
                    partType = particleCloud_.particleType(index);
                    if (partType > maxTypeCG_)
                    {
                        FatalError<< "Too few coarse-graining factors provided." << abort(FatalError);
                    }
                    cg = typeCG_[partType - 1];
                }
                ds = 2.*particleCloud_.radius(index)/cg;
                // make sure reasonable values are used
                Tp = partTemp_[index][0];
                if (Tp < TMin) Tp = TMin;
                voidfraction = prescribedVoidfraction_[cellI];
                if (voidfraction < voidfracMin) voidfraction = voidfracMin;
                else if (voidfraction > voidfracMax) voidfraction = voidfracMax;

                prefac = 4.0*constant::physicoChemical::sigma.value()*ds*Tp*Tp*Tp;
                L = partThermCond_[index][0]/prefac;
                // LIGGGGHTS counts types 1, 2, ..., C++ array starts at 0
                partThermRad_[index][0] = prefac*FE(voidfraction,typePartEmissivity_[partType - 1],L);
            }
    }

    partThermRadField_.primitiveFieldRef() = 0.0;
    particleCloud_.averagingM().resetWeightFields();
    particleCloud_.averagingM().setScalarAverage
    (
        partThermRadField_,
        partThermRad_,
        particleCloud_.particleWeights(),
        particleCloud_.averagingM().UsWeightField(),
        NULL
    );
}

scalar heatTransferInterGrain::FE(scalar voidfraction, scalar emissivity, scalar L)
{
    scalar B = 1.25 * pow((1-voidfraction)/voidfraction,1.111);
    scalar sqrt_aP = sqrt(1-voidfraction);
    scalar x = 2.0/emissivity - 1.0;
    scalar y = 1.0 / (x*L);
    scalar z = sqrt_aP / (2.0/emissivity - 1.0) * (B+1.0)/B * 1.0/(1.0 + y);
    return (1 - sqrt_aP)*voidfraction + z;
}

void heatTransferInterGrain::heatFlux(label index, scalar vol, scalar voidfraction, scalar QPartPart)
{
    double**& partHeatFlux_ = particleCloud_.getParticlePropertyRef<double**>(partHeatFluxName_);
    partHeatFlux_[index][0] = vol * QPartPart / (1.0 - voidfraction) ;
}

void heatTransferInterGrain::giveData()
{
    if (calcTotalHeatFlux_)
    {
        reduce(totalHeatFlux_, sumOp<scalar>());
        Info << "total conductive particle-particle heat flux [W] (Eulerian) = " << totalHeatFlux_ << endl;
    }

    double**& partHeatFlux_ = particleCloud_.getParticlePropertyRef<double**>(partHeatFluxName_);
    particleCloud_.dataExchangeM().giveData(partHeatFluxName_,"scalar-atom", partHeatFlux_);
}

void heatTransferInterGrain::postFlow()
{
    giveData();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

