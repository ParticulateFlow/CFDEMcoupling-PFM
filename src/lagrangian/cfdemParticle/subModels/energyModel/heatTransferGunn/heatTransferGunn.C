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
#include "heatTransferGunn.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(heatTransferGunn, 0);

addToRunTimeSelectionTable(energyModel, heatTransferGunn, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
heatTransferGunn::heatTransferGunn
(
    const dictionary& dict,
    cfdemCloudEnergy& sm
)
:
    energyModel(dict,sm),
    propsDict_(dict.subDict(typeName + "Props")),
    multiTypes_(false),
    expNusselt_(propsDict_.lookupOrDefault<bool>("expNusselt",false)),
    interpolation_(propsDict_.lookupOrDefault<bool>("interpolation",false)),
    verbose_(propsDict_.lookupOrDefault<bool>("verbose",false)),
    implicit_(propsDict_.lookupOrDefault<bool>("implicit",true)),
    QPartFluidName_(propsDict_.lookupOrDefault<word>("QPartFluidName","QPartFluid")),
    QPartFluid_
    (   IOobject
        (
            QPartFluidName_,
            sm.mesh().time().timeName(),
            sm.mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        sm.mesh(),
        dimensionedScalar("zero", dimensionSet(1,-1,-3,0,0,0,0), 0.0)
    ),
    QPartFluidCoeffName_(propsDict_.lookupOrDefault<word>("QPartFluidCoeffName","QPartFluidCoeff")),
    QPartFluidCoeff_
    (   IOobject
        (
            QPartFluidCoeffName_,
            sm.mesh().time().timeName(),
            sm.mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        sm.mesh(),
        dimensionedScalar("zero", dimensionSet(1,-1,-3,-1,0,0,0), 0.0)
    ),
    partTempField_
    (   IOobject
        (
            "partTemp",
            sm.mesh().time().timeName(),
            sm.mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        sm.mesh(),
        dimensionedScalar("zero", dimensionSet(0,0,0,1,0,0,0), 0.0),
        "zeroGradient"
    ),
    partRelTempField_
    (   IOobject
        (
            "particleRelTemp",
            sm.mesh().time().timeName(),
            sm.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        sm.mesh(),
        dimensionedScalar("zero", dimensionSet(0,0,0,0,0,0,0), 0.0)
    ),
    ReField_
    (   IOobject
        (
            "ReField",
            sm.mesh().time().timeName(),
            sm.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        sm.mesh(),
        dimensionedScalar("zero", dimensionSet(0,0,0,0,0,0,0), 0.0)
    ),
    NuField_
    (   IOobject
        (
            "NuField",
            sm.mesh().time().timeName(),
            sm.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        sm.mesh(),
        dimensionedScalar("zero", dimensionSet(0,0,0,0,0,0,0), 0.0)
    ),
    partRefTemp_("partRefTemp", dimensionSet(0,0,0,1,0,0,0), 0.0),
    calcPartTempField_(propsDict_.lookupOrDefault<bool>("calcPartTempField",false)),
    calcPartTempAve_(propsDict_.lookupOrDefault<bool>("calcPartTempAve",false)),
    partTempAve_(0.0),
    tempFieldName_(propsDict_.lookupOrDefault<word>("tempFieldName","T")),
    tempField_(sm.mesh().lookupObject<volScalarField> (tempFieldName_)),
    voidfractionFieldName_(propsDict_.lookupOrDefault<word>("voidfractionFieldName","voidfraction")),
    voidfraction_(sm.mesh().lookupObject<volScalarField> (voidfractionFieldName_)),
    maxSource_(1e30),
    velFieldName_(propsDict_.lookupOrDefault<word>("velFieldName","U")),
    U_(sm.mesh().lookupObject<volVectorField> (velFieldName_)),
    densityFieldName_(propsDict_.lookupOrDefault<word>("densityFieldName","rho")),
    rho_(sm.mesh().lookupObject<volScalarField> (densityFieldName_)),
    partTempName_(propsDict_.lookup("partTempName")),
    partTemp_(NULL),
    partHeatFluxName_(propsDict_.lookup("partHeatFluxName")),
    partHeatFlux_(NULL),
    partHeatFluxCoeff_(NULL),
    partRe_(NULL),
    partNu_(NULL),
    scaleDia_(1.),
    typeCG_(propsDict_.lookupOrDefault<scalarList>("coarseGrainingFactors",scalarList(1,1.0)))
{
    allocateMyArrays();

    if (propsDict_.found("maxSource"))
    {
        maxSource_=readScalar(propsDict_.lookup ("maxSource"));
        Info << "limiting eulerian source field to: " << maxSource_ << endl;
    }

    if (calcPartTempField_)
    {
        calcPartTempAve_ = true;
        if (propsDict_.found("partRefTemp"))
        {
            partRefTemp_.value()=readScalar(propsDict_.lookup ("partRefTemp"));
        }
        partTempField_.writeOpt() = IOobject::AUTO_WRITE;
        partRelTempField_.writeOpt() = IOobject::AUTO_WRITE;
        partTempField_.write();
        partRelTempField_.write();
        Info <<  "Particle temperature field activated." << endl;
    }

    if (!implicit_)
    {
        QPartFluidCoeff_.writeOpt() = IOobject::NO_WRITE;
    }

    if (verbose_)
    {
        ReField_.writeOpt() = IOobject::AUTO_WRITE;
        NuField_.writeOpt() = IOobject::AUTO_WRITE;
        ReField_.write();
        NuField_.write();
        if (expNusselt_)
        {
          FatalError <<"Cannot read and create NuField at the same time!\n" << abort(FatalError);
        }
    }

    if (propsDict_.found("scale") && typeCG_.size()==1)
    {
        // if "scale" is specified and there's only one single type, use "scale"
        scaleDia_=scalar(readScalar(propsDict_.lookup("scale")));
        typeCG_[0] = scaleDia_;
    }
    else if (typeCG_.size()>1) multiTypes_ = true;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

heatTransferGunn::~heatTransferGunn()
{
    particleCloud_.dataExchangeM().destroy(partTemp_,1);
    particleCloud_.dataExchangeM().destroy(partHeatFlux_,1);
    particleCloud_.dataExchangeM().destroy(partRe_,1);
    particleCloud_.dataExchangeM().destroy(partNu_,1);
    if (implicit_)
    {
        particleCloud_.dataExchangeM().destroy(partHeatFluxCoeff_,1);
    }
}

// * * * * * * * * * * * * * * * private Member Functions  * * * * * * * * * * * * * //
void heatTransferGunn::allocateMyArrays() const
{
    // get memory for 2d arrays
    double initVal=0.0;
    particleCloud_.dataExchangeM().allocateArray(partTemp_,initVal,1);  // field/initVal/with/lenghtFromLigghts
    particleCloud_.dataExchangeM().allocateArray(partHeatFlux_,initVal,1);
    if(implicit_)
    {
        particleCloud_.dataExchangeM().allocateArray(partHeatFluxCoeff_,initVal,1);
    }

    if(verbose_)
    {
        particleCloud_.dataExchangeM().allocateArray(partRe_,initVal,1);
        particleCloud_.dataExchangeM().allocateArray(partNu_,initVal,1);
    }
}

// * * * * * * * * * * * * * * * * Member Fct  * * * * * * * * * * * * * * * //

void heatTransferGunn::calcEnergyContribution()
{
   // realloc the arrays
    allocateMyArrays();

    // reset Scalar field
    QPartFluid_.primitiveFieldRef() = 0.0;

    // get DEM data
    particleCloud_.dataExchangeM().getData(partTempName_,"scalar-atom",partTemp_);

    if(calcPartTempField_)
    {
        partTempField_.primitiveFieldRef() = 0.0;
        particleCloud_.averagingM().resetWeightFields();
        particleCloud_.averagingM().setScalarAverage
        (
            partTempField_,
            partTemp_,
            particleCloud_.particleWeights(),
            particleCloud_.averagingM().UsWeightField(),
            NULL
        );

        volScalarField sumTp (particleCloud_.averagingM().UsWeightField() * partTempField_);
        dimensionedScalar aveTemp("aveTemp",dimensionSet(0,0,0,1,0,0,0), gSum(sumTp) / particleCloud_.numberOfParticles());
        partRelTempField_ = (partTempField_ - aveTemp) / (aveTemp - partRefTemp_);
        Info << "heatTransferGunn: average part. temp = " << aveTemp.value() << endl;
    }

    #ifdef compre
       const volScalarField mufField = particleCloud_.turbulence().mu();
    #else
       const volScalarField mufField = particleCloud_.turbulence().nu()*rho_;
    #endif

    if (typeCG_.size()>1 || typeCG_[0] > 1)
    {
        Info << "heatTransferGunn using scale = " << typeCG_ << endl;
    }
    else if (particleCloud_.cg() > 1)
    {
        scaleDia_=particleCloud_.cg();
        Info << "heatTransferGunn using scale from liggghts cg = " << scaleDia_ << endl;
    }

    // calc La based heat flux
    scalar voidfraction(1);
    vector Ufluid(0,0,0);
    scalar Tfluid(0);
    label cellI=0;
    vector Us(0,0,0);
    scalar ds(0);
    scalar ds_scaled(0);
    scalar scaleDia3 = typeCG_[0]*typeCG_[0]*typeCG_[0];
    scalar muf(0);
    scalar magUr(0);
    scalar Rep(0);
    scalar Pr(0);
    scalar Nup(0);
    scalar Tsum(0.0);

    scalar cg = typeCG_[0];
    label partType = 1;

    interpolationCellPoint<scalar> voidfractionInterpolator_(voidfraction_);
    interpolationCellPoint<vector> UInterpolator_(U_);
    interpolationCellPoint<scalar> TInterpolator_(tempField_);

    for(int index = 0;index < particleCloud_.numberOfParticles(); ++index)
    {
            cellI = particleCloud_.cellIDs()[index][0];
            if(cellI >= 0)
            {
                if(interpolation_)
                {
                    vector position = particleCloud_.position(index);
                    voidfraction = voidfractionInterpolator_.interpolate(position,cellI);
                    Ufluid = UInterpolator_.interpolate(position,cellI);
                    Tfluid = TInterpolator_.interpolate(position,cellI);
                }
                else
                {
                    voidfraction = voidfraction_[cellI];
                    Ufluid = U_[cellI];
                    Tfluid = tempField_[cellI];
                }

                if (voidfraction < 0.01)
                    voidfraction = 0.01;

                if (multiTypes_)
                {
                    partType = particleCloud_.particleType(index);
                    cg = typeCG_[partType - 1];
                    scaleDia3 = cg*cg*cg;
                }

                // calc relative velocity
                Us = particleCloud_.velocity(index);
                magUr = mag(Ufluid - Us);
                ds = 2.*particleCloud_.radius(index);
                ds_scaled = ds/cg;
                muf = mufField[cellI];
                Rep = ds_scaled * magUr * voidfraction * rho_[cellI]/ muf;
                Pr = max(SMALL, Cp_ * muf / kf0_);

                Nup = Nusselt(voidfraction, Rep, Pr);

                Tsum += partTemp_[index][0];
                scalar h = kf0_ * Nup / ds_scaled;
                scalar As = ds_scaled * ds_scaled * M_PI; // surface area of sphere

                // calc convective heat flux [W]
                heatFlux(index, h, As, Tfluid, scaleDia3);

                if(verbose_)
                {
                    partRe_[index][0] = Rep;
                    partNu_[index][0] = Nup;
                }

                if(particleCloud_.verbose() && index >=0 && index <2)
                {
                    Info << "partHeatFlux = " << partHeatFlux_[index][0] << endl;
                    Info << "magUr = " << magUr << endl;
                    Info << "As = " << As << endl;
                    Info << "muf = " << muf << endl;
                    Info << "Rep = " << Rep << endl;
                    Info << "Pr = " << Pr << endl;
                    Info << "Nup = " << Nup << endl;
                    Info << "voidfraction = " << voidfraction << endl;
                    Info << "partTemp_[index][0] = " << partTemp_[index][0] << endl;
                    Info << "Tfluid = " << Tfluid << endl  ;
                }
            }
    }

    // gather particle temperature sums and obtain average
    if(calcPartTempAve_)
    {
        reduce(Tsum, sumOp<scalar>());
        partTempAve_ = Tsum / particleCloud_.numberOfParticles();
        Info << "mean particle temperature = " << partTempAve_ << endl;
    }

    if(calcPartTempField_) partTempField();

    particleCloud_.averagingM().setScalarSum
    (
        QPartFluid_,
        partHeatFlux_,
        particleCloud_.particleWeights(),
        NULL
    );

    QPartFluid_.primitiveFieldRef() /= -QPartFluid_.mesh().V();

    if(implicit_)
    {
        QPartFluidCoeff_.primitiveFieldRef() = 0.0;

        particleCloud_.averagingM().setScalarSum
        (
            QPartFluidCoeff_,
            partHeatFluxCoeff_,
            particleCloud_.particleWeights(),
            NULL
        );

        QPartFluidCoeff_.primitiveFieldRef() /= -QPartFluidCoeff_.mesh().V();
    }

    if(verbose_)
    {
        ReField_.primitiveFieldRef() = 0.0;
        NuField_.primitiveFieldRef() = 0.0;
        particleCloud_.averagingM().resetWeightFields();
        particleCloud_.averagingM().setScalarAverage
        (
            ReField_,
            partRe_,
            particleCloud_.particleWeights(),
            particleCloud_.averagingM().UsWeightField(),
            NULL
        );
        particleCloud_.averagingM().resetWeightFields();
        particleCloud_.averagingM().setScalarAverage
        (
            NuField_,
            partNu_,
            particleCloud_.particleWeights(),
            particleCloud_.averagingM().UsWeightField(),
            NULL
        );
    }

    // limit source term in explicit treatment
    if(!implicit_)
    {
        forAll(QPartFluid_,cellI)
        {
            scalar EuFieldInCell = QPartFluid_[cellI];

            if(mag(EuFieldInCell) > maxSource_ )
            {
                 Pout << "limiting source term\n"  << endl  ;
                 QPartFluid_[cellI] = sign(EuFieldInCell) * maxSource_;
            }
        }
    }

    QPartFluid_.correctBoundaryConditions();
}

void heatTransferGunn::addEnergyContribution(volScalarField& Qsource) const
{
    Qsource += QPartFluid_;
}

void heatTransferGunn::addEnergyCoefficient(volScalarField& Qsource) const
{
    if(implicit_)
    {
        Qsource += QPartFluidCoeff_;
    }
}

scalar heatTransferGunn::Nusselt(scalar voidfraction, scalar Rep, scalar Pr) const
{
    return (7 - 10 * voidfraction + 5 * voidfraction * voidfraction) *
                        (1 + 0.7 * Foam::pow(Rep,0.2) * Foam::pow(Pr,0.33)) +
                        (1.33 - 2.4 * voidfraction + 1.2 * voidfraction * voidfraction) *
                        Foam::pow(Rep,0.7) * Foam::pow(Pr,0.33);
}

void heatTransferGunn::heatFlux(label index, scalar h, scalar As, scalar Tfluid, scalar cg3)
{
    scalar hAs = h * As * cg3;

    if (particleCloud_.getParticleEffVolFactors())
    {
        scalar effVolFac = particleCloud_.particleEffVolFactor(index);
        hAs *= effVolFac;
    }

    partHeatFlux_[index][0] = - hAs * partTemp_[index][0];
    if(!implicit_)
    {
        partHeatFlux_[index][0] += hAs * Tfluid;
    }
    else
    {
        partHeatFluxCoeff_[index][0] = hAs;
    }
}

void heatTransferGunn::giveData()
{
    Info << "total convective particle-fluid heat flux [W] (Eulerian) = " << gSum(QPartFluid_*1.0*QPartFluid_.mesh().V()) << endl;

    particleCloud_.dataExchangeM().giveData(partHeatFluxName_,"scalar-atom", partHeatFlux_);
}

void heatTransferGunn::postFlow()
{
    if(implicit_)
    {
        label cellI;
        scalar Tfluid(0.0);
        scalar Tpart(0.0);
        interpolationCellPoint<scalar> TInterpolator_(tempField_);

        for(int index = 0;index < particleCloud_.numberOfParticles(); ++index)
        {
                cellI = particleCloud_.cellIDs()[index][0];
                if(cellI >= 0)
                {
                    if(interpolation_)
                    {
                        vector position = particleCloud_.position(index);
                        Tfluid = TInterpolator_.interpolate(position,cellI);
                    }
                    else
                    {
                        Tfluid = tempField_[cellI];
                    }

                    Tpart = partTemp_[index][0];
                    partHeatFlux_[index][0] = (Tfluid - Tpart) * partHeatFluxCoeff_[index][0];
                }
        }
    }
    giveData();
}

scalar heatTransferGunn::aveTpart() const
{
    return partTempAve_;
}

void heatTransferGunn::partTempField()
{
    partTempField_.primitiveFieldRef() = 0.0;
    particleCloud_.averagingM().resetWeightFields();
    particleCloud_.averagingM().setScalarAverage
    (
        partTempField_,
        partTemp_,
        particleCloud_.particleWeights(),
        particleCloud_.averagingM().UsWeightField(),
        NULL
    );

    dimensionedScalar aveTemp("aveTemp",dimensionSet(0,0,0,1,0,0,0), partTempAve_);
    partRelTempField_ = (partTempField_ - aveTemp) / (aveTemp - partRefTemp_);
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

