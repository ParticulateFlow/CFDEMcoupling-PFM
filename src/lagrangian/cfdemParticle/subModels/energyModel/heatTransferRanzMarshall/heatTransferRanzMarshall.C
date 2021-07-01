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
#include "heatTransferRanzMarshall.H"
#include "addToRunTimeSelectionTable.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(heatTransferRanzMarshall, 0);

addToRunTimeSelectionTable(energyModel, heatTransferRanzMarshall, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
heatTransferRanzMarshall::heatTransferRanzMarshall
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
    calcTotalHeatFlux_(propsDict_.lookupOrDefault<bool>("calcTotalHeatFlux",true)),
    initPartTemp_(propsDict_.lookupOrDefault<bool>("initPartTemp",false)),
    Tmin_(propsDict_.lookupOrDefault<scalar>("Tmin",0.0)),
    Tmax_(propsDict_.lookupOrDefault<scalar>("Tmax",1e6)),
    totalHeatFlux_(0.0),
    NusseltScalingFactor_(1.0),
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
    partTempAve_("partTempAve", dimensionSet(0,0,0,1,0,0,0), 0.0),
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
    partHeatFluxName_(propsDict_.lookupOrDefault<word>("partHeatFluxName","convectiveHeatFlux")),
    partHeatFluxCoeffRegName_(typeName + "partHeatFluxCoeff"),
    partReRegName_(typeName + "partRe"),
    partNuRegName_(typeName + "partNu"),
    scaleDia_(1.),
    typeCG_(propsDict_.lookupOrDefault<scalarList>("coarseGrainingFactors",scalarList(1,1.0))),
    maxTypeCG_(typeCG_.size())
{
    particleCloud_.registerParticleProperty<double**>(partTempName_,1);
    particleCloud_.registerParticleProperty<double**>(partHeatFluxName_,1);
    if (implicit_)
    {
        particleCloud_.registerParticleProperty<double**>(partHeatFluxCoeffRegName_,1);
    }
    if(verbose_)
    {
        particleCloud_.registerParticleProperty<double**>(partReRegName_,1);
        particleCloud_.registerParticleProperty<double**>(partNuRegName_,1);
    }

    if (propsDict_.found("NusseltScalingFactor"))
    {
        NusseltScalingFactor_=readScalar(propsDict_.lookup ("NusseltScalingFactor"));
        Info << "NusseltScalingFactor set to: " << NusseltScalingFactor_ << endl;
    }

    if (propsDict_.found("maxSource"))
    {
        maxSource_=readScalar(propsDict_.lookup ("maxSource"));
        Info << "limiting eulerian source field to: " << maxSource_ << endl;
    }

    if (calcPartTempField_)
    {
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

    if (expNusselt_)
    {
        NuField_.writeOpt() = IOobject::AUTO_WRITE;
        NuField_.write();
        Info <<  "Using predefined Nusselt number field." << endl;
    }

    if (propsDict_.found("scale") && typeCG_.size()==1)
    {
        // if "scale" is specified and there's only one single type, use "scale"
        scaleDia_=scalar(readScalar(propsDict_.lookup("scale")));
        typeCG_[0] = scaleDia_;
    }
    else if (typeCG_.size()>1)
    {
        multiTypes_ = true;
    }

    if (initPartTemp_ && !partTempField_.headerOk())
    {
        FatalError <<"Trying to initialize particle temperatures, but no field found.\n" << abort(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

heatTransferRanzMarshall::~heatTransferRanzMarshall()
{
}

// * * * * * * * * * * * * * * * private Member Functions  * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Member Fct  * * * * * * * * * * * * * * * //

void heatTransferRanzMarshall::calcEnergyContribution()
{
    double**& partTemp_ = particleCloud_.getParticlePropertyRef<double**>(partTempName_);
    double**& partHeatFlux_ = particleCloud_.getParticlePropertyRef<double**>(partHeatFluxName_);

    // reset Scalar field
    QPartFluid_.primitiveFieldRef() = 0.0;

    if (initPartTemp_)
    {
        // if particle temperatures are to be initialized from field, do a one-time push to DEM
        initPartTemp();
        particleCloud_.dataExchangeM().giveData("Temp","scalar-atom", partTemp_);
        initPartTemp_ = false;
    }
    else
    {
        // get DEM data
        particleCloud_.dataExchangeM().getData(partTempName_,"scalar-atom",partTemp_);
    }

    #ifdef compre
       const volScalarField mufField = particleCloud_.turbulence().mu();
    #else
       const volScalarField mufField = particleCloud_.turbulence().nu()*rho_;
    #endif

    if (typeCG_.size()>1 || typeCG_[0] > 1)
    {
        Info << "heatTransferRanzMarshall using scale = " << typeCG_ << endl;
    }
    else if (particleCloud_.cg() > 1)
    {
        scaleDia_=particleCloud_.cg();
        Info << "heatTransferRanzMarshall using scale from liggghts cg = " << scaleDia_ << endl;
    }

    // calc La based heat flux
    scalar voidfraction(1);
    vector Ufluid(0,0,0);
    scalar Tfluid(0);
    label cellI=0;
    vector Us(0,0,0);
    scalar ds(0);
    scalar ds_scaled(0);
    scalar scaleDia3 = scaleDia_*scaleDia_*scaleDia_;
    scalar muf(0);
    scalar magUr(0);
    scalar Rep(0);
    scalar Pr(0);
    scalar Nup(0);
    scalar Tsum(0.0);
    scalar Nsum(0.0);

    scalar cg = scaleDia_;
    label partType = 1;

    interpolationCellPoint<scalar> voidfractionInterpolator_(voidfraction_);
    interpolationCellPoint<vector> UInterpolator_(U_);
    interpolationCellPoint<scalar> TInterpolator_(tempField_);

    totalHeatFlux_ = 0.0;

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
                    if (partType > maxTypeCG_)
                    {
                        FatalError<< "Too few coarse-graining factors provided." << abort(FatalError);
                    }
                    cg = typeCG_[partType - 1];
                    scaleDia3 = cg*cg*cg;
                }

                ds = 2.*particleCloud_.radius(index);
                ds_scaled = ds/cg;

                if (expNusselt_)
                {
                    Nup = NuField_[cellI];
                    if (Nup < 2.0)
                        Nup = 2.0;
                }
                else
                {
                    Us = particleCloud_.velocity(index);
                    magUr = mag(Ufluid - Us);
                    muf = mufField[cellI];
                    Rep = ds_scaled * magUr * voidfraction * rho_[cellI]/ muf;
                    Pr = max(SMALL, Cp_ * muf / kf0_);
                    Nup = Nusselt(voidfraction, Rep, Pr);
                }
                Nup *= NusseltScalingFactor_;

                Tsum += partTemp_[index][0];
                Nsum += 1.0;

                scalar h = kf0_ * Nup / ds_scaled;
                scalar As = ds_scaled * ds_scaled * M_PI; // surface area of sphere

                // calc convective heat flux [W]
                heatFlux(index, h, As, Tfluid, scaleDia3);
                if (calcTotalHeatFlux_ & !implicit_)
                {
                    totalHeatFlux_ += partHeatFlux_[index][0];
                }

                if(verbose_)
                {
                    double**& partRe_ = particleCloud_.getParticlePropertyRef<double**>(partReRegName_);
                    double**& partNu_ = particleCloud_.getParticlePropertyRef<double**>(partNuRegName_);
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

    if(calcPartTempField_)
    {
        reduce(Tsum, sumOp<scalar>());
        reduce(Nsum, sumOp<scalar>());
        partTempAve_.value() = Tsum / Nsum;
        partTempField();
    }

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
        double**& partHeatFluxCoeff_ = particleCloud_.getParticlePropertyRef<double**>(partHeatFluxCoeffRegName_);
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
        double**& partRe_ = particleCloud_.getParticlePropertyRef<double**>(partReRegName_);
        double**& partNu_ = particleCloud_.getParticlePropertyRef<double**>(partNuRegName_);
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

void heatTransferRanzMarshall::addEnergyContribution(volScalarField& Qsource) const
{
    Qsource += QPartFluid_;
}

void heatTransferRanzMarshall::addEnergyCoefficient(volScalarField& Qsource) const
{
    if(implicit_)
    {
        Qsource += QPartFluidCoeff_;
    }
}

scalar heatTransferRanzMarshall::Nusselt(scalar voidfraction, scalar Rep, scalar Pr) const
{
    return (2 + 0.6 * Foam::pow(Rep,0.5) * Foam::pow(Pr,0.33));
}

void heatTransferRanzMarshall::heatFlux(label index, scalar h, scalar As, scalar Tfluid, scalar cg3)
{
    scalar hAs = h * As * cg3;
    double**& partTemp_ = particleCloud_.getParticlePropertyRef<double**>(partTempName_);
    double**& partHeatFlux_ = particleCloud_.getParticlePropertyRef<double**>(partHeatFluxName_);

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
        double**& partHeatFluxCoeff_ = particleCloud_.getParticlePropertyRef<double**>(partHeatFluxCoeffRegName_);
        partHeatFluxCoeff_[index][0] = hAs;
    }
}

void heatTransferRanzMarshall::giveData()
{
   // Info << "total convective particle-fluid heat flux [W] (Eulerian) = " << gSum(QPartFluid_*1.0*QPartFluid_.mesh().V()) << endl;
    if (calcTotalHeatFlux_)
    {
        reduce(totalHeatFlux_, sumOp<scalar>());
        Info << "total convective particle-fluid heat flux [W] = " << totalHeatFlux_ << endl;
    }
    double**& partHeatFlux_ = particleCloud_.getParticlePropertyRef<double**>(partHeatFluxName_);
    particleCloud_.dataExchangeM().giveData(partHeatFluxName_,"scalar-atom", partHeatFlux_);
}

void heatTransferRanzMarshall::postFlow()
{
    if(implicit_)
    {
        label cellI;
        scalar Tfluid(0.0);
        scalar Tpart(0.0);
        interpolationCellPoint<scalar> TInterpolator_(tempField_);
        double**& partTemp_ = particleCloud_.getParticlePropertyRef<double**>(partTempName_);
        double**& partHeatFlux_ = particleCloud_.getParticlePropertyRef<double**>(partHeatFluxName_);
        double**& partHeatFluxCoeff_ = particleCloud_.getParticlePropertyRef<double**>(partHeatFluxCoeffRegName_);

        totalHeatFlux_ = 0.0;

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
                    if (calcTotalHeatFlux_)
                    {
                        totalHeatFlux_ += partHeatFlux_[index][0];
                    }
                }
        }
    }
    giveData();
}


void heatTransferRanzMarshall::partTempField()
{
        double**& partTemp_ = particleCloud_.getParticlePropertyRef<double**>(partTempName_);
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

        dimensionedScalar denom = partTempAve_ - partRefTemp_;
        if (denom.value() < SMALL && denom.value() > -SMALL) denom.value() = SMALL;
        partRelTempField_ = (partTempField_ - partTempAve_) / denom;

        Info << "heatTransferRanzMarshall: average part. temp = " << partTempAve_.value() << endl;
}

void heatTransferRanzMarshall::initPartTemp()
{
    label cellI = 0;
    scalar T = 0.0;
    double**& partTemp_ = particleCloud_.getParticlePropertyRef<double**>(partTempName_);

    for(int index = 0;index < particleCloud_.numberOfParticles(); ++index)
    {
        cellI = particleCloud_.cellIDs()[index][0];
        if(cellI >= 0)
        {
            T = partTempField_[cellI];
            if (T < Tmin_)
            {
                T = Tmin_;
            }
            else if (T > Tmax_)
            {
                T = Tmax_;
            }
            partTemp_[index][0] = T;
        }
    }
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

