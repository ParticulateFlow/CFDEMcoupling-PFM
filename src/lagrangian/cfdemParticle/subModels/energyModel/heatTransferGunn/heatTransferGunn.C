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
    expNusselt_(propsDict_.lookupOrDefault<bool>("expNusselt",false)),
    interpolation_(propsDict_.lookupOrDefault<bool>("interpolation",false)),
    verbose_(propsDict_.lookupOrDefault<bool>("verbose",false)),
    NusseltFieldName_(propsDict_.lookupOrDefault<word>("NusseltFieldName","NuField")),
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
    partTempField_
    (   IOobject
        (
            "particleTemp",
            sm.mesh().time().timeName(),
            sm.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        sm.mesh(),
        dimensionedScalar("zero", dimensionSet(0,0,0,1,0,0,0), 0.0)
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
            NusseltFieldName_,
            sm.mesh().time().timeName(),
            sm.mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        sm.mesh(),
        dimensionedScalar("zero", dimensionSet(0,0,0,0,0,0,0), 0.0)
    ),
    partRefTemp_("partRefTemp", dimensionSet(0,0,0,1,0,0,0), 0.0),
    calcPartTempField_(propsDict_.lookupOrDefault<bool>("calcPartTempField",false)),
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
    partRe_(NULL),
    partNu_(NULL)
{
     allocateMyArrays();

    if (propsDict_.found("maxSource"))
    {
        maxSource_=readScalar(propsDict_.lookup ("maxSource"));
        Info << "limiting eulerian source field to: " << maxSource_ << endl;
    }
    if (calcPartTempField_)
    {
        if (propsDict_.found("partRefTemp"))
	    partRefTemp_.value()=readScalar(propsDict_.lookup ("partRefTemp"));
	partTempField_.writeOpt() = IOobject::AUTO_WRITE;
	partRelTempField_.writeOpt() = IOobject::AUTO_WRITE;
	partTempField_.write();
	partRelTempField_.write();
	Info <<  "Particle temperature field activated." << endl;
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
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

heatTransferGunn::~heatTransferGunn()
{
    delete partTemp_;
    delete partHeatFlux_;
    delete partRe_;
    delete partNu_;
}

// * * * * * * * * * * * * * * * private Member Functions  * * * * * * * * * * * * * //
void heatTransferGunn::allocateMyArrays() const
{
    // get memory for 2d arrays
    double initVal=0.0;
    particleCloud_.dataExchangeM().allocateArray(partTemp_,initVal,1);  // field/initVal/with/lenghtFromLigghts
    particleCloud_.dataExchangeM().allocateArray(partHeatFlux_,initVal,1);
    
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
    particleCloud_.clockM().start(29,"getDEM_Tdata");
    particleCloud_.dataExchangeM().getData(partTempName_,"scalar-atom",partTemp_);
    particleCloud_.clockM().stop("getDEM_Tdata");
    
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
	volScalarField volP (1 - voidfraction_);
	volScalarField weigthedTp (volP * partTempField_);
	// average per cell-value, not per volume * cell-value
	dimensionedScalar aveTemp = weigthedTp.average() / volP.average();
	partRelTempField_ = (partTempField_ - aveTemp) / (aveTemp - partRefTemp_);
    }

    #ifdef compre
       const volScalarField mufField = particleCloud_.turbulence().mu();
    #else
       const volScalarField mufField = particleCloud_.turbulence().nu()*rho_;
    #endif
    

    // calc La based heat flux
    scalar voidfraction(1);
    vector Ufluid(0,0,0);
    scalar Tfluid(0);
    label cellI=0;
    vector Us(0,0,0);
    scalar ds(0);
    scalar muf(0);
    scalar magUr(0);
    scalar Rep(0);
    scalar Pr(0);
    scalar Nup(0);


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

                ds = 2.*particleCloud_.radius(index);
		
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
                    Rep = ds * magUr * voidfraction * rho_[cellI]/ muf;
                    Pr = max(SMALL, Cp_ * muf / kf0_);

                    Nup = Nusselt(voidfraction, Rep, Pr);
		}
                

                scalar h = kf0_ * Nup / ds;
                scalar As = ds * ds * M_PI; // surface area of sphere

                // calc convective heat flux [W]
                heatFlux(index, h, As, Tfluid);
                heatFluxCoeff(index, h, As);
		
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

    particleCloud_.averagingM().setScalarSum
    (
        QPartFluid_,
        partHeatFlux_,
        particleCloud_.particleWeights(),
        NULL
    );

    QPartFluid_.primitiveFieldRef() /= -QPartFluid_.mesh().V();

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
    
    // limit source term
    forAll(QPartFluid_,cellI)
    {
        scalar EuFieldInCell = QPartFluid_[cellI];

        if(mag(EuFieldInCell) > maxSource_ )
        {
	     Pout << "limiting source term\n"  << endl  ;
             QPartFluid_[cellI] = sign(EuFieldInCell) * maxSource_;
        }
    }
    
    QPartFluid_.correctBoundaryConditions();
    
    giveData(0);

}

void heatTransferGunn::addEnergyContribution(volScalarField& Qsource) const
{
    Qsource += QPartFluid_;
}

scalar heatTransferGunn::Nusselt(scalar voidfraction, scalar Rep, scalar Pr) const
{
    scalar Nup(0.0);
    Nup = (7 - 10 * voidfraction + 5 * voidfraction * voidfraction) *
                        (1 + 0.7 * Foam::pow(Rep,0.2) * Foam::pow(Pr,0.33)) +
			(1.33 - 2.4 * voidfraction + 1.2 * voidfraction * voidfraction) *
			Foam::pow(Rep,0.7) * Foam::pow(Pr,0.33);
                        
    return Nup;
}

void heatTransferGunn::heatFlux(label index, scalar h, scalar As, scalar Tfluid)
{
    partHeatFlux_[index][0] = h * As * (Tfluid - partTemp_[index][0]);
}

void heatTransferGunn::heatFluxCoeff(label index, scalar h, scalar As)
{
   //no heat transfer coefficient in explicit model   
}

void heatTransferGunn::giveData(int call)
{
    if(call == 0)
    {
        Info << "total convective particle-fluid heat flux [W] (Eulerian) = " << gSum(QPartFluid_*1.0*QPartFluid_.mesh().V()) << endl;
	
        particleCloud_.clockM().start(30,"giveDEM_Tdata");
        particleCloud_.dataExchangeM().giveData(partHeatFluxName_,"scalar-atom", partHeatFlux_);
	particleCloud_.clockM().stop("giveDEM_Tdata");
    }
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

