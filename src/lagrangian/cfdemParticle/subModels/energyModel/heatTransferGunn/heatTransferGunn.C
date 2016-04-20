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
    interpolation_(propsDict_.lookupOrDefault<bool>("interpolation",false)),
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
        dimensionedScalar("zero", dimensionSet(1,1,-3,0,0,0,0), 0.0)
    ),
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
    partHeatFlux_(NULL)
{
     allocateMyArrays();

    if (propsDict_.found("maxSource"))
    {
        maxSource_=readScalar(propsDict_.lookup ("maxSource"));
        Info << "limiting eulerian source field to: " << maxSource_ << endl;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

heatTransferGunn::~heatTransferGunn()
{
    delete partTemp_;
    delete partHeatFlux_;
}

// * * * * * * * * * * * * * * * private Member Functions  * * * * * * * * * * * * * //
void heatTransferGunn::allocateMyArrays() const
{
    // get memory for 2d arrays
    double initVal=0.0;
    particleCloud_.dataExchangeM().allocateArray(partTemp_,initVal,1);  // field/initVal/with/lenghtFromLigghts
    particleCloud_.dataExchangeM().allocateArray(partHeatFlux_,initVal,1);
}

// * * * * * * * * * * * * * * * * Member Fct  * * * * * * * * * * * * * * * //

void heatTransferGunn::calcEnergyContribution()
{
   // realloc the arrays
    allocateMyArrays();

    // reset Scalar field
    QPartFluid_.internalField() = 0.0;

    // get DEM data
    particleCloud_.dataExchangeM().getData(partTempName_,"scalar-atom",partTemp_);

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

                // calc relative velocity
                Us = particleCloud_.velocity(index);
                magUr = mag(Ufluid - Us);
                ds = 2.*particleCloud_.radius(index);
                muf = mufField[cellI];
                Rep = ds * magUr * voidfraction * rho_[cellI]/ muf;
                Pr = max(SMALL, Cp_ * muf / kf0_);

                Nup = (7 - 10 * voidfraction + 5 * voidfraction * voidfraction) *
                        (1 + 0.7 * Foam::pow(Rep,0.2) * Foam::pow(Pr,0.33)) +
			(1.33 - 2.4 * voidfraction + 1.2 * voidfraction * voidfraction) *
			Foam::pow(Rep,0.7) * Foam::pow(Pr,0.33);                      
                

                scalar h = kf0_ * Nup / ds;
                scalar As = ds * ds * M_PI; // surface area of sphere

                // calc convective heat flux [W]
                scalar partHeatFlux = h * As * (Tfluid - partTemp_[index][0]);
                partHeatFlux_[index][0] = partHeatFlux;
		
		if(particleCloud_.verbose() && index >=0 && index <2)
                {
                    Info << "partHeatFlux = " << partHeatFlux << endl;
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

    QPartFluid_.internalField() /= -QPartFluid_.mesh().V();

    // limit source term
    forAll(QPartFluid_,cellI)
    {
        scalar EuFieldInCell = QPartFluid_[cellI];

        if(mag(EuFieldInCell) > maxSource_ )
        {
             QPartFluid_[cellI] = sign(EuFieldInCell) * maxSource_;
        }
    }
    
    QPartFluid_.correctBoundaryConditions();

    Info << "total convective particle-fluid heat flux [W] (Eulerian) = " << gSum(QPartFluid_*1.0*QPartFluid_.mesh().V()) << endl;

    // give DEM data
    particleCloud_.dataExchangeM().giveData(partHeatFluxName_,"scalar-atom", partHeatFlux_); 
}

void heatTransferGunn::addEnergyContribution(volScalarField& Qsource) const
{
    Qsource += QPartFluid_;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
