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
#include "YagiWallHTImplicit.H"
#include "addToRunTimeSelectionTable.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(YagiWallHTImplicit, 0);

addToRunTimeSelectionTable(energyModel, YagiWallHTImplicit, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
YagiWallHTImplicit::YagiWallHTImplicit
(
    const dictionary& dict,
    cfdemCloudEnergy& sm
)
:
    YagiWallHT(dict,sm),
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
    partHeatFluxCoeff_(NULL)
{
    allocateMyArrays();

    // no limiting necessary for implicit heat transfer
    maxSource_ = 1e30;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

YagiWallHTImplicit::~YagiWallHTImplicit()
{
    particleCloud_.dataExchangeM().destroy(partHeatFluxCoeff_,1);
}

// * * * * * * * * * * * * * * * private Member Functions  * * * * * * * * * * * * * //
void YagiWallHTImplicit::allocateMyArrays() const
{
//    YagiWallHT::allocateMyArrays();
    double initVal=0.0;
    particleCloud_.dataExchangeM().allocateArray(partHeatFluxCoeff_,initVal,1);
}
// * * * * * * * * * * * * * * * * Member Fct  * * * * * * * * * * * * * * * //

void YagiWallHTImplicit::calcEnergyContribution()
{
    allocateMyArrays();

    YagiWallHT::calcEnergyContribution();

    QPartFluidCoeff_.primitiveFieldRef() = 0.0;

    particleCloud_.averagingM().setScalarSum
    (
        QPartFluidCoeff_,
        partHeatFluxCoeff_,
        particleCloud_.particleWeights(),
        NULL
    );

    QPartFluidCoeff_.primitiveFieldRef() /= -QPartFluidCoeff_.mesh().V();

//    QPartFluidCoeff_.correctBoundaryConditions();

}

void YagiWallHTImplicit::addEnergyCoefficient(volScalarField& Qsource) const
{
    Qsource += QPartFluidCoeff_;
}

void YagiWallHTImplicit::heatFlux(label index, scalar h, scalar As, scalar Tfluid)
{
    partHeatFlux_[index][0] = -h * As * partTemp_[index][0];
}

void YagiWallHTImplicit::heatFluxCoeff(label index, scalar h, scalar As)
{
   partHeatFluxCoeff_[index][0] = h * As;
}

void YagiWallHTImplicit::giveData(int call)
{
    if(call == 1)
    {
        //Info << "total convective particle-fluid heat flux [W] (Eulerian) = " << gSum(QPartFluid_*1.0*QPartFluid_.mesh().V()) << endl;

        particleCloud_.dataExchangeM().giveData(partHeatFluxName_,"scalar-atom", partHeatFlux_);
    }
}

void YagiWallHTImplicit::postFlow()
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
                    Tfluid = tempField_[cellI];

                Tpart = partTemp_[index][0];
                partHeatFlux_[index][0] = (Tfluid - Tpart) * partHeatFluxCoeff_[index][0];
            }
    }

    giveData(1);
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
