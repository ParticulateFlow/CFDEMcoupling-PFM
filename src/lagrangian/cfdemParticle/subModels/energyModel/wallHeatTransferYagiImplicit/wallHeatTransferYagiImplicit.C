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
#include "wallHeatTransferYagiImplicit.H"
#include "addToRunTimeSelectionTable.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(wallHeatTransferYagiImplicit, 0);

addToRunTimeSelectionTable(energyModel, wallHeatTransferYagiImplicit, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
wallHeatTransferYagiImplicit::wallHeatTransferYagiImplicit
(
    const dictionary& dict,
    cfdemCloudEnergy& sm
)
:
    wallHeatTransferYagi(dict,sm),
    QWallFluidCoeffName_(propsDict_.lookupOrDefault<word>("QWallFluidCoeffName","QWallFluidCoeff")),
    QWallFluidCoeff_
    (   IOobject
        (
            QWallFluidCoeffName_,
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

wallHeatTransferYagiImplicit::~wallHeatTransferYagiImplicit()
{
    particleCloud_.dataExchangeM().destroy(partHeatFluxCoeff_,1);
}

// * * * * * * * * * * * * * * * private Member Functions  * * * * * * * * * * * * * //
void wallHeatTransferYagiImplicit::allocateMyArrays() const
{
//    wallHeatTransferYagi::allocateMyArrays();
    double initVal=0.0;
    particleCloud_.dataExchangeM().allocateArray(partHeatFluxCoeff_,initVal,1);
}
// * * * * * * * * * * * * * * * * Member Fct  * * * * * * * * * * * * * * * //

void wallHeatTransferYagiImplicit::calcEnergyContribution()
{
    allocateMyArrays();

    QWallFluidCoeff_.primitiveFieldRef() = 0.0;

    wallHeatTransferYagi::calcEnergyContribution();

    QWallFluidCoeff_.primitiveFieldRef() /= QWallFluidCoeff_.mesh().V();

//    QWallFluidCoeff_.correctBoundaryConditions();

}

void wallHeatTransferYagiImplicit::addEnergyCoefficient(volScalarField& Qsource) const
{
    Qsource += QWallFluidCoeff_;
}

void wallHeatTransferYagiImplicit::heatFlux(label faceCelli, scalar H, scalar area, scalar Twall, scalar Tfluid)
{
    QWallFluid_[faceCelli] += H * area * Twall;
}

void wallHeatTransferYagiImplicit::heatFluxCoeff(label faceCelli, scalar H, scalar area)
{
    QWallFluidCoeff_[faceCelli] -= H * area;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
