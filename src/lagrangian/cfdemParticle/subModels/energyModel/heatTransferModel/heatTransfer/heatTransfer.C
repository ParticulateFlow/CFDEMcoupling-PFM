/*---------------------------------------------------------------------------*\
    CFDEMcoupling - Open Source CFD-DEM coupling

    CFDEMcoupling is part of the CFDEMproject
    www.cfdem.com
                                Christoph Goniva, christoph.goniva@cfdem.com
                                Copyright 2009-2012 JKU Linz
                                Copyright 2012-     DCS Computing GmbH, Linz
-------------------------------------------------------------------------------
License
    This file is part of CFDEMcoupling.

    CFDEMcoupling is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 3 of the License, or (at your
    option) any later version.

    CFDEMcoupling is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with CFDEMcoupling; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Description
    This code is designed to realize coupled CFD-DEM simulations using LIGGGHTS
    and OpenFOAM(R). Note: this code is not part of OpenFOAM(R) (see DISCLAIMER).
\*---------------------------------------------------------------------------*/

#include "error.H"
#include "heatTransfer.H"
#include "addToRunTimeSelectionTable.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(heatTransfer, 0);

addToRunTimeSelectionTable(heatTransferModel, heatTransfer, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
heatTransfer::heatTransfer
(
    const dictionary& dict,
    cfdemCloudEnergy& sm
)
:
    heatTransferModel(dict,sm),
    propsDict_(dict.subDict(typeName + "Props")),
    tempFieldName_(propsDict_.lookup("tempFieldName")),
    tempField_(sm.mesh().lookupObject<volScalarField> (tempFieldName_)),
    voidfractionFieldName_(propsDict_.lookup("voidfractionFieldName")),
    voidfraction_(sm.mesh().lookupObject<volScalarField> (voidfractionFieldName_)),
    maxSource_(1e30),
    velFieldName_(propsDict_.lookup("velFieldName")),
    U_(sm.mesh().lookupObject<volVectorField> (velFieldName_)),
    partTempName_(propsDict_.lookup("partTempName")),
    partTemp_(NULL),
    partHeatFluxName_(propsDict_.lookup("partHeatFluxName")),
    partHeatFlux_(NULL),
    lambda_(readScalar(propsDict_.lookup("lambda"))),
    lambdaPart_(readScalar(propsDict_.lookup("lambdaPart"))),
    Cp_(readScalar(propsDict_.lookup("Cp")))
{
   if (propsDict_.found("maxSource"))
    {
        maxSource_=readScalar(propsDict_.lookup ("maxSource"));
        Info << "limiting eulerian source field to: " << maxSource_ << endl;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

heatTransfer::~heatTransfer()
{}

// * * * * * * * * * * * * * * * * Member Fct  * * * * * * * * * * * * * * * //
void heatTransfer::allocateMyArrays() const
{
    // get memory for 2d arrays
    double initVal=0.0;
    particleCloud_.dataExchangeM().allocateArray(partTemp_,initVal,1);  // field/initVal/with/lenghtFromLigghts
    particleCloud_.dataExchangeM().allocateArray(partHeatFlux_,initVal,1);
}

void heatTransfer::partFluidHeatTransfer(volScalarField& EuField) const
{
        // realloc the arrays
    allocateMyArrays();

    // reset Scalar field
    EuField.internalField() = 0.0;

    // get DEM data
    particleCloud_.dataExchangeM().getData(partTempName_,"scalar-atom",partTemp_);

    const volScalarField& nufField = particleCloud_.forceM(0).forceSubM(0).nuField();
    const volScalarField& rhoField = particleCloud_.forceM(0).forceSubM(0).rhoField();

    // calc La based heat flux
    scalar voidfraction(1);
    vector Ufluid(0,0,0);
    scalar Tfluid(0);
    label cellI=0;
    vector Us(0,0,0);
    scalar ds(0);
    scalar nuf(0);
    scalar magUr(0);
    scalar Rep(0);
    scalar Pr(0);
    scalar Nup(0);
    const scalar n = 3.5; // model parameter (found suitable for 3-mm polymer pellets when modelling dilute flows)

    interpolationCellPoint<scalar> voidfractionInterpolator_(voidfraction_);
    interpolationCellPoint<vector> UInterpolator_(U_);
    interpolationCellPoint<scalar> TInterpolator_(tempField_);

    for(int index = 0;index < particleCloud_.numberOfParticles(); ++index)
    {
        //if(particleCloud_.regionM().inRegion()[index][0])
        //{
            cellI = particleCloud_.cellIDs()[index][0];
            if(cellI >= 0)
            {
                if(particleCloud_.forceM(0).forceSubM(0).interpolation())
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
                nuf = nufField[cellI];
                Rep = ds * magUr / nuf;
                Pr = max(SMALL, Cp_ * nuf * rhoField[cellI] / lambda_);

                if (Rep < 200.)
                {
                    Nup = 2. + 0.6 * pow(voidfraction,n) * sqrt(Rep) * pow(Pr,0.33);
                }
                else if (Rep < 1500.)
                {
                    Nup = 2. + (0.5 * sqrt(Rep) + 0.02 * pow(Rep,0.8)) * pow(voidfraction,n) * pow(Pr,0.33);
                }
                else
                {
                    Nup = 2. + 0.000045 * pow(voidfraction,n) * pow(Rep,1.8);
                }

                scalar h = lambda_ * Nup / ds;
                scalar As = ds * ds * M_PI; // surface area of sphere

                // calc convective heat flux [W]
                scalar partHeatFlux = h * As * (Tfluid - partTemp_[index][0]);
                partHeatFlux_[index][0] = partHeatFlux;


                if(particleCloud_.forceM(0).forceSubM(0).verbose() && index >=0 && index <2)
                {
                    Info << "partHeatFlux = " << partHeatFlux << endl;
                    Info << "magUr = " << magUr << endl;
                    Info << "As = " << As << endl;
                    Info << "nuf = " << nuf << endl;
                    Info << "Rep = " << Rep << endl;
                    Info << "Pr = " << Pr << endl;
                    Info << "Nup = " << Nup << endl;
                    Info << "voidfraction = " << voidfraction << endl;
                    Info << "partTemp_[index][0] = " << partTemp_[index][0] << endl;
                    Info << "Tfluid = " << Tfluid << endl  ;
                }
            }
        //}
    }

    particleCloud_.averagingM().setScalarSum
    (
        EuField,
        partHeatFlux_,
        particleCloud_.particleWeights(),
        NULL
    );

    // scale with -1/(Vcell*rho)
    EuField.internalField() /= -EuField.mesh().V();

    // limit source term
    forAll(EuField,cellI)
    {
        scalar EuFieldInCell = EuField[cellI];

        if(mag(EuFieldInCell) > maxSource_ )
        {
             EuField[cellI] = sign(EuFieldInCell) * maxSource_;
        }
    }

   // Info << "total convective particle-fluid heat flux [W] (Eulerian) = " << gSum(EuField*EuField.mesh().V()) << endl;

    // give DEM data
    particleCloud_.dataExchangeM().giveData(partHeatFluxName_,"scalar-atom", partHeatFlux_);
}

scalar heatTransfer::partThermCond() const
{
    return lambdaPart_;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
