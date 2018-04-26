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

#include "BeetstraDragPoly.H"
#include "addToRunTimeSelectionTable.H"
#include "averagingModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(BeetstraDragPoly, 0);

addToRunTimeSelectionTable
(
    forceModel,
    BeetstraDragPoly,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
BeetstraDragPoly::BeetstraDragPoly
(
    const dictionary& dict,
    cfdemCloud& sm
)
:
    BeetstraDrag(dict,sm),
    fines_(propsDict_.lookupOrDefault<bool>("fines",false)),
    dSauter_(sm.mesh().lookupObject<volScalarField> ("dSauter"))
{
    if (fines_)
    {
        volScalarField& alphaSt(const_cast<volScalarField&>(sm.mesh().lookupObject<volScalarField> ("alphaSt")));
        alphaSt_.set(&alphaSt);
        volScalarField& dSauterMix(const_cast<volScalarField&>(sm.mesh().lookupObject<volScalarField> ("dSauterMix")));
        dSauterMix_.set(&dSauterMix);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

BeetstraDragPoly::~BeetstraDragPoly()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void BeetstraDragPoly::adaptVoidfraction(double& voidfraction, label cellI) const
{
    voidfraction -= alphaSt_()[cellI];
    if (voidfraction < minVoidfraction_) voidfraction = minVoidfraction_;
}

scalar BeetstraDragPoly::effDiameter(double d, double voidfraction, label cellI, label index) const
{
    scalar dS = dSauter_[cellI];
    scalar effD = d*d/dS;
    if (fines_)
    {
        scalar pureVoidfraction = voidfraction_[cellI];
        scalar dSmix = dSauterMix_()[cellI];
        effD *= pureVoidfraction / voidfraction * (1 - voidfraction) / (1 - pureVoidfraction);
        effD *= dS * dS / (dSmix * dSmix);
    }
    if (particleCloud_.getParticleEffVolFactors()) 
    {
        scalar effVolFac = particleCloud_.particleEffVolFactor(index);
        effD *= effVolFac;
    }
    return effD;
    
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
