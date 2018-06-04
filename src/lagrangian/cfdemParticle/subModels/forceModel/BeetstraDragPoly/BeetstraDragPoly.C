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
    dFine_(1.0)
{
    // if fines are present, take mixture dSauter, otherwise normal dSauter
    if (fines_)
    {
        dFine_ = readScalar(propsDict_.lookup("dFine"));
        volScalarField& alphaP(const_cast<volScalarField&>(sm.mesh().lookupObject<volScalarField> ("alphaP")));
        alphaP_.set(&alphaP);
        volScalarField& alphaSt(const_cast<volScalarField&>(sm.mesh().lookupObject<volScalarField> ("alphaSt")));
        alphaSt_.set(&alphaSt);
        volScalarField& dSauter(const_cast<volScalarField&>(sm.mesh().lookupObject<volScalarField> ("dSauterMix")));
        dSauter_.set(&dSauter);
    }
    else
    {
        volScalarField& dSauter(const_cast<volScalarField&>(sm.mesh().lookupObject<volScalarField> ("dSauter")));
        dSauter_.set(&dSauter);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

BeetstraDragPoly::~BeetstraDragPoly()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void BeetstraDragPoly::adaptVoidfraction(double& voidfraction, label cellI) const
{
    if (fines_) voidfraction -= alphaSt_()[cellI];
    if (voidfraction < minVoidfraction_) voidfraction = minVoidfraction_;
}

scalar BeetstraDragPoly::effDiameter(double d, label cellI, label index) const
{
    scalar dS = dSauter_()[cellI];
    scalar effD = d*d / dS + 0.064*d*d*d*d / (dS*dS*dS);
    
    if (fines_)
    {
        scalar fineCorr = dFine_*dFine_ / dS + 0.064*dFine_*dFine_*dFine_*dFine_ / (dS*dS*dS);
        fineCorr *= d*d*d / (dFine_*dFine_*dFine_) * alphaSt_()[cellI] / alphaP_()[cellI];
        effD += fineCorr;
    }
    
    if (particleCloud_.getParticleEffVolFactors()) 
    {
        scalar effVolFac = particleCloud_.particleEffVolFactor(index);
        effD *= effVolFac;
    }

    return effD;
}

scalar BeetstraDragPoly::meanSauterDiameter(double d, label cellI) const
{
    return dSauter_()[cellI];
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
