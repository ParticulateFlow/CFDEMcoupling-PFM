/*---------------------------------------------------------------------------*\
    CFDEMcoupling - Open Source CFD-DEM coupling

    Copyright (C) 2023  Behrad Esgandari, JKU Linz, Austria
-------------------------------------------------------------------------------
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

\*---------------------------------------------------------------------------*/

#include "error.H"
#include "gravityEff.H"
#include "mathExtra.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //
defineTypeNameAndDebug(gravityEff, 0);

addToRunTimeSelectionTable(otherForceModel,gravityEff, dictionary);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
// Construct from components
gravityEff::gravityEff
(
    const dictionary& dict,
    cfdemCloud& sm
)
:
    otherForceModel(dict,sm),
    propsDict_(dict.subDict(typeName + "Props")),
    voidfractionFieldName_(propsDict_.lookupOrDefault<word>("voidfractionFieldName","voidfraction")),
    voidfraction_(sm.mesh().lookupObject<volScalarField> (voidfractionFieldName_)),
    rhoGFieldName_(propsDict_.lookupOrDefault<word>("rhoGFieldName","rho")),
    rhoG_(sm.mesh().lookupObject<volScalarField> (rhoGFieldName_)),
    gravityFieldName_(propsDict_.lookupOrDefault<word>("gravityFieldName","g")),
    g_(sm.mesh().lookupObject<uniformDimensionedVectorField> (gravityFieldName_)),
    rhoPart_("rhoPart",dimensionSet(1,-3,0,0,0),0.0),
    unity_(sm.mesh().lookupObject<volScalarField> ("unity"))
{
    if (propsDict_.found("rhoPart"))
        rhoPart_.value()=readScalar(propsDict_.lookup ("rhoPart"));
    else
        FatalError <<"Please specify density of solid phase.\n" << abort(FatalError);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

gravityEff::~gravityEff()
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
tmp<volVectorField> gravityEff::exportForceField()
{
    tmp<volVectorField> tsource
    (
        new volVectorField
        (
            IOobject
            (
                "grav",
                particleCloud_.mesh().time().timeName(),
                particleCloud_.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            particleCloud_.mesh(),
            dimensionedVector
            (
                "zero",dimensionSet(1, -2, -2, 0, 0),vector::zero
            )
        )
    );

    if (propsDict_.found("rhoPart"))
    {
        //Mixture (gas and solid) density based on whole domain
        dimensionedScalar rhoMix_ = fvc::domainIntegrate((1.0-voidfraction_) * rhoPart_ + voidfraction_ * rhoG_) / fvc::domainIntegrate(unity_);
        volVectorField& source = tsource.ref();
        source = (rhoG_ - rhoMix_) * voidfraction_ * g_;
    }

    return tsource;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
