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
#include "SyamlalThermCond.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(SyamlalThermCond, 0);

addToRunTimeSelectionTable
(
    thermCondModel,
    SyamlalThermCond,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
SyamlalThermCond::SyamlalThermCond
(
    const dictionary& dict,
    cfdemCloud& sm
)
:
    thermCondModel(dict,sm),
    propsDict_(dict.subDict(typeName + "Props")),
    voidfractionFieldName_(propsDict_.lookupOrDefault<word>("voidfractionFieldName","voidfraction")),
    voidfraction_(sm.mesh().lookupObject<volScalarField> (voidfractionFieldName_)),
    rhoFieldName_(propsDict_.lookupOrDefault<word>("rhoFieldName","rho")),
    rho_(sm.mesh().lookupObject<volScalarField> (rhoFieldName_)),
    wallLayerFactor_(1.0)
{
    if (propsDict_.found("wallLayerFactor"))
         wallLayerFactor_=readScalar(propsDict_.lookup ("wallLayerFactor"));
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

SyamlalThermCond::~SyamlalThermCond()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<volScalarField> SyamlalThermCond::thermCond() const
{
    tmp<volScalarField> tvf
    (
        new volScalarField
        (
            IOobject
            (
                "tmpThCond",
                voidfraction_.instance(),
                voidfraction_.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            voidfraction_.mesh(),
	    dimensionedScalar("zero", dimensionSet(1,1,-3,-1,0,0,0), 0.0)
        )
    );
    volScalarField& svf = tvf();
    
    svf = (1-sqrt(1-voidfraction_)) / (voidfraction_) * kf0_;
    
    forAll(voidfraction_.boundaryField(), patchi)
        svf.boundaryField()[patchi] *= wallLayerFactor_;
  
    return tvf;
}

tmp<volScalarField> SyamlalThermCond::thermDiff() const
{
    return thermCond()/(rho_*Cp_);
}



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
