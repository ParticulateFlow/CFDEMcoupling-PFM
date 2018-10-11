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
    wallQFactorName_(propsDict_.lookupOrDefault<word>("wallQFactorName","wallQFactor")),
    wallQFactor_
    (   IOobject
        (
            wallQFactorName_,
            sm.mesh().time().timeName(),
            sm.mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        sm.mesh(),
        dimensionedScalar("zero", dimensionSet(0,0,0,0,0,0,0), 1.0)
    ),
    hasWallQFactor_(false)
{
    if (wallQFactor_.headerOk())
    {
        Info << "Found field for scaling wall heat flux.\n" << endl;
        hasWallQFactor_ = true;
    }
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

    volScalarField& svf = tvf.ref();


    forAll(svf,cellI)
    {
        if (1-voidfraction_[cellI] < SMALL) svf[cellI] = kf0_.value();
        else if (voidfraction_[cellI] < SMALL) svf[cellI] = 0.0;
        else svf[cellI] = (1-sqrt(1-voidfraction_[cellI]+SMALL)) / (voidfraction_[cellI]) * kf0_.value();
    }

    // if a wallQFactor field is present, use it to scale heat transport through a patch
    if (hasWallQFactor_)
    {
        wallQFactor_.correctBoundaryConditions();
        forAll(wallQFactor_.boundaryField(), patchi)
            svf.boundaryFieldRef()[patchi] *= wallQFactor_.boundaryField()[patchi];
    }

    return tvf;
}

tmp<volScalarField> SyamlalThermCond::thermDiff() const
{
    return thermCond()/(rho_*Cp_);
}



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
