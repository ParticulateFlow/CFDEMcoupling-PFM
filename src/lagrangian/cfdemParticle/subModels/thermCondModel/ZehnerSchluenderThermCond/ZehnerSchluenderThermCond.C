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
#include "ZehnerSchluenderThermCond.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(ZehnerSchluenderThermCond, 0);

addToRunTimeSelectionTable
(
    thermCondModel,
    ZehnerSchluenderThermCond,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
ZehnerSchluenderThermCond::ZehnerSchluenderThermCond
(
    const dictionary& dict,
    cfdemCloud& sm
)
:
    thermCondModel(dict,sm),
    propsDict_(dict.subDict(typeName + "Props")),
    voidfractionFieldName_(propsDict_.lookupOrDefault<word>("voidfractionFieldName","voidfraction")),
    voidfraction_(sm.mesh().lookupObject<volScalarField> (voidfractionFieldName_)),
    ks0_(transportProperties_.lookup("ks")),
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

ZehnerSchluenderThermCond::~ZehnerSchluenderThermCond()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<volScalarField> ZehnerSchluenderThermCond::thermCond() const
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
            dimensionedScalar("zero", dimensionSet(1,1,-3,-1,0,0,0), 0.0),
            "zeroGradient"
        )
    );

    volScalarField& svf = tvf.ref();

    scalar A = ks0_.value()/kf0_.value();
    scalar B = 0.0;
    scalar C = 0.0;
    scalar k = 0.0;
    scalar OnemBoA = 0.0;
    scalar voidfraction = 0.0;
    scalar w = 7.26e-3;

    forAll(svf, cellI)
    {
        // debugging
        Pout << "calculating field in cell " << cellI << endl;
        voidfraction = voidfraction_[cellI];
        if(voidfraction > 1.0 - SMALL) svf[cellI] = 0.0;
        else
        {
            B = 1.25 * Foam::pow((1 - voidfraction) / voidfraction, 1.11);
            OnemBoA = 1.0 - B/A;
            C = (A - 1) / (OnemBoA * OnemBoA) * B/A * log(A/B) - (B - 1)/OnemBoA - 0.5 * (B + 1);
            C *= 2.0 / OnemBoA;
            k = Foam::sqrt(1 - voidfraction) * (w * A + (1 - w) * C) * kf0_.value();
            svf[cellI] = k / (1 - voidfraction);
        }
    }

    // debugging
    Pout << "patch types of svf boundary: " << svf.boundaryField().types() << endl;
    svf.correctBoundaryConditions();

    // if a wallQFactor field is present, use it to scale heat transport through a patch
    if (hasWallQFactor_)
    {
        wallQFactor_.correctBoundaryConditions();
        forAll(wallQFactor_.boundaryField(), patchi)
            svf.boundaryFieldRef()[patchi] *= wallQFactor_.boundaryField()[patchi];
    }

    return tvf;
}

tmp<volScalarField> ZehnerSchluenderThermCond::thermDiff() const
{
    FatalError << "ZehnerSchluenderThermCond does not provide thermal diffusivity." << abort(FatalError);
}



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
