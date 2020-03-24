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
#include "SyamlalDiffCoeff.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(SyamlalDiffCoeff, 0);

addToRunTimeSelectionTable
(
    diffCoeffModel,
    SyamlalDiffCoeff,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
SyamlalDiffCoeff::SyamlalDiffCoeff
(
    const dictionary& dict,
    cfdemCloud& sm
)
:
    diffCoeffModel(dict,sm),
    propsDict_(dict.subDict(typeName + "Props")),
    voidfractionFieldName_(propsDict_.lookupOrDefault<word>("voidfractionFieldName","voidfraction")),
    voidfraction_(sm.mesh().lookupObject<volScalarField> (voidfractionFieldName_))
{

}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

SyamlalDiffCoeff::~SyamlalDiffCoeff()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<volScalarField> SyamlalDiffCoeff::diffCoeff() const
{
	const volScalarField& D0Field_ = D0Field();

    tmp<volScalarField> tvf
    (
        new volScalarField
        (
            IOobject
            (
                "tmpDiffCoeff",
                voidfraction_.instance(),
                voidfraction_.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            voidfraction_.mesh(),
            dimensionedScalar("zero", dimensionSet(0,2,-1,0,0,0,0), 0.0)
        )
    );

    volScalarField& svf = tvf.ref();

    forAll(svf,cellI)
    {
        scalar D0 = D0Field_[cellI];
      
        if (1-voidfraction_[cellI] < SMALL) svf[cellI] = D0;
        else if (voidfraction_[cellI] < SMALL) svf[cellI] = 0.0;
        else svf[cellI] = (1-sqrt(1-voidfraction_[cellI]+SMALL)) / (voidfraction_[cellI]) * D0;
    }
    return tvf;
}



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
