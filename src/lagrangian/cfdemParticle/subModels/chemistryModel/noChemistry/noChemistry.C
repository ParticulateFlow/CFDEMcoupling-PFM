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
#include "noChemistry.H"
#include "addToRunTimeSelectionTable.H"

#include "IFstream.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(noChemistry, 0);

addToRunTimeSelectionTable(chemistryModel, noChemistry, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
noChemistry::noChemistry
(
    const dictionary& dict,
    cfdemCloudEnergy& sm
)
:
    chemistryModel(dict,sm)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

noChemistry::~noChemistry()
{}

// * * * * * * * * * * * * * * * private Member Functions  * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Member Fct  * * * * * * * * * * * * * * * //

void noChemistry::execute()
{}

void noChemistry::tryout(const label i)
{
    mi_.set
            (
                i,
                new volScalarField
                (
                    IOobject
                    (
                        "empty1",
                        particleCloud_.mesh().time().timeName(),
                        particleCloud_.mesh(),
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                     ),
                     particleCloud_.mesh(),
                     dimensionedScalar("zero",dimless,0.0)
                )
            );
}

tmp<volScalarField> noChemistry ::Smi(const label i) const
{
    return tmp<volScalarField> (mi_[i]);
}

  /*  tmp<volScalarField> mi_[i]
    (
        new volScalarField
        (
            IOobject
            (
                "empty1",
                particleCloud_.mesh().time().timeName(),
                particleCloud_.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
             ),
             particleCloud_.mesh(),
             dimensionedScalar
             (
              "zero",
              dimMass/dimTime,//dimensionSet(0,2,-1,0,0,0,0),
              0.0
             )
        )
    );

    return mi_[i];
}*/

tmp<volScalarField> noChemistry::Sm() const
{
    tmp<volScalarField> sm_
    (
        new volScalarField
        (
            IOobject
            (
                "empty2",
                particleCloud_.mesh().time().timeName(),
                particleCloud_.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
             ),
             particleCloud_.mesh(),
             dimensionedScalar("zero",dimless,0.0)
        )
    );

    return sm_;
}


/*tmp<Foam::fvScalarMatrix> noChemistry::Smi(const label i) const
tmp<Foam::fvScalarMatrix> noChemistry::Sm() const*/


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
