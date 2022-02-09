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
#include "expParticleForces.H"
#include "mathExtra.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(expParticleForces, 0);

addToRunTimeSelectionTable(otherForceModel,expParticleForces, dictionary);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
expParticleForces::expParticleForces
(
    const dictionary& dict,
    cfdemCloud& sm
)
:
    otherForceModel(dict,sm)//,
    //propsDict_(dict.subDict(typeName + "Props"))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

expParticleForces::~expParticleForces()
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<volVectorField> expParticleForces::exportForceField()
{
    tmp<volVectorField> tsource
    (
        new volVectorField
        (
            IOobject
            (
                "expParticleForces",
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

    volVectorField& source = tsource.ref();

    // negative sign in sum because force on particles = - force on fluid
    for(int i=0; i<particleCloud_.nrMomCoupleModels(); i++)
        if (particleCloud_.momCoupleM(i).type() == "explicitCouple")
            source -= particleCloud_.momCoupleM(i).expMomSource();

    return tsource;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
