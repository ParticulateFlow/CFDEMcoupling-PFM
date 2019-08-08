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
                        M.Efe Kinaci, JKU Linz, Austria

\*---------------------------------------------------------------------------*/

#include "error.H"
#include "reactantPerParticle.H"
#include "addToRunTimeSelectionTable.H"

#include "dataExchangeModel.H"
#include "IFstream.H"



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(reactantPerParticle, 0);

addToRunTimeSelectionTable
(
        chemistryModel,
        reactantPerParticle,
        dictionary
);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
reactantPerParticle::reactantPerParticle
(
    const dictionary& dict,
    cfdemCloudEnergy& sm
)
:
    chemistryModel(dict,sm),
    propsDict_(dict.subDict(typeName + "Props")),
    mesh_(sm.mesh()),
    verbose_(propsDict_.lookupOrDefault<bool>("verbose",false)),
    reactantPerParticle_(NULL),
    voidfractionFieldName_(propsDict_.lookupOrDefault<word>("voidfractionFieldName","voidfraction")),
    voidfraction_(sm.mesh().lookupObject<volScalarField>(voidfractionFieldName_)),
    particlesPerCell_
    (   IOobject
        (
            "particlesPerCell",
            sm.mesh().time().timeName(),
            sm.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        sm.mesh(),
        dimensionedScalar("zero", dimensionSet(0,0,0,0,0), 0)
    ),
    loopCounter_(-1),
    Nevery_(propsDict_.lookupOrDefault<label>("Nevery",1))
{
    particleCloud_.checkCG(false);
    allocateMyArrays();
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

reactantPerParticle::~reactantPerParticle()
{
    particleCloud_.dataExchangeM().destroy(reactantPerParticle_,1);
}

// * * * * * * * * * * * * * * * private Member Functions  * * * * * * * * * * * * * //
void reactantPerParticle::allocateMyArrays() const
{
    double initVal=0.0;
    if (particleCloud_.dataExchangeM().maxNumberOfParticles() > 0)
    {
        // get memory for 2d arrays
        particleCloud_.dataExchangeM().allocateArray(reactantPerParticle_,initVal,1,"nparticles");
    }
}

void reactantPerParticle::reAllocMyArrays() const
{
    double initVal=0.0;
    particleCloud_.dataExchangeM().allocateArray(reactantPerParticle_,initVal,1);
}

// * * * * * * * * * * * * * * * * Member Fct  * * * * * * * * * * * * * * * //

void reactantPerParticle::execute()
{
    loopCounter_++;
    if (loopCounter_ % Nevery_ != 0)
    {
        return;
    }
    // realloc the arrays
    reAllocMyArrays();

    particlesPerCell_ *= 0.0;

    label  cellI=0;
    scalar voidfraction(1);
    scalar cellvolume(0.0);
    scalar particlesPerCell(1.0);

    // first create particles per cell field
    for (int index=0; index<particleCloud_.numberOfParticles(); ++index)
    {
        cellI=particleCloud_.cellIDs()[index][0];
        if (cellI >= 0)
        {
            particlesPerCell_[cellI] += 1.0;
        }
    }

    // no fill array and communicate it
    for (int index=0; index<particleCloud_.numberOfParticles(); ++index)
    {
        cellI=particleCloud_.cellIDs()[index][0];
        if (cellI >= 0)
        {
            voidfraction    =   voidfraction_[cellI];
            cellvolume      =   mesh_.V()[cellI];
            particlesPerCell=   particlesPerCell_[cellI];
            reactantPerParticle_[index][0] = voidfraction * cellvolume / particlesPerCell;
        }

        if (verbose_) Info << "reactantPerParticle_" << reactantPerParticle_[index][0] << endl;
    }

        // give DEM data
        particleCloud_.dataExchangeM().giveData("reactantPerParticle", "scalar-atom", reactantPerParticle_);

        Info << "give data done" << endl;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
