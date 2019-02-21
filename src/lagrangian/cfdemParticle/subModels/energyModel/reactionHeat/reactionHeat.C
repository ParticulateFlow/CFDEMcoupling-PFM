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
#include "reactionHeat.H"
#include "addToRunTimeSelectionTable.H"
#include "dataExchangeModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(reactionHeat, 0);

addToRunTimeSelectionTable(energyModel, reactionHeat, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
reactionHeat::reactionHeat
(
    const dictionary& dict,
    cfdemCloudEnergy& sm
)
:
    energyModel(dict,sm),
    propsDict_(dict.subDict(typeName + "Props")),
    interpolation_(propsDict_.lookupOrDefault<bool>("interpolation",false)),
    verbose_(propsDict_.lookupOrDefault<bool>("verbose",false)),
    mesh_(sm.mesh()),
    maxSource_(1e30),
    reactionHeatName_(propsDict_.lookupOrDefault<word>("reactionHeatName","reactionHeat")),
    reactionHeat_(NULL),
    reactionHeatField_
    (
        IOobject
        (
            "reactionHeatField",
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero", dimensionSet(1,-1,-3,0,0,0,0),0.0)
    )
{
    allocateMyArrays();

    if(propsDict_.found("maxsource"))
    {
        maxSource_ = readScalar(propsDict_.lookup("maxSource"));
        Info << "limiting eulerian field to:" << maxSource_ << endl;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

reactionHeat::~reactionHeat()
{
    particleCloud_.dataExchangeM().destroy(reactionHeat_,1);
}

// * * * * * * * * * * * * * * * private Member Functions  * * * * * * * * * * * * * //
void reactionHeat::allocateMyArrays() const
{
    // get memory for 2d arrays
    double initVal=0.0;
    particleCloud_.dataExchangeM().allocateArray(reactionHeat_,initVal,1);
}

// * * * * * * * * * * * * * * * * Member Fct  * * * * * * * * * * * * * * * //

void reactionHeat::calcEnergyContribution()
{
   // realloc the arrays
    allocateMyArrays();

    particleCloud_.dataExchangeM().getData(reactionHeatName_,"scalar-atom",reactionHeat_);

    for(int index = 0;index < particleCloud_.numberOfParticles(); ++index)
    {
        if (verbose_ && index>=0 && index < 2)
        {
            Pout << "reactionHeat = " << reactionHeat_[index][0] << endl;
        }
    }

    reactionHeatField_.primitiveFieldRef() = 0.0;
    reactionHeatField_.boundaryFieldRef() = 0.0;

    particleCloud_.averagingM().setScalarSum
    (
        reactionHeatField_,
        reactionHeat_,
        particleCloud_.particleWeights(),
        NULL
    );

    reactionHeatField_.primitiveFieldRef() /= (reactionHeatField_.mesh().V());

    forAll(reactionHeatField_,cellI)
    {
        scalar EuFieldInCell = reactionHeatField_[cellI];

        if(mag(EuFieldInCell) > maxSource_ )
        {
            reactionHeatField_[cellI] = sign(EuFieldInCell) * maxSource_;
        }
    }

    reactionHeatField_.correctBoundaryConditions();
}

void reactionHeat::addEnergyContribution(volScalarField& Qsource) const
{
    Qsource += reactionHeatField_;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
