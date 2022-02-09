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
    multiTypes_(false),
    partReactantName_("reactantPerParticle"),
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
    Nevery_(propsDict_.lookupOrDefault<label>("Nevery",1)),
    maxTypeCG_(1),
    typeCG_(propsDict_.lookupOrDefault<scalarList>("coarseGrainingFactors",scalarList(1,1.0))),
    scaleDia_(1.)
{
    particleCloud_.registerParticleProperty<double**>(partReactantName_,1);

    particleCloud_.checkCG(true);
    if (propsDict_.found("scale") && typeCG_.size()==1)
    {
        // if "scale" is specified and there's only one single type, use "scale"
        scaleDia_=scalar(readScalar(propsDict_.lookup("scale")));
        typeCG_[0] = scaleDia_;
    }

    if (typeCG_.size()>1)
    {
        multiTypes_ = true;
        maxTypeCG_ = typeCG_.size();
    }
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

reactantPerParticle::~reactantPerParticle()
{
}

// * * * * * * * * * * * * * * * private Member Functions  * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Member Fct  * * * * * * * * * * * * * * * //

void reactantPerParticle::execute()
{
    loopCounter_++;
    if (loopCounter_ % Nevery_ != 0)
    {
        return;
    }

    particlesPerCell_ *= 0.0;

    label  cellI=0;
    scalar voidfraction = 1.0;
    scalar cellvolume = 0.0;
    scalar particlesPerCell = 1.0;
    scalar cg = 1.0;
    label partType = 1;

    if (particleCloud_.cg() > 1)
    {
        scaleDia_ = particleCloud_.cg();
    }
    scalar scaleDia3 = scaleDia_*scaleDia_*scaleDia_;

    double**& reactantPerParticle_ = particleCloud_.getParticlePropertyRef<double**>(partReactantName_);

    // first create particles per cell field
    for (int index=0; index<particleCloud_.numberOfParticles(); ++index)
    {
        cellI=particleCloud_.cellIDs()[index][0];
        if (cellI >= 0)
        {

            if (multiTypes_)
            {
                partType = particleCloud_.particleType(index);
                if (partType > maxTypeCG_)
                {
                    FatalError<< "Too few coarse-graining factors provided." << abort(FatalError);
                }
                cg = typeCG_[partType - 1];
                scaleDia3 = cg*cg*cg;
            }
            particlesPerCell_[cellI] += scaleDia3;
        }
    }

    // fill array and communicate it
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
    particleCloud_.dataExchangeM().giveData(partReactantName_, "scalar-atom", reactantPerParticle_);

    Info << "give data done" << endl;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
