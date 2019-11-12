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

#include "particleDeformation.H"
#include "addToRunTimeSelectionTable.H"
#include "dataExchangeModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(particleDeformation, 0);

addToRunTimeSelectionTable
(
    forceModel,
    particleDeformation,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
particleDeformation::particleDeformation
(
    const dictionary& dict,
    cfdemCloud& sm
)
:
    forceModel(dict,sm),
    propsDict_(dict.subDict(typeName + "Props")),
    initialExec_(true),
    refFieldName_(propsDict_.lookup("refFieldName")),
    refField_(),
    partType_(propsDict_.lookupOrDefault<label>("partType",0)),
    lowerBound_(readScalar(propsDict_.lookup ("lowerBound"))),
    upperBound_(readScalar(propsDict_.lookup ("upperBound"))),
    partDeformations_(NULL)
{
    allocateMyArrays();

    // init force sub model
    setForceSubModels(propsDict_);

    // define switches which can be read from dict
    forceSubM(0).setSwitchesList(SW_VERBOSE,true); // activate search for verbose switch
    forceSubM(0).setSwitchesList(SW_INTERPOLATION,true); // activate search for interpolate switch

    // read those switches defined above, if provided in dict
    forceSubM(0).readSwitches();

    particleCloud_.checkCG(false);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

particleDeformation::~particleDeformation()
{
    particleCloud_.dataExchangeM().destroy(partDeformations_,1);
}

// * * * * * * * * * * * * * * * private Member Functions  * * * * * * * * * * * * * //
void particleDeformation::allocateMyArrays() const
{
    // get memory for 2d arrays
    double initVal = 0.0;
    particleCloud_.dataExchangeM().allocateArray(partDeformations_,initVal,1);
}
// * * * * * * * * * * * * * * * public Member Functions  * * * * * * * * * * * * * //

void particleDeformation::setForce() const
{
    if (initialExec_)
    {
        init();
        initialExec_ = false;
    }
    // realloc the arrays
    allocateMyArrays();

    label cellI = 0;
    label partType = -1;
    scalar refFieldValue = 0.0;
    scalar deformationDegree = 0.0;
   
    interpolationCellPoint<scalar> refFieldInterpolator_(refField_());

    for(int index = 0; index < particleCloud_.numberOfParticles(); ++index)
    {
        cellI = particleCloud_.cellIDs()[index][0];
        partType = particleCloud_.particleType(index);
        if (cellI >= 0 && partType == partType_)
        {
            if (forceSubM(0).interpolation())
            {
                vector position = particleCloud_.position(index);
                refFieldValue = refFieldInterpolator_.interpolate(position,cellI);
            }
            else
            {
                refFieldValue = refField_()[cellI];
            }

            if (refFieldValue <= lowerBound_)
            {
                deformationDegree = 0.0;
            }
            else if (refFieldValue >= upperBound_)
            {
                deformationDegree = 1.0;
            }
            else
            {
                deformationDegree = (refFieldValue - lowerBound_) / (upperBound_ - lowerBound_);
            }

            partDeformations_[index][0] = deformationDegree;


            if(forceSubM(0).verbose() && index >= 0 && index < 2)
            {
                Info << "refFieldValue = " << refFieldValue << endl;
                Info << "deformationDegree = " << deformationDegree << endl;
            }
        }
    }

    // give DEM data
    particleCloud_.dataExchangeM().giveData("partDeformations","scalar-atom", partDeformations_);
}

void particleDeformation::init() const
{
    // check if ref field with corresponding name has been read by some other class or if it needs to be newly created
    if (particleCloud_.mesh().foundObject<volScalarField> (refFieldName_))
    {
   //     volScalarField& refField(particleCloud_.mesh().lookupObject<volScalarField> (refFieldName_));
        volScalarField& refField(const_cast<volScalarField&>(particleCloud_.mesh().lookupObject<volScalarField> (refFieldName_)));
        refField_.set(&refField);
    }
    else
    {
        // get start time to read field from
        scalar tstart = particleCloud_.mesh().time().startTime().value();
        refField_.set
        (
            new volScalarField
            (
                IOobject
                (
                    refFieldName_,
                    particleCloud_.mesh().time().timeName(tstart),
                    particleCloud_.mesh(),
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                ),
                particleCloud_.mesh()
            )
        );
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
