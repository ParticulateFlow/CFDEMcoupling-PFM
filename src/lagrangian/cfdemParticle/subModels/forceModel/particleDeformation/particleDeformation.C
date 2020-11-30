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
#include "OFstream.H"

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
    voidfractionFieldName_(propsDict_.lookupOrDefault<word>("voidfractionFieldName","voidfraction")),
    voidfraction_(sm.mesh().lookupObjectRef<volScalarField> (voidfractionFieldName_)),
    defaultDeformCellsName_(propsDict_.lookupOrDefault<word>("defaultDeformCellsName","none")),
    defaultDeformCells_(),
    existDefaultDeformCells_(false),
    defaultDeformation_(propsDict_.lookupOrDefault<scalar>("defaultDeformation",1.0)),
    existBackgroundVoidage_(false),
    backgroundVoidage_(propsDict_.lookupOrDefault<scalar>("backgroundVoidage",-1.0)),
    backgroundRef_(propsDict_.lookupOrDefault<scalar>("backgroundRef",0.35)),
    partTypes_(propsDict_.lookupOrDefault<labelList>("partTypes",labelList(1,-1))),
    lowerBounds_(propsDict_.lookupOrDefault<scalarList>("lowerBounds",scalarList(1,-1.0))),
    upperBounds_(propsDict_.lookupOrDefault<scalarList>("upperBounds",scalarList(1,-1.0))),
    partDeformationsName_("partDeformations"),
    controlForceOnDefPart_(propsDict_.lookupOrDefault<bool>("controlForceOnDefPart",false)),
    controlFieldName_(propsDict_.lookupOrDefault<word>("controlFieldName","")),
    controlField_(NULL),
    controlPoint_(propsDict_.lookupOrDefault<vector>("controlPoint",vector::zero)),
    controlCell_(-1),
    controlTargetValue_(propsDict_.lookupOrDefault<scalar>("controlTargetValue",0.0)),
    controlCouplingStrength_(propsDict_.lookupOrDefault<scalar>("controlCouplingStrength",0.0))
{
    particleCloud_.registerParticleProperty<double**>(partDeformationsName_,1);

    // init force sub model
    setForceSubModels(propsDict_);

    // define switches which can be read from dict
    forceSubM(0).setSwitchesList(SW_VERBOSE,true); // activate search for verbose switch
    forceSubM(0).setSwitchesList(SW_INTERPOLATION,true); // activate search for interpolate switch

    // read those switches defined above, if provided in dict
    forceSubM(0).readSwitches();

    particleCloud_.checkCG(false);

    if(defaultDeformCellsName_ != "none")
    {
       defaultDeformCells_.set(new cellSet(particleCloud_.mesh(),defaultDeformCellsName_));
       existDefaultDeformCells_ = true;
       Info << "particleDeformation: default deformation of " << defaultDeformation_ << " in cellSet " << defaultDeformCells_().name() <<
        " with " << defaultDeformCells_().size() << " cells." << endl;
       if (defaultDeformation_ < 0.0 || defaultDeformation_ > 1.0)
       {
           defaultDeformation_ = min(max(defaultDeformation_,0.0),1.0);
           Info << "Resetting defaultDeformation to range [0;1]" << endl;
       }
    }

    if (backgroundVoidage_ <= 1.0 && backgroundVoidage_ >= 0.01)
    {
        existBackgroundVoidage_ = true;
    }

    // check if only single value instead of list was provided
    if (propsDict_.found("partType"))
    {
        partTypes_[0] = readLabel(propsDict_.lookup("partType"));
    }

    if (propsDict_.found("lowerBound"))
    {
        lowerBounds_[0] = readScalar(propsDict_.lookup("lowerBound"));
    }

    if (propsDict_.found("upperBound"))
    {
        upperBounds_[0] = readScalar(propsDict_.lookup("upperBound"));
    }

    if (partTypes_.size() != lowerBounds_.size() || partTypes_.size() != upperBounds_.size())
    {
        FatalError << "Inconsistent number of particle types and/or bounds provided." << abort(FatalError);
    }

    Info << "partTypes: " << partTypes_ << endl;
    Info << "lowerBounds: " << lowerBounds_ << endl;
    Info << "upperBounds: " << upperBounds_ << endl;

    if (controlForceOnDefPart_)
    {
        controlCell_ = sm.mesh().findCell(controlPoint_);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

particleDeformation::~particleDeformation()
{
}

// * * * * * * * * * * * * * * * private Member Functions  * * * * * * * * * * * * * //

bool particleDeformation::defaultDeformCell(label cell) const
{
    if (!existDefaultDeformCells_) return false;
    else return defaultDeformCells_()[cell];
}
// * * * * * * * * * * * * * * * public Member Functions  * * * * * * * * * * * * * //

void particleDeformation::setForce() const
{
    if (initialExec_)
    {
        init();
        initialExec_ = false;
    }

    scalar controlCurrValue = 0.0;
    scalar scaleFactor = 0.0;
    if (controlForceOnDefPart_)
    {
        controlCurrValue = (*controlField_)[controlCell_];
        reduce(controlCurrValue, sumOp<scalar>());
        scaleFactor = 1.0 + controlCouplingStrength_ * (controlTargetValue_ - controlCurrValue);
        Info << "particleDeformation: scaleFactor = " << scaleFactor << endl;
    }

    double**& partDeformations_ = particleCloud_.getParticlePropertyRef<double**>(partDeformationsName_);

    label cellI = 0;
    label partType = -1;
    scalar refFieldValue = 0.0;
    scalar deformationDegree = 0.0;

    interpolationCellPoint<scalar> refFieldInterpolator_(refField_());

    for(int index = 0; index < particleCloud_.numberOfParticles(); ++index)
    {
        cellI = particleCloud_.cellIDs()[index][0];
        partType = particleCloud_.particleType(index);
        label listIndex = getListIndex(partType);
        if (cellI >= 0 && listIndex >= 0)
        {
            if (defaultDeformCell(cellI))
            {
                deformationDegree = defaultDeformation_;
            }
            else
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

                if (refFieldValue <= lowerBounds_[listIndex])
                {
                    deformationDegree = 0.0;
                }
                else if (refFieldValue >= upperBounds_[listIndex])
                {
                    deformationDegree = 1.0;
                }
                else
                {
                    deformationDegree = (refFieldValue - lowerBounds_[listIndex]) / (upperBounds_[listIndex] - lowerBounds_[listIndex]);
                }
            }

            partDeformations_[index][0] = deformationDegree;

            if (existBackgroundVoidage_ && deformationDegree > 0.1)
            {
                voidfraction_[cellI] = backgroundRef_ + deformationDegree * (backgroundVoidage_ - backgroundRef_);
            }

            if (controlForceOnDefPart_)
            {
                particleCloud_.scaleForce(index,scaleFactor);
            }

            if(forceSubM(0).verbose() && index >= 0 && index < 2)
            {
                Info << "refFieldValue = " << refFieldValue << endl;
                Info << "deformationDegree = " << deformationDegree << endl;
            }
        }
    }

    // give DEM data
    particleCloud_.dataExchangeM().giveData(partDeformationsName_,"scalar-atom", partDeformations_);
}

void particleDeformation::init() const
{
    // check if ref field with corresponding name has been read by some other class or if it needs to be newly created
    if (particleCloud_.mesh().foundObject<volScalarField> (refFieldName_))
    {
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
                    IOobject::AUTO_WRITE
                ),
                particleCloud_.mesh()
            )
        );
    }

    if (controlForceOnDefPart_)
    {
        volScalarField& controlField(const_cast<volScalarField&>(particleCloud_.mesh().lookupObject<volScalarField> (controlFieldName_)));
        controlField_ = &controlField;
    }
}

label particleDeformation::getListIndex(label testElement) const
{
    for(label ind = 0; ind<partTypes_.size(); ind++)
    {
        if (partTypes_[ind] == testElement) return ind;
    }
    return -1;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
