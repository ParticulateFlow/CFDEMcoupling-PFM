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

Description
    calculates the correlation between particle momentum and diameter

SourceFiles
    pdCorrelation.C

Contributing Author
    2018    Paul Kieckhefen, TUHH 
\*---------------------------------------------------------------------------*/

#include "error.H"

#include "pdCorrelation.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(pdCorrelation, 0);

addToRunTimeSelectionTable
(
    forceModel,
    pdCorrelation,
    dictionary
);

const dimensionSet dimMomentum(dimForce*dimTime);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
pdCorrelation::pdCorrelation
(
    const dictionary& dict,
    cfdemCloud& sm
)
:
    forceModel(dict,sm),
    propsDict_(dict.subDict(typeName + "Props")),
    d_(nullptr),
    p_(nullptr),
    d2_(nullptr),
    pd_(nullptr),
    cg3_(nullptr),
    dField_
    (   IOobject
        (
            "dField",
            sm.mesh().time().timeName(),
            sm.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        sm.mesh(),
        dimensionedScalar("zero", dimLength, 0),
        "zeroGradient"
    ),
    pField_
    (   IOobject
        (
            "pField",
            sm.mesh().time().timeName(),
            sm.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        sm.mesh(),
        dimensionedVector("zero", dimMomentum, vector::zero)
    ),
    d2Field_
    (   IOobject
        (
            "d2Field",
            sm.mesh().time().timeName(),
            sm.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        sm.mesh(),
        dimensionedScalar("zero", dimLength*dimLength, 0)
    ),
    pdField_
    (   IOobject
        (
            "pdCorrelation",
            sm.mesh().time().timeName(),
            sm.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        sm.mesh(),
        dimensionedVector("zero", dimMomentum*dimLength, vector::zero),
        "zeroGradient"
    ),
    cg3Field_
    (   IOobject
        (
            "cg3Field",
            sm.mesh().time().timeName(),
            sm.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        sm.mesh(),
        dimensionedScalar("zero", dimless, 0)
    ),
    typeCG_
    (
        propsDict_.lookupOrDefault
        (
            "coarseGrainingFactors",
            scalarList(1, sm.cg())
        )
    ),
    particleDensities_
    (
        propsDict_.lookupOrDefault
        (
            "particleDensities",
            scalarList(1, -1.)
        )
    ),
    constantCG_(typeCG_.size() < 2),
    CG_(!constantCG_ || typeCG_[0] > 1. + SMALL),
    runOnWriteOnly_(propsDict_.lookupOrDefault("runOnWriteOnly", false))
{
    if ((particleDensities_[0] < 0) && !particleCloud_.getParticleDensities())
    {
        FatalError<< "Please set the particle density either in LIGGGHTS"
                  << "or the pdCorrelationPropsDict."
                  << "In the first case, set getParticleDensities in the couplingProperties."
                  << abort(FatalError);
    }

    allocateMyArrays();

    dField_.write();
    pdField_.write();


    // init force sub model
    setForceSubModels(propsDict_);
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

pdCorrelation::~pdCorrelation()
{
    particleCloud_.dataExchangeM().destroy(cg3_, 1);
    particleCloud_.dataExchangeM().destroy(d_,  1);
    particleCloud_.dataExchangeM().destroy(p_,  3);
    particleCloud_.dataExchangeM().destroy(d2_, 1);
    particleCloud_.dataExchangeM().destroy(pd_, 3);
}

// * * * * * * * * * * * * * * * private Member Functions  * * * * * * * * * * * * * //
void pdCorrelation::allocateMyArrays() const
{
    // get memory for 2d arrays
    double initVal = 0.0;
    particleCloud_.dataExchangeM().allocateArray(d_,  initVal, 1);
    particleCloud_.dataExchangeM().allocateArray(p_,  initVal, 3);
    particleCloud_.dataExchangeM().allocateArray(d2_, initVal, 1);
    particleCloud_.dataExchangeM().allocateArray(pd_, initVal, 3);
    particleCloud_.dataExchangeM().allocateArray(cg3_, initVal, 1);
}
// * * * * * * * * * * * * * * * public Member Functions  * * * * * * * * * * * * * //

void pdCorrelation::setForce() const
{
    const fvMesh& mesh = particleCloud_.mesh();

    if (runOnWriteOnly_ && !mesh.write()) return; // skip if it's not write time

    allocateMyArrays();

    const Switch densityFromList
    (
        particleCloud_.getParticleTypes()
        && !particleCloud_.getParticleDensities()
    );
    const scalar piOverSix = constant::mathematical::pi / 6.;

    dField_.primitiveFieldRef()   = 0.0;
    pField_.primitiveFieldRef()   = vector::zero;
    d2Field_.primitiveFieldRef()  = 0.0;
    pdField_.primitiveFieldRef()  = vector::zero;
    cg3Field_.primitiveFieldRef() = 0.;

    label  celli(0);
    scalar dp(0.);
    scalar cg(typeCG_[0]);
    scalar rhop(particleDensities_[0]);
    label  typei(0);
    scalar particleV(0.);
    for(int pi = 0; pi < particleCloud_.numberOfParticles(); ++pi)
    {
        celli = particleCloud_.cellIDs()[pi][0];
        if (celli < 0) continue;

        dp = particleCloud_.d(pi);
        if (particleCloud_.getParticleTypes())
        {
            typei = particleCloud_.particleTypes()[pi][0] - 1;
            rhop  = densityFromList ? particleDensities_[typei]
                                    : particleCloud_.particleDensity(pi);
            cg    = constantCG_     ? typeCG_[0] : typeCG_[typei];
        }

        particleV  = piOverSix * rhop * dp * dp * dp;
        d_ [pi][0] = cg * cg * dp;
        d2_[pi][0] = cg * dp * dp;
        p_ [pi][0] = particleV * particleCloud_.velocities()[pi][0];
        p_ [pi][1] = particleV * particleCloud_.velocities()[pi][1];
        p_ [pi][2] = particleV * particleCloud_.velocities()[pi][2];
        pd_[pi][0] = p_[pi][0] * dp / cg;
        pd_[pi][1] = p_[pi][1] * dp / cg;
        pd_[pi][2] = p_[pi][2] * dp / cg;

        cg3_[pi][0] = CG_ ? cg * cg * cg : 1.;
    }


    particleCloud_.averagingM().setScalarSum
    (
        dField_,
        d_,
        particleCloud_.particleWeights(),
        nullptr
    );
    particleCloud_.averagingM().setVectorSum
    (
        pField_,
        p_,
        particleCloud_.particleWeights(),
        nullptr
    );
    particleCloud_.averagingM().setScalarSum
    (
        d2Field_,
        d2_,
        particleCloud_.particleWeights(),
        nullptr
    );
    particleCloud_.averagingM().setVectorSum
    (
        pdField_,
        pd_,
        particleCloud_.particleWeights(),
        nullptr
    );

    particleCloud_.averagingM().setScalarSum
    (
        cg3Field_,
        cg3_,
        particleCloud_.particleWeights(),
        nullptr
    );

    scalar oneOverCG3(1.);
    const scalar oneMinusSmall(1.-SMALL);
    forAll(cg3Field_, celli)
    {
        if (cg3Field_[celli] < oneMinusSmall) continue;
        oneOverCG3 = 1./cg3Field_[celli];

        dField_ [celli] *= oneOverCG3;
        pField_ [celli] *= oneOverCG3;
        d2Field_[celli] *= oneOverCG3;
        pdField_[celli] *= oneOverCG3;
    }

    scalar denom(0.);
    forAll(dField_, celli)
    {
        if (dField_[celli] < SMALL)
        {
            pdField_[celli] = vector::zero;
            continue;
        }
        denom = d2Field_[celli] - dField_[celli] * dField_[celli];
        if (denom < SMALL)
        {
            pdField_[celli] = vector::zero;
            continue;
        }

        pdField_[celli] -= dField_[celli] * pField_[celli];
        pdField_[celli] /= denom;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
