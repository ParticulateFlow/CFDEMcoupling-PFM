/*---------------------------------------------------------------------------*\
    CFDEMcoupling - Open Source CFD-DEM coupling

    CFDEMcoupling is part of the CFDEMproject
    www.cfdem.com
                                Christoph Goniva, christoph.goniva@cfdem.com
                                Copyright 2009-2012 JKU Linz
                                Copyright 2012-     DCS Computing GmbH, Linz
-------------------------------------------------------------------------------
License
    This file is part of CFDEMcoupling.

    CFDEMcoupling is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 3 of the License, or (at your
    option) any later version.

    CFDEMcoupling is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with CFDEMcoupling; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Description
    This code is designed to realize coupled CFD-DEM simulations using LIGGGHTS
    and OpenFOAM(R). Note: this code is not part of OpenFOAM(R) (see DISCLAIMER).
\*---------------------------------------------------------------------------*/

#include "error.H"

#include "explicitCouple.H"
#include "addToRunTimeSelectionTable.H"
#include "forceModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(explicitCouple, 0);

addToRunTimeSelectionTable
(
    momCoupleModel,
    explicitCouple,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
explicitCouple::explicitCouple
(
    const dictionary& dict,
    cfdemCloud& sm
)
:
    momCoupleModel(dict,sm),
    propsDict_(dict.subDict(typeName + "Props")),
    fPrev_
    (   IOobject
        (
            "fPrev",
            sm.mesh().time().timeName(),
            sm.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        sm.mesh(),
        dimensionedVector("zero", dimensionSet(1,-2,-2,0,0), vector::zero), // N/m3
        "zeroGradient"
    ),
    fNext_
    (   IOobject
        (
            "fNext",
            sm.mesh().time().timeName(),
            sm.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        sm.mesh(),
        dimensionedVector("zero", dimensionSet(1,-2,-2,0,0), vector::zero), // N/m3
        "zeroGradient"
    ),
    fLimit_(1e10,1e10,1e10)
{
    if (propsDict_.found("fLimit"))
    {
        fLimit_=vector(propsDict_.lookup ("fLimit"));
        Info << "explicit momentum exchange field is limited to : " << fLimit_ << endl;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

explicitCouple::~explicitCouple()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
tmp<volVectorField> explicitCouple::expMomSource()
{
    const scalar tsf = particleCloud_.dataExchangeM().timeStepFraction();

    // calc Ksl

    // update KslNext in first subTS
    // NOTE: without following if we could update f every subTS (based on current values) and use this value
    if(tsf < particleCloud_.mesh().time().deltaT().value()/particleCloud_.dataExchangeM().couplingTime() + 0.000001 )
    {
        // calc fNext
        forAll(fNext_,cellI)
        {
            fNext_[cellI] = arrayToField(cellI);

            // limiter
            for (int i=0;i<3;i++)
            {
                scalar magF = mag(fNext_[cellI][i]);
                if (magF > fLimit_[i]) fNext_[cellI][i] *= fLimit_[i]/magF;
            }
        }
    }

    return tmp<volVectorField>
    (
        new volVectorField("f_explicitCouple", (1. - tsf) * fPrev_ + tsf * fNext_)
    );
}

void explicitCouple::resetMomSourceField()
{
    fPrev_.ref() = fNext_.ref();
    fNext_.primitiveFieldRef() = vector::zero;
}

inline vector explicitCouple::arrayToField(label cellI) const
{
    return particleCloud_.forceM(0).expParticleForces()[cellI] / particleCloud_.mesh().V()[cellI];
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
