/*---------------------------------------------------------------------------*\
    CFDEMcoupling academic - Open Source CFD-DEM coupling
    
    Contributing authors:
    Thomas Lichtenegger
    Copyright (C) 2015- Johannes Kepler University, Linz
-------------------------------------------------------------------------------
License
    This file is part of CFDEMcoupling academic.

    CFDEMcoupling academic is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    CFDEMcoupling academic is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with CFDEMcoupling academic.  If not, see <http://www.gnu.org/licenses/>.
\*---------------------------------------------------------------------------*/

#include "error.H"

#include "directedDiffusiveRelaxation.H"
#include "addToRunTimeSelectionTable.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(directedDiffusiveRelaxation, 0);

addToRunTimeSelectionTable
(
    forceModel,
    directedDiffusiveRelaxation,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
directedDiffusiveRelaxation::directedDiffusiveRelaxation
(
    const dictionary& dict,
    cfdemCloud& sm
)
:
    forceModel(dict,sm),
    propsDict_(dict.subDict(typeName + "Props")),
    interpolate_(propsDict_.lookupOrDefault<bool>("interpolation", false)),
    measureDiff_(propsDict_.lookupOrDefault<bool>("measureDiff", false)),
    recErrorFile_("recurrenceError"),
    voidfractionFieldName_(propsDict_.lookupOrDefault<word>("voidfractionFieldName","voidfraction")),
    voidfraction_(sm.mesh().lookupObject<volScalarField> (voidfractionFieldName_)),
    voidfractionRecFieldName_(propsDict_.lookupOrDefault<word>("voidfractionRecFieldName","voidfractionRec")),
    voidfractionRec_(sm.mesh().lookupObject<volScalarField> (voidfractionRecFieldName_)),
    critVoidfraction_(propsDict_.lookupOrDefault<scalar>("critVoidfraction", 1.0)),
    D0_(readScalar(propsDict_.lookup("D0"))),
    maxDisplacement_(propsDict_.lookupOrDefault<scalar>("maxDisplacement", -1.0)),
    correctedField_
    (   IOobject
        (
            "correctedField",
            sm.mesh().time().timeName(),
            sm.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        voidfraction_
    ),
    relaxStream_
    (   IOobject
        (
            "relaxStream",
            sm.mesh().time().timeName(),
            sm.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        sm.mesh(),
        dimensionedVector("zero",dimensionSet(0,1,-1,0,0), vector::zero)
    ),
    DField_
    (   IOobject
        (
            "DField",
            sm.mesh().time().timeName(),
            sm.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        sm.mesh(),
        dimensionedScalar("zero", dimensionSet(0,2,-1,0,0,0,0), 0.0),
        "zeroGradient"
    ),
    dtCFD_(voidfraction_.mesh().time().deltaTValue()),
    dtDEM_(particleCloud_.dataExchangeM().DEMts()),
    timeFac_(1.0),
    ignoreReg_(propsDict_.lookupOrDefault<bool>("ignoreRegion",false)),
    ignoreDirection_(propsDict_.lookupOrDefault<vector>("ignoreDirection",vector::zero)),
    ignorePoint_(propsDict_.lookupOrDefault<vector>("ignorePoint",vector::zero)),
    relaxForT_(propsDict_.lookupOrDefault<scalar>("relaxForT", -1.0))
{
    if(ignoreReg_)
    {
        if(mag(ignoreDirection_) < SMALL)
        {
            FatalError <<"ignoreDirection vector has very little norm, please use larger one\n" << abort(FatalError);
        }
        Info << "directedDiffusiveRelaxation: ignoring fluctuations below plane specified by point " <<
        ignorePoint_ << " and normal vector " << ignoreDirection_ << endl;
    }

    if (dtDEM_ > dtCFD_) timeFac_ = dtCFD_ / dtDEM_;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

directedDiffusiveRelaxation::~directedDiffusiveRelaxation()
{
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void directedDiffusiveRelaxation::setForce() const
{
    if (relaxForT_ >= 0.0 && particleCloud_.mesh().foundObject<IOdictionary>("lastJumpTime"))
    {
        const IOdictionary& lastJumpTimeDict(particleCloud_.mesh().lookupObject<IOdictionary>("lastJumpTime"));
        scalar lastJumpTime = lastJumpTimeDict.lookupOrDefault("lastJumpTime",-1e5);
        scalar currTime = particleCloud_.mesh().time().timeOutputValue();
        if (currTime - lastJumpTime > relaxForT_) return;
    }
  
    relax(D0_);

    relaxStream_ = DField_ * fvc::grad(correctedField_ - voidfractionRec_);

   // explicit evaluation?
   // relaxStream_ = DField_ * fvc::grad(voidfraction_ - voidfractionRec_);

    vector position(0,0,0);
    scalar voidfraction(0.0);
    vector flucU(0,0,0);
    label cellI=0;

    interpolationCellPoint<scalar> voidfractionInterpolator_(voidfraction_);
    interpolationCellPoint<vector> relaxStreamInterpolator_(relaxStream_);

    scalar maxVel = maxDisplacement_ / dtDEM_;

    for(int index = 0;index <  particleCloud_.numberOfParticles(); ++index)
    {
            cellI = particleCloud_.cellIDs()[index][0];
            if (cellI > -1) // particle found
            {
                if(ignoreReg_)
                {
                    // check if cell center is below "ignore plane"
                    vector cellCenter = particleCloud_.mesh().C()[cellI];
                    scalar aboveBelow = ignoreDirection_ & (ignorePoint_ - cellCenter);
                    if(aboveBelow > 0) continue;
                }

     //           if(voidfractionRec_[cellI] < critVoidfraction_)
                {
                    if( interpolate_ )
                    {
                        position = particleCloud_.position(index);
                        voidfraction = voidfractionInterpolator_.interpolate(position,cellI);
                        flucU = relaxStreamInterpolator_.interpolate(position,cellI);
                    }
                    else
                    {
                        voidfraction = voidfraction_[cellI];
                        flucU = relaxStream_[cellI];
                    }

                    if (voidfraction > 1.0-SMALL) voidfraction = 1.0 - SMALL;
                    flucU /= (1-voidfraction);
                    flucU *= timeFac_;

                    scalar magFlucU = mag(flucU);
                    if (maxVel > 0.0 && magFlucU > maxVel)
                    {
                        flucU *= maxVel / magFlucU;
                    }
                    
                    // write particle based data to global array
                    for(int j=0;j<3;j++)
                    {
                        particleCloud_.particleFlucVels()[index][j] += flucU[j];
                    }
                }
            }
    }
    
    if (measureDiff_)
    {
        dimensionedScalar diff( fvc::domainIntegrate( sqr( voidfraction_ - voidfractionRec_ ) ) );
        scalar t = particleCloud_.mesh().time().timeOutputValue(); 
        recErrorFile_ << t << "\t" << diff.value() << endl;
    }
}

void directedDiffusiveRelaxation::relax(scalar D0) const
{
    correctedField_.oldTime() = voidfraction_;
    forAll(DField_, cellI)
    {
        if(voidfraction_[cellI] > 0.99) DField_[cellI] = 0.0;
        else DField_[cellI] = D0;
    }
    
    solve
    (
        fvm::ddt(correctedField_)
       -fvm::laplacian(DField_, correctedField_)
       ==
       -fvc::laplacian(DField_, voidfractionRec_)
    );
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
