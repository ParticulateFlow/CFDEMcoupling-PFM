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
    D1_(readScalar(propsDict_.lookup("D1"))),
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
    dt_(particleCloud_.dataExchangeM().DEMts()),
    ignoreReg_(propsDict_.lookupOrDefault<bool>("ignoreRegion",false)),
    ignoreDirection_(propsDict_.lookupOrDefault<vector>("ignoreDirection",vector::zero)),
    ignorePoint_(propsDict_.lookupOrDefault<vector>("ignorePoint",vector::zero)),
    vfluc_(NULL)
{
    allocateMyArrays();
    
    if(ignoreReg_)
    {
        if(mag(ignoreDirection_) < SMALL)
        {
            FatalError <<"ignoreDirection vector has very little norm, please use larger one\n" << abort(FatalError);
        }
        Info << "directedDiffusiveRelaxation: ignoring fluctuations below plane specified by point " <<
        ignorePoint_ << " and normal vector " << ignoreDirection_ << endl;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

directedDiffusiveRelaxation::~directedDiffusiveRelaxation()
{
    delete vfluc_;
}


// * * * * * * * * * * * * * * * private Member Functions  * * * * * * * * * * * * * //
void directedDiffusiveRelaxation::allocateMyArrays() const
{
    // get memory for 2d arrays
    double initVal=0.0;
    particleCloud_.dataExchangeM().allocateArray(vfluc_,initVal,3); 
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void directedDiffusiveRelaxation::setForce() const
{

    relax(D0_,D1_);
    
    volVectorField relaxStream = DField_ * fvc::grad(correctedField_ - voidfractionRec_);
    
  // volVectorField relaxStream = DField_ * fvc::grad(voidfraction_ - voidfractionRec_);
    
     // realloc the arrays
    allocateMyArrays();
    
    vector position(0,0,0);
    scalar voidfraction(0.0);
    vector flucU(0,0,0);
    label cellI=0;
   
    interpolationCellPoint<scalar> voidfractionInterpolator_(voidfraction_);
    interpolationCellPoint<vector> relaxStreamInterpolator_(relaxStream);
    
    scalar dtDEM = particleCloud_.dataExchangeM().DEMts();
    scalar dtCFD = voidfraction_.mesh().time().deltaTValue();
    
    // if DEM time step > CFD time step, scale velocity down
    scalar timeFac = 1.0;
    if (dtDEM > dtCFD) timeFac = dtCFD / dtDEM;

    
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
                        flucU = relaxStream[cellI];
                    }
                    
                    if (voidfraction > 1.0-SMALL) voidfraction = 1.0 - SMALL;
                    flucU /= (1-voidfraction);
		    flucU *= timeFac;
                    // write particle based data to global array 
                    for(int i = 0; i < 3; i++)
                    {
                        vfluc_[index][i]=flucU[i];
                    }
                }
            }
    }
    
    particleCloud_.dataExchangeM().giveData("vfluc","vector-atom", vfluc_);
    
    if (measureDiff_)
    {
        dimensionedScalar diff( fvc::domainIntegrate( sqr( voidfraction_ - voidfractionRec_ ) ) );
        scalar t = particleCloud_.mesh().time().timeOutputValue(); 
        recErrorFile_ << t << "\t" << diff.value() << endl;
    }
}

void directedDiffusiveRelaxation::relax(scalar D0, scalar D1) const
{
    correctedField_.oldTime() = voidfraction_;
    forAll(DField_, cellI)
    {
        if(voidfraction_[cellI] > 0.99) DField_[cellI] = 0.0;
        else DField_[cellI] = D0 + D1*(voidfractionRec_[cellI] - voidfraction_[cellI]);
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
