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

#include "FinesFields.H"
#include "addToRunTimeSelectionTable.H"
#include "averagingModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
FinesFields::FinesFields
(
    const dictionary& dict,
    cfdemCloud& sm
)
:
    particleCloud_(sm),
    propsDict_(dict.subDict(typeName + "Props")),
    velFieldName_(propsDict_.lookup("velFieldName")),
    U_(sm.mesh().lookupObject<volVectorField> (velFieldName_)),
    voidfractionFieldName_(propsDict_.lookup("voidfractionFieldName")),
    voidfraction_(sm.mesh().lookupObject<volScalarField> (voidfractionFieldName_)),
    UsFieldName_(propsDict_.lookup("granVelFieldName")),
    UsField_(sm.mesh().lookupObject<volVectorField> (UsFieldName_)),
    alphaSt_
    (   IOobject
        (
            "alphaSt",
            sm.mesh().time().timeName(),
            sm.mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        sm.mesh()
    ),
    alphaDyn_
    (   IOobject
        (
            "alphaDyn",
            sm.mesh().time().timeName(),
            sm.mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        sm.mesh()
    ),
    uDyn_
    (   IOobject
        (
            "uDyn",
            sm.mesh().time().timeName(),
            sm.mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        sm.mesh()
    ),
    Sds_
    (   IOobject
        (
            "Sds",
            sm.mesh().time().timeName(),
            sm.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        sm.mesh(),
	dimensionedScalar("zero", dimensionSet(0,0,0,0,0), 0)
    ),
    estimatedVelFrac_(0.2)
{
    if (propsDict_.found("estimatedVelFrac"))
    {
        estimatedVelFrac_=readScalar(propsDict_.lookup ("estimatedVelFrac"));
    }
    Sds_.write();
    // init force sub model
 //   setForceSubModels(propsDict_);
    // define switches which can be read from dict


    

}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

FinesFields::~FinesFields()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //



void FinesFields::updateUDyn() 
{
  volVectorField U0 = estimatedVelFrac_*U_ + (1-estimatedVelFrac_)*UsField_;
  
  // more implementation to follow
  uDyn_ = U0;
  
}

void FinesFields::calcSource() 
{
    Sds_=0;
  
  
}

void FinesFields::updateFields() 
{
  
    surfaceScalarField phiSt(linearInterpolate(UsField_) & particleCloud_.mesh().Sf());
    surfaceScalarField phiDyn(linearInterpolate(uDyn_) & particleCloud_.mesh().Sf());
  
    fvScalarMatrix alphaStEqn
    (
          fvm::ddt(alphaSt_)
	+ fvm::div(phiSt,alphaSt_)
	==
	Sds_
    );
    fvScalarMatrix alphaDynEqn
    (
        fvm::ddt(alphaDyn_)
	+ fvm::div(phiDyn,alphaDyn_)
	==
	-Sds_
    );
    alphaStEqn.solve();
    alphaDynEqn.solve();
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
