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
    // this is probably really bad
    voidfraction_(const_cast<volScalarField&>(sm.mesh().lookupObject<volScalarField> (voidfractionFieldName_))),
    UsFieldName_(propsDict_.lookup("granVelFieldName")),
    UsField_(sm.mesh().lookupObject<volVectorField> (UsFieldName_)),
    dSauter_(sm.mesh().lookupObject<volScalarField> ("dSauter")),
    dSauterMix_
    (   IOobject
        (
            "dSauterMix",
            sm.mesh().time().timeName(),
            sm.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        sm.mesh(),
	dimensionedScalar("zero", dimensionSet(0,1,0,0,0), 0)
    ),
    dHydMix_
    (   IOobject
        (
            "dHydMix",
            sm.mesh().time().timeName(),
            sm.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        sm.mesh(),
	dimensionedScalar("zero", dimensionSet(0,1,0,0,0), 0)
    ),
    alphaP_
    (   IOobject
        (
            "alphaP",
            sm.mesh().time().timeName(),
            sm.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        sm.mesh(),
	dimensionedScalar("zero", dimensionSet(0,0,0,0,0), 0)
    ),
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
    Froude_
    (   IOobject
        (
            "Froude",
            sm.mesh().time().timeName(),
            sm.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        sm.mesh(),
	dimensionedScalar("zero", dimensionSet(0,0,0,0,0), 0)
    ),
    estimatedVelFrac_(readScalar(propsDict_.lookup ("estimatedVelFrac"))),
    dFine_("dFine", dimensionSet(0,1,0,0,0), 0.0),
    depRate_(readScalar(propsDict_.lookup ("depRate"))),
    rhoDyn_(readScalar(propsDict_.lookup ("rhoDyn"))),
    nCrit_(readScalar(propsDict_.lookup ("nCrit"))),
    alphaMax_(readScalar(propsDict_.lookup ("alphaMax"))),
    critVoidfraction_(readScalar(propsDict_.lookup ("critVoidfraction"))),
    g("g",dimensionSet(0,1,-2,0,0),9.81)
{
    dFine_.value()=readScalar(propsDict_.lookup ("dFine"));
    Sds_.write();
    // init force sub model
 //   setForceSubModels(propsDict_);
    // define switches which can be read from dict


    

}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

FinesFields::~FinesFields()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void FinesFields::update()
{
    // the sequence is important: for the source terms, the Sauter mean diameter of the large particles
    // is needed, for the force calculation, the static hold-up's contribution needs to be included
    alphaP_ = 1.0 - voidfraction_;
    calcSource();
    updateDSauter();
    dHydMix_ = 2*(1 - alphaP_ - alphaSt_) / (3*(alphaP_ + alphaSt_) ) * dSauterMix_;
    Froude_ = alphaDyn_*mag(uDyn_ - UsField_) / Foam::sqrt(dHydMix_*g_); 
    updateUDyn();
    integrateFields();
    updateVoidfraction();
}

void FinesFields::updateDSauter() 
{
  dSauterMix_ = (alphaP_ + alphaSt_) / ( alphaP_/dSauter_ + alphaSt_/dFine_ );
  
  // manipulate voidfraction field
  
}

void FinesFields::calcSource() 
{
    Sds_=0;
    label cellI=0;
    scalar f(0.0);
    scalar critpore(0.0);
    scalar deltaAlpha(0.0);
    
    forAll(Sds_,cellI)
    {
        critpore = nCrit_*dFine_.value()/dSauter_[cellI];
        // mean pore diameter \approx 0.2 mean particle diameter, width +- 0.075 particle diameters
        f = (critpore*critpore*critpore - 0.0019531) / (0.020797 - 0.0019531);
	if (f<0)
	{
	    f=0.0;    
	}
	else if (f>1.0)
	{
	    f=1.0;
        }
	
	// at this point, voidfraction is still calculated from the true particle sizes
	deltaAlpha = f * (alphaMax_ - alphaP_[cellI]) - alphaSt_[cellI];
	
	if (deltaAlpha < 0)
	{
	    Sds_[cellI] = deltaAlpha;
	}
	else if (deltaAlpha > alphaDyn_[cellI])
	{
	    Sds_[cellI] = alphaDyn_[cellI];
	}
	else
	{
	    Sds_[cellI] = depRate_ * deltaAlpha;
	}
    }
  
}

void FinesFields::updateUDyn() 
{
  
  
  // more implementation to follow

  
}

void FinesFields::integrateFields() 
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

void FinesFields::updateVoidfraction() 
{
  forAll(voidfraction_,cellI)
  {
      if(voidfraction_[cellI] - alphaSt_[cellI] > critVoidfraction_)
      {
	  voidfraction_[cellI] -= alphaSt_[cellI];
      }
      else
      {
	  voidfraction_[cellI] = critVoidfraction_;
      }
  }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
