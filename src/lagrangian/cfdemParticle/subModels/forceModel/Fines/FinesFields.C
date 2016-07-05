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
    velFieldName_(propsDict_.lookupOrDefault<word>("velFieldName","U")),
    U_(sm.mesh().lookupObject<volVectorField> (velFieldName_)),
    voidfractionFieldName_(propsDict_.lookupOrDefault<word>("voidfractionFieldName","voidfraction")),
    // this is probably really bad
    voidfraction_(const_cast<volScalarField&>(sm.mesh().lookupObject<volScalarField> (voidfractionFieldName_))),
    UsFieldName_(propsDict_.lookupOrDefault<word>("granVelFieldName","Us")),
    UsField_(sm.mesh().lookupObject<volVectorField> (UsFieldName_)),
    pFieldName_(propsDict_.lookupOrDefault<word>("pFieldName","p")),
    p_(sm.mesh().lookupObject<volScalarField> (pFieldName_)),
    rhoGFieldName_(propsDict_.lookupOrDefault<word>("rhoGFieldName","rho")),
    rhoG_(sm.mesh().lookupObject<volScalarField> (rhoGFieldName_)),
    dSauter_(sm.mesh().lookupObject<volScalarField> ("dSauter")),   
    alphaG_
    (   IOobject
        (
            "alphaG",
            sm.mesh().time().timeName(),
            sm.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        sm.mesh(),
        dimensionedScalar("zero", dimensionSet(0,0,0,0,0), 0)
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
    DragCoeff_
    (   IOobject
        (
            "DragCoeff",
            sm.mesh().time().timeName(),
            sm.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        sm.mesh(),
        dimensionedScalar("zero", dimensionSet(1,-3,-1,0,0), 0)
    ),
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
    FanningCoeff_
    (   IOobject
        (
            "FanningCoeff",
            sm.mesh().time().timeName(),
            sm.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        sm.mesh(),
        dimensionedScalar("zero", dimensionSet(1,-3,-1,0,0), 0)
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
    dFine_("dFine", dimensionSet(0,1,0,0,0), 0.0),
    nuAve_("nuAve",dimensionSet(0,2,-1,0,0),1.6e-5),
    rhoFine_("rhoFine",dimensionSet(1,-3,0,0,0),0.0),
    g_("g",dimensionSet(0,1,-2,0,0),vector(0,0,-9.81)),
    alphaMax_(readScalar(propsDict_.lookup ("alphaMax"))),
    critVoidfraction_(readScalar(propsDict_.lookup ("critVoidfraction"))),
    depRate_(readScalar(propsDict_.lookup ("depRate"))),
    exponent_(-1.33),
    nCrit_(readScalar(propsDict_.lookup ("nCrit"))),
    prefactor_(14.98) 
{
    dFine_.value()=readScalar(propsDict_.lookup ("dFine"));
    Sds_.write();

    if (propsDict_.found("prefactor"))
        prefactor_=readScalar(propsDict_.lookup ("prefactor"));
    if (propsDict_.found("exponent"))
        exponent_=readScalar(propsDict_.lookup ("exponent"));
    if (propsDict_.found("dFine"))
        dFine_.value()=readScalar(propsDict_.lookup ("dFine"));
    else
        FatalError <<"Please specify dFine.\n" << abort(FatalError);
    if (propsDict_.found("rhoFine"))
        rhoFine_.value()=readScalar(propsDict_.lookup ("rhoFine"));
    else
        FatalError <<"Please specify rhoFine.\n" << abort(FatalError);  
    if (propsDict_.found("nuAve"))
        nuAve_.value()=readScalar(propsDict_.lookup ("nuAve"));
    

}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

FinesFields::~FinesFields()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void FinesFields::update()
{
    updateAlphaP();
    updateAlphaG();
    calcSource();
    updateDSauter();
    updateDHydMix();
    updateFroude();
    updateFanningCoeff();
    updateUDyn();
    integrateFields();
    // update voidfraction, probably really bad...
    voidfraction_ = alphaG_;
}


void FinesFields::calcSource() 
{
    Sds_=0;
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


void FinesFields::updateAlphaG() 
{
  forAll(alphaG_,cellI)
  {
      if(voidfraction_[cellI] - alphaSt_[cellI] > critVoidfraction_)
      {
	  alphaG_[cellI] = voidfraction_[cellI] - alphaSt_[cellI];
      }
      else
      {
	  alphaG_[cellI] = critVoidfraction_;
      }
  }
}


void FinesFields::updateAlphaP() 
{
    alphaP_ = 1.0 - voidfraction_;
}


void FinesFields::updateDHydMix() 
{
    dHydMix_ = 2*(1 - alphaP_ - alphaSt_) / (3*(alphaP_ + alphaSt_) ) * dSauterMix_;
}


void FinesFields::updateDragCoeff()
{
    volScalarField beta = 0.75 * rhoG_ * alphaDyn_ / dFine_ * mag(U_ - uDyn_) * Foam::pow(alphaG_,-4.65);
    volScalarField Ref = dFine_ * alphaG_ / nuAve_ * mag(U_ - uDyn_);
    scalar Cd(0.0);
    scalar Ref1(0.0);
    forAll(DragCoeff_,cellI)
    {
        Ref1 = Ref[cellI];
        if(Ref1 <= SMALL)
	    Cd = 24.0 / SMALL;
        else if(Ref1 <= 1.0)
	    Cd = 24.0 / Ref1;
	else if(Ref1 <= 1000)
	    Cd = 24 * (1.0 + 0.15 * Foam::pow(Ref1,0.687) ) / Ref1;
	else
	    Cd = 0.44;
	DragCoeff_ = Cd * beta[cellI];
    }
}


void FinesFields::updateDSauter() 
{
  dSauterMix_ = (alphaP_ + alphaSt_) / ( alphaP_/dSauter_ + alphaSt_/dFine_ );  
}


void FinesFields::updateFanningCoeff()
{
    FanningCoeff_ = alphaDyn_ * rhoFine_ * mag(uDyn_ - UsField_) * prefactor_ * Foam::pow(Froude_, exponent_) / (2 * dHydMix_);
}


void FinesFields::updateFroude() 
{
    Froude_ = alphaDyn_ * mag(uDyn_ - UsField_) / Foam::sqrt(dHydMix_*mag(g_));
}


void FinesFields::updateUDyn() 
{
    volVectorField num = rhoFine_ * alphaDyn_ * g_ + FanningCoeff_ * UsField_ + DragCoeff_ * U_;
    if (p_.dimensions()==dimensionSet(0,2,-2,0,0) )
        num -= alphaDyn_ * rhoG_ * fvc::grad(p_) ;
    else
        num -= alphaDyn_ * fvc::grad(p_) ;
    volScalarField denom = FanningCoeff_ + DragCoeff_;
  
    uDyn_ = num / denom;  
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
