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
#include "smoothingModel.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

//defineTypeNameAndDebug(FinesFields, 0);

//defineRunTimeSelectionTable(FinesFields, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
FinesFields::FinesFields
(
    const dictionary& dict,
    cfdemCloud& sm
)
:
    particleCloud_(sm),
    propsDict_(dict.subDict("FinesFieldsProps")),
    smoothing_(propsDict_.lookupOrDefault<bool>("smoothing",false)),
    verbose_(propsDict_.lookupOrDefault<bool>("verbose",false)),
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
    deltaAlpha_
    (   IOobject
        (
            "deltaAlpha",
            sm.mesh().time().timeName(),
            sm.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        sm.mesh(),
        dimensionedScalar("zero", dimensionSet(0,0,0,0,0), 0.0)
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
        dimensionedScalar("zero", dimensionSet(0,1,0,0,0), 0),
        "zeroGradient"
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
        dimensionedScalar("zero", dimensionSet(0,1,0,0,0), 0),
        "zeroGradient"
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
        dimensionedScalar("zero", dimensionSet(0,0,-1,0,0), 0)
    ),
    massFluxDyn_
    (   IOobject
        (
            "massFluxDyn",
            sm.mesh().time().timeName(),
            sm.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        sm.mesh(),
        dimensionedScalar("zero", dimensionSet(1,0,-1,0,0), 0)
        //dimensionedVector("zero", dimensionSet(1,-2,-1,0,0), vector::zero)
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
    diffCoeff_("diffCoeff",dimensionSet(0,2,-1,0,0,0,0),0.0),
    nuAve_("nuAve",dimensionSet(0,2,-1,0,0),1.6e-5),
    rhoFine_("rhoFine",dimensionSet(1,-3,0,0,0),0.0),
    g_("g",dimensionSet(0,1,-2,0,0),vector(0,0,-9.81)),
    alphaDynMax_(0.1),
    alphaMax_(readScalar(propsDict_.lookup ("alphaMax"))),
    critVoidfraction_(readScalar(propsDict_.lookup ("critVoidfraction"))),
    depRate_(readScalar(propsDict_.lookup ("depRate"))),
    exponent_(-1.33),
    nCrit_(readScalar(propsDict_.lookup ("nCrit"))),
    poresizeWidth_(readScalar(propsDict_.lookup ("poresizeWidth"))),
    prefactor_(10.5),
    ratioHydraulicPore_(1.5)
{
    Sds_.write();

    if (propsDict_.found("prefactor"))
        prefactor_=readScalar(propsDict_.lookup ("prefactor"));
    if (propsDict_.found("exponent"))
        exponent_=readScalar(propsDict_.lookup ("exponent"));
    if (propsDict_.found("dFine"))
        dFine_.value()=readScalar(propsDict_.lookup ("dFine"));
    else
        FatalError <<"Please specify dFine.\n" << abort(FatalError);
    if (propsDict_.found("diffCoeff"))
        diffCoeff_.value()=readScalar(propsDict_.lookup ("diffCoeff"));
    if (propsDict_.found("rhoFine"))
        rhoFine_.value()=readScalar(propsDict_.lookup ("rhoFine"));
    else
        FatalError <<"Please specify rhoFine.\n" << abort(FatalError);
    if (propsDict_.found("nuAve"))
        nuAve_.value()=readScalar(propsDict_.lookup ("nuAve"));
    if (propsDict_.found("alphaDynMax"))
        alphaDynMax_=readScalar(propsDict_.lookup ("alphaDynMax"));

    if(verbose_)
    {
      alphaG_.writeOpt() = IOobject::AUTO_WRITE;
      alphaG_.write();
      alphaP_.writeOpt() = IOobject::AUTO_WRITE;
      alphaP_.write();
      deltaAlpha_.writeOpt() = IOobject::AUTO_WRITE;
      deltaAlpha_.write();
      dHydMix_.writeOpt() = IOobject::AUTO_WRITE;
      dHydMix_.write();
      DragCoeff_.writeOpt() = IOobject::AUTO_WRITE;
      DragCoeff_.write();
      dSauterMix_.writeOpt() = IOobject::AUTO_WRITE;
      dSauterMix_.write();
      FanningCoeff_.writeOpt() = IOobject::AUTO_WRITE;
      FanningCoeff_.write();
      Froude_.writeOpt() = IOobject::AUTO_WRITE;
      Froude_.write();
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

FinesFields::~FinesFields()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void FinesFields::update()
{
    if(verbose_)  Info << "FinesFields: Updating alphaP.\n" << endl;
    updateAlphaP();
    if(verbose_)  Info << "FinesFields: Updating alphaG.\n" << endl;
    updateAlphaG();
    if(verbose_)  Info << "FinesFields: Calculating source terms.\n" << endl;
    calcSource();
    if(verbose_)  Info << "FinesFields: Updating dSauter.\n" << endl;
    updateDSauter();
    if(verbose_)  Info << "FinesFields: Updating dHydMix.\n" << endl;
    updateDHydMix();
    if(verbose_)  Info << "FinesFields: Updating Froude.\n" << endl;
    updateFroude();
    if(verbose_)  Info << "FinesFields: Updating FanningCoeff.\n" << endl;
    updateFanningCoeff();
    if(verbose_)  Info << "FinesFields: Updating DragCoeff.\n" << endl;
    updateDragCoeff();
    if(verbose_)  Info << "FinesFields: Updating uDyn.\n" << endl;
    updateUDyn();
    if(verbose_)  Info << "FinesFields: Integrating alphas.\n" << endl;
    integrateFields();
    if(verbose_)  Info << "FinesFields: Update finished.\n" << endl;
}


void FinesFields::calcSource()
{
    Sds_.primitiveFieldRef() = 0;
    deltaAlpha_.primitiveFieldRef() = 0.0;
    scalar f(0.0);
    scalar critpore(0.0);
    scalar dmean(0.0);
    scalar d1(0.0);
    scalar d2(0.0);

    forAll(Sds_,cellI)
    {
        // calculate everything in units auf dSauter
        critpore = nCrit_*dFine_.value()/dSauter_[cellI];
        // pore size from hydraulic radius
        dmean = 2 * (1 - alphaP_[cellI]) / ( (1 + poresizeWidth_*poresizeWidth_/3) * 3 * alphaP_[cellI] );
        // Sweeney and Martin, Acta Materialia 51 (2003): ratio of hydraulic to pore throat radius
        dmean /= ratioHydraulicPore_;
        d1 = dmean * (1 - poresizeWidth_);
        d2 = dmean * (1 + poresizeWidth_);

        f = (critpore*critpore*critpore - d1 * d1 * d1) / (d2 * d2 * d2 - d1 * d1 * d1);
        if (f < 0)
        {
            f = 0.0;
        }
        else if (f > 1.0)
        {
            f = 1.0;
        }

        // at this point, voidfraction is still calculated from the true particle sizes
        deltaAlpha_[cellI] = f * (alphaMax_ - alphaP_[cellI]) - alphaSt_[cellI];
        // too much volume occupied: release it (50% per time step)
        if (deltaAlpha_[cellI] < 0.0)
        {
            Sds_[cellI] = 0.5*deltaAlpha_[cellI];
        }
        // volume too occupy available: deposit at most 80% of dyn hold up
        else if (depRate_ * deltaAlpha_[cellI] > 0.8 * alphaDyn_[cellI])
        {
            Sds_[cellI] = 0.8 * alphaDyn_[cellI];
        }
        else
        {
            Sds_[cellI] = depRate_ * deltaAlpha_[cellI];
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
        - fvm::laplacian(diffCoeff_,alphaDyn_)
        ==
        -Sds_
    );
    alphaStEqn.solve();
    alphaDynEqn.solve();

    if(smoothing_)
        particleCloud_.smoothingM().smoothen(alphaDyn_);

    // limit hold-ups, should be done more elegantly

    scalar alphaStErr(0.0);
    scalar alphaDynErr1(0.0);
    scalar alphaDynErr2(0.0);
    forAll(alphaSt_, cellI)
    {
        if (alphaSt_[cellI] < 0.0)
        {
            alphaStErr += alphaSt_[cellI] * particleCloud_.mesh().V()[cellI];
            alphaSt_[cellI] = 0.0;
        }

        if (alphaDyn_[cellI] < 0.0)
        {
            alphaDynErr1 += alphaDyn_[cellI] * particleCloud_.mesh().V()[cellI];
            alphaDyn_[cellI] = 0.0;
        }
        else if (alphaDyn_[cellI] > alphaDynMax_)
        {
            alphaDynErr2 +=  (alphaDyn_[cellI] - alphaDynMax_) * particleCloud_.mesh().V()[cellI];
            alphaDyn_[cellI] = alphaDynMax_;
        }
    }

    if (verbose_)
    {
        Sout << "[" << Pstream::myProcNo() << "] " << "amount of alphaSt added because of positivity requirement: " << -alphaStErr << endl;
        Sout << "[" << Pstream::myProcNo() << "] " << "amount of alphaDyn added because of positivity requirement: " << -alphaDynErr1 << endl;
        Sout << "[" << Pstream::myProcNo() << "] " << "amount of alphaDyn removed because of max. value: " << -alphaDynErr2 << endl;
    }

    alphaSt_.correctBoundaryConditions();
    alphaDyn_.correctBoundaryConditions();

    massFluxDyn_ = rhoFine_ * fvc::interpolate(alphaDyn_) * phiDyn;
}


void FinesFields::updateAlphaG()
{
    alphaG_ = max(voidfraction_ - alphaSt_ - alphaDyn_, critVoidfraction_);
}


void FinesFields::updateAlphaP()
{
    alphaP_ = 1.0 - voidfraction_ + SMALL;
}


void FinesFields::updateDHydMix()
{
    forAll(dHydMix_,cellI)
    {
        scalar aPSt =  alphaP_[cellI] + alphaSt_[cellI];
        if(aPSt < SMALL || aPSt > 1 - SMALL)
            dHydMix_[cellI] = SMALL;
        else
            dHydMix_[cellI] = 2*(1 - aPSt) / (3*aPSt ) * dSauterMix_[cellI];
    }
    dHydMix_.correctBoundaryConditions();
}


void FinesFields::updateDragCoeff()
{
    volScalarField beta = 0.75 * rhoG_ * alphaDyn_ / dFine_ * mag(U_ - uDyn_) * Foam::pow(alphaG_,-4.65);
    volScalarField Ref = dFine_ * alphaG_ / nuAve_ * mag(U_ - uDyn_);
    scalar Cd(0.0);
    scalar Ref1(0.0);

    // calculate drag coefficient for cells
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
        DragCoeff_[cellI] = Cd * beta[cellI];
    }

    // calculate drag coefficient for faces
    forAll(DragCoeff_.boundaryField(), patchI)
        forAll(DragCoeff_.boundaryField()[patchI], faceI)
        {
            Ref1 = Ref.boundaryField()[patchI][faceI];
            if(Ref1 <= SMALL)
                Cd = 24.0 / SMALL;
            else if(Ref1 <= 1.0)
                Cd = 24.0 / Ref1;
            else if(Ref1 <= 1000)
                Cd = 24 * (1.0 + 0.15 * Foam::pow(Ref1,0.687) ) / Ref1;
            else
                Cd = 0.44;
            DragCoeff_.boundaryFieldRef()[patchI][faceI] = Cd * beta.boundaryFieldRef()[patchI][faceI];
        }

    DragCoeff_ = max( DragCoeff_, dimensionedScalar("SMALL", dimensionSet(1,-3,-1,0,0), SMALL) );
}


void FinesFields::updateDSauter()
{
    forAll(dSauterMix_,cellI)
    {
        scalar aP = alphaP_[cellI];
        scalar aSt = alphaSt_[cellI];
        if(aSt < SMALL)
            dSauterMix_[cellI] = dSauter_[cellI];
        else if(aP < SMALL)
            dSauterMix_[cellI] = dFine_.value();
        else
            dSauterMix_[cellI] = (aP + aSt) / (aP / dSauter_[cellI] + aSt / dFine_.value() );
    }
    dSauterMix_.correctBoundaryConditions();
}


void FinesFields::updateFanningCoeff()
{
    FanningCoeff_ = alphaDyn_ * rhoFine_ * mag(uDyn_ - UsField_) * prefactor_ * Foam::pow(Froude_, exponent_) / (2 * dHydMix_);
    // FanningCoeff_ = max( FanningCoeff_, dimensionedScalar("SMALL", dimensionSet(1,-3,-1,0,0), SMALL) );
    // FanningCoeff_ = min( FanningCoeff_, dimensionedScalar("LARGE", dimensionSet(1,-3,-1,0,0), 1e5) );
}


void FinesFields::updateFroude()
{
    // seems like different authors use different conventions for the Froude number
    // Chen et al. (1994) define it in terms of a superficial velocity,
    // Shibata et a. (1991) in terms of the actual velocity
    Froude_ = max( alphaG_ * mag(uDyn_ - UsField_) / Foam::sqrt(dHydMix_*mag(g_)), 5e-3 );
  //  Froude_ = min( Froude_, 1e0 );
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

    // limit uDyn for stability reasons
    forAll(uDyn_,cellI)
    {
        scalar mU(mag(U_[cellI]));
        scalar muDyn(mag(uDyn_[cellI]));
        if(muDyn > mU && muDyn > SMALL)
        {
            uDyn_[cellI] *= mU / muDyn;
        }
    }

    forAll(uDyn_.boundaryField(), patchI)
        forAll(uDyn_.boundaryField()[patchI], faceI)
        {
            scalar mU(mag(U_.boundaryField()[patchI][faceI]));
            scalar muDyn(mag(uDyn_.boundaryField()[patchI][faceI]));
            if(muDyn > mU && muDyn > SMALL)
            {
                uDyn_.boundaryFieldRef()[patchI][faceI] *= mU / muDyn;
            }
        }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
