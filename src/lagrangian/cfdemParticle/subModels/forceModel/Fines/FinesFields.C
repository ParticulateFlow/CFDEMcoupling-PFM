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
    clogKin_(propsDict_.lookupOrDefault<bool>("kineticClogging",false)),
    clogStick_(propsDict_.lookupOrDefault<bool>("stickyClogging",false)),
    movingBed_(propsDict_.lookupOrDefault<bool>("movingBed",true)),
    useDepositionLength_(false),
    fluxFieldName_(propsDict_.lookupOrDefault<word>("fluxFieldName","phi")),
    phi_(sm.mesh().lookupObject<surfaceScalarField> (fluxFieldName_)),
    velFieldName_(propsDict_.lookupOrDefault<word>("velFieldName","U")),
    U_(sm.mesh().lookupObject<volVectorField> (velFieldName_)),
    voidfractionFieldName_(propsDict_.lookupOrDefault<word>("voidfractionFieldName","voidfraction")),
    voidfraction_(sm.mesh().lookupObjectRef<volScalarField> (voidfractionFieldName_)),
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
    alphaMax_(propsDict_.lookupOrDefault<scalar>("alphaMax",0.95)),
    alphaMinClog_(propsDict_.lookupOrDefault<scalar>("alphaMinClog",0.3)),
    critVoidfraction_(propsDict_.lookupOrDefault<scalar>("critVoidfraction", 0.05)),
    deltaT_(voidfraction_.mesh().time().deltaTValue()),
    depositionLength_(0.0),
    exponent_(-1.33),
    nCrit_(0.0),
    poresizeWidth_(0.0),
    prefactor_(10.5),
    ratioHydraulicPore_(1.5),
    tauDeposition_(0.0),
    tauRelease_(readScalar(propsDict_.lookup ("tauRelease"))),
    uBindHigh_(propsDict_.lookupOrDefault<scalar>("uBindHigh",0.0)),
    uBindLow_(propsDict_.lookupOrDefault<scalar>("uBindLow",0.0)),
    uMin_(0.001)
{
    Sds_.write();

    if (propsDict_.found("prefactor"))
    {
        prefactor_=readScalar(propsDict_.lookup ("prefactor"));
    }
    if (propsDict_.found("exponent"))
    {
        exponent_=readScalar(propsDict_.lookup ("exponent"));
    }
    if (propsDict_.found("dFine"))
    {
        dFine_.value()=readScalar(propsDict_.lookup ("dFine"));
    }
    else
    {
        FatalError <<"Please specify dFine.\n" << abort(FatalError);
    }
    if (propsDict_.found("diffCoeff"))
    {
        diffCoeff_.value()=readScalar(propsDict_.lookup ("diffCoeff"));
    }
    if (propsDict_.found("rhoFine"))
    {
        rhoFine_.value()=readScalar(propsDict_.lookup ("rhoFine"));
    }
    else
    {
        FatalError <<"Please specify rhoFine.\n" << abort(FatalError);
    }
    if (propsDict_.found("nuAve"))
    {
        nuAve_.value()=readScalar(propsDict_.lookup ("nuAve"));
    }
    if (propsDict_.found("alphaDynMax"))
    {
        alphaDynMax_=readScalar(propsDict_.lookup ("alphaDynMax"));
    }

    if (clogKin_)
    {
        nCrit_ = readScalar(propsDict_.lookup ("nCrit"));
        poresizeWidth_ = readScalar(propsDict_.lookup ("poresizeWidth"));
    }

    if (clogStick_)
    {
        if (uBindHigh_ - uBindLow_ < SMALL)
        {
            FatalError <<"No reasonable values for critical binding velocities.\n" << abort(FatalError);
        }

        uReconstructed_.set
        (
            new volVectorField
            (
                IOobject
                (
                    "uReconstructed",
                    sm.mesh().time().timeName(),
                    sm.mesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                U_//sm.mesh()
            )
        );

    }

    if (propsDict_.found("depositionLength") && propsDict_.found("tauDeposition"))
    {
        FatalError <<"You cannot specify both a deposition length and time.\n" << abort(FatalError);
    }
    else if (propsDict_.found("depositionLength"))
    {
        depositionLength_ = readScalar(propsDict_.lookup ("depositionLength"));
        useDepositionLength_ = true;
    }
    else if (propsDict_.found("tauDeposition"))
    {
        tauDeposition_ = readScalar(propsDict_.lookup ("tauDeposition"));
        useDepositionLength_ = false;
    }
    else
    {
        FatalError <<"Specify either a deposition length or time.\n" << abort(FatalError);
    }

    if (verbose_)
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
        if (clogStick_)
        {
            uReconstructed_().writeOpt() = IOobject::AUTO_WRITE;
            uReconstructed_().write();
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

FinesFields::~FinesFields()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void FinesFields::update()
{
    if (verbose_)  Info << "FinesFields: Updating alphaP.\n" << endl;
    updateAlphaP();
    if (verbose_)  Info << "FinesFields: Updating alphaG.\n" << endl;
    updateAlphaG();
    if (clogStick_)
    {
        if (verbose_)  Info << "FinesFields: Updating uReconstructed.\n" << endl;
        updateUReconstructed();
    }
    if (verbose_)  Info << "FinesFields: Calculating source terms.\n" << endl;
    calcSource();
    if (verbose_)  Info << "FinesFields: Updating dSauter.\n" << endl;
    updateDSauter();
    if (verbose_)  Info << "FinesFields: Updating dHydMix.\n" << endl;
    updateDHydMix();
    if (verbose_)  Info << "FinesFields: Updating Froude.\n" << endl;
    updateFroude();
    if (verbose_)  Info << "FinesFields: Updating FanningCoeff.\n" << endl;
    updateFanningCoeff();
    if (verbose_)  Info << "FinesFields: Updating DragCoeff.\n" << endl;
    updateDragCoeff();
    if (verbose_)  Info << "FinesFields: Updating uDyn.\n" << endl;
    updateUDyn();
    if (verbose_)  Info << "FinesFields: Integrating alphas.\n" << endl;
    integrateFields();
    if (verbose_)  Info << "FinesFields: Update finished.\n" << endl;
}


void FinesFields::calcSource()
{
    if (!clogKin_ & !clogStick_) return;
    Sds_.primitiveFieldRef() = 0;
    deltaAlpha_.primitiveFieldRef() = 0.0;
    scalar fKin = 0.0;
    scalar fStick = 0.0;
    scalar critpore = 0.0;
    scalar dmean = 0.0;
    scalar d1 = 0.0;
    scalar d2 = 0.0;
    scalar magU = 0.0;
    scalar tauDeposition = 0.0;

    forAll(Sds_,cellI)
    {
        fKin = 0.0;
        fStick = 0.0;
        if (clogKin_ && alphaP_[cellI] > alphaMinClog_) // no kinetic cloggig in dilute regions
        {
            // calculate everything in units auf dSauter
            critpore = nCrit_*dFine_.value()/dSauter_[cellI];
            // pore size from hydraulic radius
            dmean = 2 * (1 - alphaP_[cellI]) / ( (1 + poresizeWidth_*poresizeWidth_/3) * 3 * alphaP_[cellI] );
            // Sweeney and Martin, Acta Materialia 51 (2003): ratio of hydraulic to pore throat radius
            dmean /= ratioHydraulicPore_;
            d1 = dmean * (1 - poresizeWidth_);
            d2 = dmean * (1 + poresizeWidth_);

            fKin = (critpore*critpore*critpore - d1 * d1 * d1) / (d2 * d2 * d2 - d1 * d1 * d1);
            if (fKin < 0) fKin = 0.0;
            else if (fKin > 1.0) fKin = 1.0;
        }

        if (clogStick_)
        {
            magU = mag(uReconstructed_()[cellI]); // use U reconstructed from phi to suppress oscillations at interfaces
           // fStick = 1.0 / ( 1.0 + magU/uBind_) * alphaP_[cellI] / 0.65;
            if (magU < uBindLow_) fStick = 1.0;
            else if (magU > uBindHigh_) fStick = 0.0;
            else fStick = 1.0 - (magU - uBindLow_) / (uBindHigh_ - uBindLow_);
      //      if (fStick > 1.0) fStick = 1.0;
            fStick *= alphaP_[cellI] / 0.65;
        }

        // at this point, voidfraction is still calculated from the true particle sizes
        deltaAlpha_[cellI] = max(fKin,fStick) * (alphaMax_ - alphaP_[cellI]) - alphaSt_[cellI];
        // too much volume occupied: release it
        if (deltaAlpha_[cellI] < 0.0)
        {
            Sds_[cellI] = deltaAlpha_[cellI] / tauRelease_;
        }
        // volume to occupy available
        else
        {
            if (useDepositionLength_)
            {
                tauDeposition = depositionLength_ / max(mag(U_[cellI]),uMin_);
            }
            else
            {
                tauDeposition = tauDeposition_;
            }
            if (tauDeposition < 4.0*deltaT_) tauDeposition = 4.0*deltaT_; // do not deposit all dyn hold-up within less than 3 time steps
            Sds_[cellI] = alphaDyn_[cellI] / tauDeposition;
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
        ==
        Sds_
    );
    if (movingBed_) alphaStEqn += fvm::div(phiSt,alphaSt_);

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

    if (smoothing_)
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


void FinesFields::updateAlphaG() // called after updateAlphaP() - correct voidfraction by removing space occupied by fines
{
    alphaG_ = max(voidfraction_ - alphaSt_ - alphaDyn_, critVoidfraction_);
    voidfraction_ = alphaG_;
}


void FinesFields::updateAlphaP() // called first in the update cycle - voidfraction_ is current with DEM data
{
    alphaP_ = 1.0 - voidfraction_ + SMALL;
}


void FinesFields::updateDHydMix()
{
    forAll(dHydMix_,cellI)
    {
        scalar aPSt =  alphaP_[cellI] + alphaSt_[cellI];
        if (aPSt < SMALL || aPSt > 1 - SMALL)
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
        if (Ref1 <= SMALL)
            Cd = 24.0 / SMALL;
        else if (Ref1 <= 1.0)
            Cd = 24.0 / Ref1;
        else if (Ref1 <= 1000)
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
            if (Ref1 <= SMALL)
                Cd = 24.0 / SMALL;
            else if (Ref1 <= 1.0)
                Cd = 24.0 / Ref1;
            else if (Ref1 <= 1000)
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
        if (aSt < SMALL)
            dSauterMix_[cellI] = dSauter_[cellI];
        else if (aP < SMALL)
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
        if (muDyn > mU && muDyn > SMALL)
        {
            uDyn_[cellI] *= mU / muDyn;
        }
    }

    forAll(uDyn_.boundaryField(), patchI)
        forAll(uDyn_.boundaryField()[patchI], faceI)
        {
            scalar mU(mag(U_.boundaryField()[patchI][faceI]));
            scalar muDyn(mag(uDyn_.boundaryField()[patchI][faceI]));
            if (muDyn > mU && muDyn > SMALL)
            {
                uDyn_.boundaryFieldRef()[patchI][faceI] *= mU / muDyn;
            }
        }
}

void FinesFields::updateUReconstructed()
{
    if (phi_.dimensions() == dimensionSet(1, 0, -1, 0, 0)) // compressible
    {
        uReconstructed_() = fvc::reconstruct(phi_) / (rhoG_ * voidfraction_);
    }
    else
    {
        uReconstructed_() = fvc::reconstruct(phi_) / voidfraction_;
    }
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
