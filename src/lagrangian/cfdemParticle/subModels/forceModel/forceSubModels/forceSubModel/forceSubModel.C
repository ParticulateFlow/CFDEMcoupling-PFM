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
#include "forceSubModel.H"
#include "forceModel.H"
#include "mathExtra.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(forceSubModel, 0);

defineRunTimeSelectionTable(forceSubModel, dictionary);

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
forceSubModel::forceSubModel
(
    const dictionary& dict,
    cfdemCloud& sm,
    forceModel& fm
)
:
    dict_(dict),
    particleCloud_(sm),
    forceModel_(fm),
    nrDefaultSwitches_(SW_MAX),
    switchesNameList_(nrDefaultSwitches_),
    switchesList_(nrDefaultSwitches_),
    switches_(nrDefaultSwitches_),
    nu_
    (
        IOobject
        (
            "scalarViscosity",
            sm.mesh().time().timeName(),
            sm.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        sm.mesh(),
        dimensionedScalar("nu0", dimensionSet(0, 2, -1, 0, 0), 1.)
    ),
    divTau_
    (
        IOobject
        (
            "divTau",
            sm.mesh().time().timeName(),
            sm.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        sm.mesh(),
        dimensionedVector("divTau", dimensionSet(1, -2, -2, 0, 0), vector::zero)
    ),
    IBDragPerV_
    (
        IOobject
        (
            "IBDragPerV",
            sm.mesh().time().timeName(),
            sm.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        sm.mesh(),
        dimensionedVector("IBDragPerV", dimensionSet(1, -2, -2, 0, 0), vector::zero)
    ),
    densityFieldName_(dict_.lookupOrDefault<word>("densityFieldName","rho")),
    rho_(sm.mesh().lookupObject<volScalarField> (densityFieldName_))
{
    // init standard switch list
    switchesNameList_[SW_TREAT_FORCE_EXPLICIT] = "treatForceExplicit";
    switchesNameList_[SW_TREAT_FORCE_DEM] = "treatForceDEM";
    switchesNameList_[SW_IMPL_FORCE_DEM] = "implForceDEM";
    switchesNameList_[SW_VERBOSE] = "verbose";
    switchesNameList_[SW_INTERPOLATION] = "interpolation";
    switchesNameList_[SW_FILTERED_DRAG_MODEL] ="useFilteredDragModel";
    switchesNameList_[SW_PARCEL_SIZE_DEPENDENT_FILTERED_DRAG] = "useParcelSizeDependentFilteredDrag";
    switchesNameList_[SW_IMPL_FORCE_DEM_ACCUMULATED] = "implForceDEMaccumulated";
    switchesNameList_[SW_SCALAR_VISCOSITY] = "scalarViscosity";

    // sanity check of what is defined above
    if(switchesNameList_.size() != nrDefaultSwitches_)
        FatalError<< "please check the nr of switches defined in forceSubModel class." << abort(FatalError);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

forceSubModel::~forceSubModel()
{}

// * * * * * * * * * * * * * * * * Member Fct  * * * * * * * * * * * * * * * //
void forceSubModel::partToArray
(
    label index,
    vector& dragTot,
    const vector& dragEx,
    const vector& Ufluid,
    scalar Cd
) const
{
    // forces for CFD
    if(!switches_[SW_TREAT_FORCE_DEM])// !treatDEM
    {
        if(switches_[SW_TREAT_FORCE_EXPLICIT]) // treatExplicit
        {
            for(int j=0;j<3;j++)
                particleCloud_.expForces()[index][j] += dragTot[j];
        }
        else   //implicit treatment, taking explicit force contribution into account
        {
            for(int j=0;j<3;j++)
            {
                particleCloud_.impForces()[index][j] += dragTot[j] - dragEx[j]; //only consider implicit part!
                particleCloud_.expForces()[index][j] += dragEx[j];
            }
        }
    }

    // forces for DEM
    if(switches_[SW_IMPL_FORCE_DEM]) // implForceDEM
    {
        for(int j=0;j<3;j++)
            particleCloud_.fluidVels()[index][j] = Ufluid[j];

        particleCloud_.Cds()[index][0] = Cd;
    }
    else
    {
        for(int j=0;j<3;j++)
            particleCloud_.DEMForces()[index][j] += dragTot[j];
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void forceSubModel::explicitCorr
(
    vector& dragImplicit,
    vector& dragExplicit,
    scalar& dragCoefficient,
    vector& Ufluid,
    const vector& Ucell,
    vector& Us,
    const vector& UsCell,
    bool verbose,
    label index
) const
{
    dragExplicit=vector::zero;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void forceSubModel::readSwitches()
{
    Info << "\nreading switches for forceSubModel:" << myType() << endl;
    forAll(switchesNameList_,i)
    {
        if(switchesList_[i]) //check if switch is required
        {
            Info << "  looking for " << switchesNameList_[i] << " ..." << endl;
            if (dict_.found(switchesNameList_[i]))
                switches_[i] = Switch(dict_.lookup(switchesNameList_[i]));

            Info << "\t" << switchesNameList_[i] << " = " << switches_[i] << endl;
        }
    }
    Info << endl;

    if(switches_[SW_IMPL_FORCE_DEM]) // implForceDEM=true
    {
        // communicate implForceDEM to particleCloud
        particleCloud_.impDEMdrag_ = true;

        // do sanity check
        // This can work if the accumulator is used, but is explicitely applied on the CFD side
        // Sanity check is therefore not necessary here
        /*
        if(switches_[SW_TREAT_FORCE_EXPLICIT]) // treatExplicit=true
        {
            FatalError << "Please check your settings, treatExplicit together with implForceDEM does not work!."
                       << abort(FatalError);
        }
        */
    }

    if(switches_[SW_IMPL_FORCE_DEM_ACCUMULATED]) // implForceDEMaccumulated=true
    {
        // sanity check for implForceDEMaccumulated
        if(!switches_[SW_IMPL_FORCE_DEM]) //implForceDEM=false
        {
            Warning<< "please check your settings, implForceDEMaccumulated without implForceDEM does not work! (using implForceDEMaccumulated=false)" << endl;
            switches_[SW_VERBOSE] = false;
        }else
        {
            particleCloud_.impDEMdragAcc_ = true;
        }
    }

    if(switches_[SW_SCALAR_VISCOSITY]) // scalarViscosity=true
    {
        Info << "Using a constant viscosity for this force model." << endl;
        dimensionedScalar  nu0_("nu", dimensionSet(0, 2, -1, 0, 0), dict_.lookup("nu"));
        nu_ = volScalarField
        (
            IOobject
            (
                "scalarViscosity",
                particleCloud_.mesh().time().timeName(),
                particleCloud_.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            particleCloud_.mesh(),
            nu0_
        );
    }

    // look for old nomenclature
    if (dict_.found("treatExplicit") || dict_.found("treatDEM") || dict_.found("implDEM"))
        FatalError<< "You are using an old nomenclature for force model settings, please have a look at the forceSubModel doc." << abort(FatalError);

    // look for old nomenclature
    if (dict_.found("verbose"))
        Warning<< "Please make sure you use the new nomenclature for verbose force model settings, please have a look at the forceSubModel doc." << endl;

    //if (dict_.found("interpolation"))
    //    FatalError<< "Please make sure you use the new nomenclature for interpolation in force model settings, please have a look at the forceSubModel doc." << endl;
}

const volScalarField& forceSubModel::nuField() const
{
    #ifdef compre
        nu_=particleCloud_.turbulence().mu() / rho_;
        return nu_;
    #else
        if(switches_[SW_SCALAR_VISCOSITY]) // scalarViscosity=true
            return nu_;
        else
            return particleCloud_.turbulence().nu();
    #endif
}

const volScalarField& forceSubModel::muField() const
{
    #ifdef compre
        return particleCloud_.turbulence().mu();
    #else
        // passing the ref to nu*rho will not work->generate a mu_ field like nu_
        FatalError << "implementation not complete!" << abort(FatalError);

        if(switches_[SW_SCALAR_VISCOSITY]) // scalarViscosity=true
            return nu_*rho_;
        else
            return particleCloud_.turbulence().nu()*rho_;
    #endif
}

const volScalarField& forceSubModel::rhoField() const
{
    return rho_;
}

const volVectorField& forceSubModel::divTauField(const volVectorField& U) const
{
    // calc div(Tau)
    #ifdef compre
        const volScalarField& mu_ = muField();
        divTau_ = -fvc::laplacian(mu_, U) - fvc::div(mu_*dev(fvc::grad(U)().T()));
        return divTau_;
    #else
        const volScalarField& nu_ = nuField();
        const volScalarField& rho_ = rhoField();
        divTau_ = -fvc::laplacian(nu_*rho_, U)- fvc::div(nu_*rho_*dev(fvc::grad(U)().T()));
        return divTau_;
    #endif
}

const volVectorField& forceSubModel::IBDragPerV(const volVectorField& U,const volScalarField& p) const
{
    #ifdef compre
        IBDragPerV_ = muField()*fvc::laplacian(U)-fvc::grad(p);
    #else
        IBDragPerV_ = rhoField()*(nuField()*fvc::laplacian(U)-fvc::grad(p));
    #endif
    return IBDragPerV_;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
