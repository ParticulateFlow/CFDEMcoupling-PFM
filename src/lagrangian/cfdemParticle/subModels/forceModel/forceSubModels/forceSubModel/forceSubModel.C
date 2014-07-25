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
    nrDefaultSwitches_(3),
    switchesNameList_(wordList(nrDefaultSwitches_)),
    switchesList_(List<Switch>(nrDefaultSwitches_)),
    switches_(List<Switch>(nrDefaultSwitches_))
{
    // init switches lists
    switchesNameList_[0]="treatExplicit";
    switchesNameList_[1]="treatDEM";
    switchesNameList_[2]="implDEM";
    for(int i=0;i<switchesList_.size();i++)
    {
        switchesList_[i]=false;
        switches_[i]=false;
    }

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
    label& index,
    vector& dragTot,
    const vector& dragEx,
    const vector& Ufluid,
    scalar Cd
) const
{
    // forces for CFD
    if(!switches_[1])// !treatDEM
    {
        if(switches_[0]) // treatExplicit
        {
            for(int j=0;j<3;j++)
                myForceM().expForces()[index][j] += dragTot[j];
        }    
        else   //implicit treatment, taking explicit force contribution into account
        {
            for(int j=0;j<3;j++) 
            { 
                myForceM().impForces()[index][j] += dragTot[j] - dragEx[j]; //only consider implicit part!
                myForceM().expForces()[index][j] += dragEx[j];
            }
        }
    }

    // forces for DEM
    if(switches_[2]) // implDEM
    {
        for(int j=0;j<3;j++)
            myForceM().fluidVel()[index][j]=Ufluid[j];

        myForceM().Cds()[index][0]=Cd;
    }
    else
    {
        for(int j=0;j<3;j++) 
            myForceM().DEMForces()[index][j] += dragTot[j];
    }
}

void forceSubModel::explicitInterpCorr
(
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

void forceSubModel::readSwitches() const
{
    Info << "\nforceSubModel:" << myType() << endl;
    forAll(switchesNameList_,i)
    {
        if(switchesList_[i] > 0+SMALL)
        {
            Info << "  looking for " << switchesNameList_[i] << " ..." << endl;
            if (dict_.found(switchesNameList_[i]))
                switches_[i]=Switch(dict_.lookup(switchesNameList_[i]));
                
            Info << "\t" << switchesNameList_[i] << " = " << switches_[i] << endl;
        }        
    }
    Info << endl;

    if(switches_[2]) // implDEM=true
    {
        // communicate implDEM to particleCloud
        particleCloud_.impDEMdrag_=true;

        // do sanity check
        if(switches_[0]) // treatExplicit=true
        {
            Warning<< "please check your settings, treatExplicit together with implDEM does not work! (using treatExplicit=false)" << endl;
            switches_[0]=false;
        }
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
