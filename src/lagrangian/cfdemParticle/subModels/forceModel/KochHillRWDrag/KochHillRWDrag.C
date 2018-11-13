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

#include "KochHillRWDrag.H"
#include "addToRunTimeSelectionTable.H"
#include "dataExchangeModel.H"

//#include <mpi.h>
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(KochHillRWDrag, 0);

addToRunTimeSelectionTable
(
    forceModel,
    KochHillRWDrag,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
KochHillRWDrag::KochHillRWDrag
(
    const dictionary& dict,
    cfdemCloud& sm
)
:
    forceModel(dict,sm),
    propsDict_(dict.subDict(typeName + "Props")),
    verbose_(propsDict_.found("verbose")),
    velFieldName_(propsDict_.lookup("velFieldName")),
    U_(sm.mesh().lookupObject<volVectorField> (velFieldName_)),
    voidfractionFieldName_(propsDict_.lookup("voidfractionFieldName")),
    voidfraction_(sm.mesh().lookupObject<volScalarField> (voidfractionFieldName_)),
    UsFieldName_(propsDict_.lookupOrDefault("granVelFieldName",word("Us"))),
    UsField_(sm.mesh().lookupObject<volVectorField> (UsFieldName_)),
    interpolation_(propsDict_.found("interpolation")),
    scale_(1.),
    randomTauE_(propsDict_.found("randomTauE")),
    partTime_(NULL),
    partUfluct_(NULL),
    RanGen_(label(0))
{

    // init force sub model
    setForceSubModels(propsDict_);

    // define switches which can be read from dict
    forceSubM(0).setSwitchesList(SW_TREAT_FORCE_EXPLICIT,true); // activate treatExplicit switch
    forceSubM(0).setSwitchesList(SW_IMPL_FORCE_DEM,true); // activate implDEM switch

    // read those switches defined above, if provided in dict
    forceSubM(0).readSwitches();

    particleCloud_.checkCG(true);

    if (propsDict_.found("scale"))
        scale_=scalar(readScalar(propsDict_.lookup("scale")));

    if (propsDict_.found("cl"))
        cl_= readScalar(propsDict_.lookup("cl"));
    if (propsDict_.found("rhoP"))
        rhoP_= readScalar(propsDict_.lookup("rhoP"));


//    if (particleCloud_.dataExchangeM().maxNumberOfParticles() > 0)
//    {
        //allocate memory
        particleCloud_.dataExchangeM().allocateArray(partTime_,0.,1);
        particleCloud_.dataExchangeM().allocateArray(partUfluct_,0.,3);
//    }

    //Pout << "RW-TEST: maxNumberOfParticles() == " << particleCloud_.dataExchangeM().maxNumberOfParticles() << endl; // TEST-Output

}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

KochHillRWDrag::~KochHillRWDrag()
{
    particleCloud_.dataExchangeM().destroy(partTime_, 1);
    particleCloud_.dataExchangeM().destroy(partUfluct_, 3);
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void KochHillRWDrag::setForce() const
{

    // realloc the arrays
    reAllocArrays();

    if (scale_ > 1.0)
    {
        Info << "KochHillRW using scale = " << scale_ << endl;
    }
    else if (particleCloud_.cg() > 1.0)
    {
        scale_=particleCloud_.cg();
        Info << "KochHillRW using scale from liggghts cg = " << scale_ << endl;
    }

    const volScalarField& nufField = forceSubM(0).nuField();
    const volScalarField& mufField = forceSubM(0).muField();
    const volScalarField& rhoField = forceSubM(0).rhoField();

    vector position(0,0,0);
    scalar voidfraction(1);
    vector Ufluid(0,0,0);
    vector drag(0,0,0);
    vector dragExplicit(0,0,0);
    scalar dragCoefficient(0);
    label cellI = 0;

    vector Us(0,0,0);
    vector Ur(0,0,0);
    scalar ds(0);
    scalar nuf(0);
    scalar rho(0);
    scalar magUr(0);
    scalar Rep(0);
    scalar Vs(0);
    scalar volumefraction(0);
    scalar betaP(0);

    // [Random walk model]
    /* TODO: 1. read cl (ansys man) from dictionary or replace by c_mu from cfd-online
     *       2. get particle density for particle relaxation time
     */
    scalar k(0);
    scalar epsilon(0);
    scalar timeE(0);
    scalar timeCr(0);
    scalar lengthE(0);
    scalar tauPart(0);
    scalar mu(0);

    // [Random walk model]
    /* TODO: 1. k and epsilon inside of compiler if
     *       2. Check which turbulence model
     *       3. Interpolation of k and epsilon?
     *       4. Time step check!
     */
    const volScalarField& kField = particleCloud_.turbulence().k(); // TODO: inside of compiler if?
    const volScalarField& epsilonField = particleCloud_.turbulence().epsilon(); // TODO: check which turbulence model!
    scalar t = particleCloud_.mesh().time().value();
    scalar deltaT = particleCloud_.mesh().time().deltaT().value();

    //word  test = particleCloud_.turbulence().turbulenceModelName;

    interpolationCellPoint<scalar> voidfractionInterpolator_(voidfraction_);
    interpolationCellPoint<vector> UInterpolator_(U_);

    //Info << "RW-TEST: We are in setForce() at t = " << t << endl; // TEST-Output

    for (int index = 0; index<particleCloud_.numberOfParticles(); ++index)
    {
        //if (mask[index][0])
        //{
            cellI = particleCloud_.cellIDs()[index][0];
            drag = vector::zero;
            dragExplicit = vector::zero;
            dragCoefficient = 0;
            betaP = 0;
            Vs = 0;
            Ufluid = vector::zero;

            // Pout << "RW-TEST: cellI = " << cellI << endl; // TEST-Output
            if (cellI > -1) // particle Found
            {
                if (interpolation_)
                {
                    position = particleCloud_.position(index);
                    voidfraction = voidfractionInterpolator_.interpolate(position,cellI);
                    Ufluid = UInterpolator_.interpolate(position,cellI);
                    //Ensure interpolated void fraction to be meaningful
                    // Info << " --> voidfraction: " << voidfraction << endl;
                    if (voidfraction > 1.0) voidfraction = 1.0;
                    else if (voidfraction < 0.4) voidfraction = 0.4;
                }
                else
                {
                    voidfraction = particleCloud_.voidfraction(index);
                    Ufluid = U_[cellI];
                }

                Us = particleCloud_.velocity(index);

                // [Random walk model]
                /* TODO: 1. cl_ or c_mu?! Hardcoded c_mu = 0.09 in the eddy length calculation.
                 */

                // -------------------------
                // TEST: Time boundaries
                // Can be removed when partTime_ is synced
                k = kField[cellI];
                epsilon = epsilonField[cellI];
                timeE = 2.*cl_*k/(epsilon+SMALL);

                if (partTime_[index][0] > t+10.*timeE)
                    partTime_[index][0] = t;
                // -------------------------

                ds = 2.*particleCloud_.radius(index);

                //Pout << "RW-TEST: t = " << t << " partTime_ = " << partTime_[index][0] << endl; // TEST-Output
                if (t >= partTime_[index][0])
                {
                    mu = mufField[cellI];
                    k = kField[cellI];
                    epsilon = epsilonField[cellI];

                    // Pout << "RW-TEST: mu = " << mu << endl;  // TEST-Output

                    // calculate the eddy-lifetime and the particle crossing time
                    tauPart = rhoP_*ds*ds/(18.*mu+SMALL);
                    lengthE = pow(0.09,0.63)*pow(k,1.5)/(epsilon+SMALL); //c_mu from turbulence?!

                    // two possible ways to calculate the eddy-life-time
                    if (randomTauE_)
                    {
                        timeE = - cl_*k/(epsilon+SMALL)*log(RanGen_.scalar01());
                    }
                    else
                    {
                        timeE = 2.*cl_*k/(epsilon+SMALL);
                        // timeE = lengthE / sqrt(2.*k/3.); // cfd-online version
                    }

                    // particle crossing time and determine the min. time step
                    Ur = Ufluid-Us;
                    magUr = mag(Ur);

                    scalar threshold = lengthE/(tauPart*magUr+SMALL);
                    scalar minDeltaT;
                    if (threshold < 1.)
                    {
                        timeCr = -tauPart*log(1.-threshold);
                        minDeltaT =  min(timeE,timeCr);
                    }
                    else
                    {
                        minDeltaT = timeE;
                    }
                    //Pout << "RW-TEST: timeE = " << timeE << " timeCR = " << timeCr << endl; // TEST-Output

                    // calculate time step of next update
                    partTime_[index][0] = t + minDeltaT;
                    if (minDeltaT < deltaT)
                        Warning << "Random Walk Model: Simulation time step (" << deltaT << ") is bigger than the particle eddy interaction time (" << minDeltaT <<")! " << endl;

                    // update turbulent velocity part and
                    // modify current fluid velocity
                    for (int dim=0; dim<3; dim++)
                    {
                        partUfluct_[index][dim] = RanGen_.GaussNormal()*sqrt(2.*k/3.);
                        //Pout << "RW-TEST: Ufluid[" << dim << "] = " << Ufluid[dim] << " Ufluct = " << partUfluct_[index][dim] << " k = " << k << endl; // TEST-Output
                        Ufluid[dim] = Ufluid[dim] + partUfluct_[index][dim];
                    }
                }
                else
                {
                    // no update of the turbulent velocity part
                    // modify current fluid velocity
                    for(int dim=0; dim<3; dim++)
                        Ufluid[dim] = Ufluid[dim] + partUfluct_[index][dim];
                }

                Ur = Ufluid-Us;
                magUr = mag(Ur);

                // -----------------

                scalar ds_scaled = ds/scale_;
                nuf = nufField[cellI];
                rho = rhoField[cellI];
                Rep = 0.;
                Vs = ds*ds*ds*M_PI/6.; // sphere volume
                volumefraction = 1.-voidfraction+SMALL;

                if (magUr > 0.)
                {
                    // calc particle Re Nr
                    Rep = ds_scaled*voidfraction*magUr/(nuf+SMALL);

                    // calc model coefficient F0
                    scalar F0 = 0.;
                    if(volumefraction < 0.4)
                    {
                        F0 = (1. + 3.*sqrt(volumefraction/2.) + (135./64.)*volumefraction*log(volumefraction)
                              + 16.14*volumefraction
                             )/
                             (1. + 0.681*volumefraction - 8.48*sqr(volumefraction)
                              + 8.16*volumefraction*volumefraction*volumefraction
                             );
                    }
                    else
                    {
                        F0 = 10.*volumefraction/(voidfraction*voidfraction*voidfraction);
                    }

                    // calc model coefficient F3
                    scalar F3 = 0.0673 + 0.212*volumefraction + 0.0232/pow(voidfraction,5);

                    //Calculate F in the formulation of van der Hoef et al. (JFM 528:233-254)
                    scalar F = voidfraction * (F0 + 0.5*F3*Rep);

                    // calc drag model coefficient betaP
                    betaP = (18.*nuf*rho/(ds_scaled*ds_scaled))*voidfraction*F;

                    // calc particle's drag
                    dragCoefficient = Vs*betaP;//*scaleDrag_;
                    if (modelType_ == "B")
                        dragCoefficient /= voidfraction;

                    drag = dragCoefficient * Ur;

                    // explicitCorr
                    forceSubM(0).explicitCorr(drag,dragExplicit,dragCoefficient,Ufluid,U_[cellI],Us,UsField_[cellI],verbose_);
                }

                if (verbose_ && index >= 0 && index < 2)
                {
                    Pout << "index = " << index << endl;
                    Pout << "Us = " << Us << endl;
                    Pout << "Ur = " << Ur << endl;
                    Pout << "ds = " << ds << endl;
                    Pout << "ds/scale = " << ds_scaled << endl;
                    Pout << "rho = " << rho << endl;
                    Pout << "nuf = " << nuf << endl;
                    Pout << "voidfraction = " << voidfraction << endl;
                    Pout << "Rep = " << Rep << endl;
                    Pout << "drag = " << drag << endl;
                }
            }
            // write particle based data to global array
            forceSubM(0).partToArray(index,drag,dragExplicit,Ufluid,dragCoefficient);
        //}
    }
}


void KochHillRWDrag::reAllocArrays() const
{
    if (particleCloud_.numberOfParticlesChanged())
    {
        particleCloud_.dataExchangeM().allocateArray(partTime_,0.0,1);  // field/initVal/with/lenghtFromLigghts
        particleCloud_.dataExchangeM().allocateArray(partUfluct_,0.0,3);
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
