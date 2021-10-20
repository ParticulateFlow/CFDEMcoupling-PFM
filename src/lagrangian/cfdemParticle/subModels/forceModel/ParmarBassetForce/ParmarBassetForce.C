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

#include "ParmarBassetForce.H"
#include "addToRunTimeSelectionTable.H"
#include "smoothingModel.H"
#include "constDiffSmoothing.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

#define NOTONCPU 9999

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(ParmarBassetForce, 0);

addToRunTimeSelectionTable
(
    forceModel,
    ParmarBassetForce,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
ParmarBassetForce::ParmarBassetForce
(
    const dictionary& dict,
    cfdemCloud& sm
)
:
    forceModel(dict,sm),
    propsDict_(dict.subDict(typeName + "Props")),
    velFieldName_(propsDict_.lookup("velFieldName")),
    U_(sm.mesh().lookupObject<volVectorField> (velFieldName_)),
    UsFieldName_(propsDict_.lookup("granVelFieldName")),
    Us_(sm.mesh().lookupObject<volVectorField> (UsFieldName_)),
    nInt_(readLabel(propsDict_.lookup("nIntegral"))),
    discOrder_(readLabel(propsDict_.lookup("discretisationOrder"))),
    nHist_(nInt_+discOrder_+1),
    ddtUrelHist_(nHist_,NULL),                     // UrelHist_[ndt in past][particle ID][dim]
    rHist_(nHist_,NULL),                           // rHist_[ndt in past][particle ID][0]
    FHist_(2,List<double**>(2*discOrder_+1,NULL)), // FHist_[k={1,2}-1][ndt in past][particle ID][dim]
    gH0_(NULL),
    tRef_(NULL),
    mRef_(NULL),
    lRef_(NULL),
    Urel_
    (   IOobject
        (
            "Urel",
            sm.mesh().time().timeName(),
            sm.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        sm.mesh(),
        dimensionedVector("zero", dimensionSet(0,1,-1,0,0,0,0), vector::zero)
    ),
    ddtUrel_
    (   IOobject
        (
            "ddtUrel",
            sm.mesh().time().timeName(),
            sm.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        sm.mesh(),
        dimensionedVector("zero", dimensionSet(0,1,-2,0,0,0,0), vector::zero)
    ),
    smoothingModel_
    (
        smoothingModel::New
        (
            propsDict_,
            sm
        )
    )
{

    // init force sub model
    setForceSubModels(propsDict_);
    // define switches which can be read from dict
    forceSubM(0).setSwitchesList(SW_TREAT_FORCE_EXPLICIT,true); // activate treatExplicit switch
    forceSubM(0).setSwitchesList(SW_INTERPOLATION,true); // activate search for interpolate switch
    forceSubM(0).readSwitches();

    //Extra switches/settings
    particleCloud_.checkCG(true);

    if (discOrder_ < 1 || discOrder_ > 2)
        FatalError << "Parmar Basset Force: Discretisation order > 2 not implemented!" << abort(FatalError);

    //Append the field names to be probed
    particleCloud_.probeM().initialize(typeName, typeName+".logDat");
    particleCloud_.probeM().vectorFields_.append("ParmarBassetForce"); //first entry must the be the force
    particleCloud_.probeM().vectorFields_.append("Urel");
    particleCloud_.probeM().vectorFields_.append("ddtUrel");
    //
    particleCloud_.probeM().vectorFields_.append("UrelNoSmooth");
    particleCloud_.probeM().vectorFields_.append("ddtUrelNoSmooth");
    //
    particleCloud_.probeM().vectorFields_.append("Fshort");
    particleCloud_.probeM().vectorFields_.append("Flong1");
    particleCloud_.probeM().vectorFields_.append("Flong2");
    particleCloud_.probeM().scalarFields_.append("ReRef");
    particleCloud_.probeM().scalarFields_.append("tRef");
    particleCloud_.probeM().scalarFields_.append("mRef");
    particleCloud_.probeM().scalarFields_.append("lRef");
    particleCloud_.probeM().scalarFields_.append("r");
    particleCloud_.probeM().scalarFields_.append("t0");
    particleCloud_.probeM().scalarFields_.append("nKnown");
    particleCloud_.probeM().writeHeader();

    Urel_ = Us_ - U_;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

ParmarBassetForce::~ParmarBassetForce()
{
    particleCloud_.dataExchangeM().destroy(gH0_,    1);
    particleCloud_.dataExchangeM().destroy(tRef_,   1);
    particleCloud_.dataExchangeM().destroy(mRef_,   1);
    particleCloud_.dataExchangeM().destroy(lRef_,   1);

    for (int i=0; i<nHist_; i++)
    {
        particleCloud_.dataExchangeM().destroy(ddtUrelHist_[i],3);
        particleCloud_.dataExchangeM().destroy(rHist_      [i],1);
    }

    for (int i=0; i<2*discOrder_+1; i++)
        for (int k=0; k<2; k++)
            particleCloud_.dataExchangeM().destroy(FHist_[k][i],3);

}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void ParmarBassetForce::setForce() const
{
    // allocate arrays
    if(particleCloud_.numberOfParticlesChanged())
        reAllocArrays();

    vector position(0,0,0);
    vector Urel(0,0,0);
    vector ddtUrel(0,0,0);

    scalar t0min    = 0.00;

    const volScalarField& nufField = forceSubM(0).nuField();
    const volScalarField& rhoField = forceSubM(0).rhoField();

    double c[2][4][3]   = {{{0.004045,0.011788,-0.851409},
                            {0.184350,0.009672,-0.425408},
                            {0.000746,0.014390,-1.262288},
                            {0.021708,0.016611,-0.867352}},

                           {{0.536783,0.005453,-1.430163},
                            {1.683212,0.003598,-0.773872},
                            {0.316350,0.004367,-1.423649},
                            {0.420423,0.000677,-0.670551}}
                           };

    double chi[2][4][2] = {{{1.464787,  0.422501},
                            {0.785093, -0.078301},
                            {0.585753, -0.118737},
                            {-0.077407,-0.016665}},

                           {{1.001435, -0.003566},
                            {0.415345, -0.136697},
                            {0.536055, -0.205041},
                            {-0.032671,-0.021047}}
                           };

    #include "setupProbeModel.H"

    Urel_ = Us_ - U_;
    ddtUrel_ = fvc::ddt(Us_) - fvc::ddt(U_) - (Us_ & fvc::grad(U_));

    //
    volVectorField UrelNoSmooth_ = Urel_;
    volVectorField ddtUrelNoSmooth_ = ddtUrel_;
    //

    smoothingM().smoothen(Urel_);
    smoothingM().smoothen(ddtUrel_);

    interpolationCellPoint<vector> UrelInterpolator_(Urel_);
    interpolationCellPoint<vector> ddtUrelInterpolator_(ddtUrel_);

    //
    interpolationCellPoint<vector> UrelNoSmoothInterpolator_(UrelNoSmooth_);
    interpolationCellPoint<vector> ddtUrelNoSmoothInterpolator_(ddtUrelNoSmooth_);
    //

    for(int index = 0;index <  particleCloud_.numberOfParticles(); index++)
    {
            vector ParmarBassetForce(0,0,0);
            vector Fshort(0,0,0);
            label  cellI = particleCloud_.cellIDs()[index][0];

            if (cellI > -1) // particle Found
            {
                // particle position (m)
                position = particleCloud_.position(index);

                // relative velocity (m/s)
                if(forceSubM(0).interpolation())
                    Urel = UrelInterpolator_.interpolate(position,cellI);
                else
                    Urel = Urel_[cellI];

                // acceleration (m/s2)
                if(forceSubM(0).interpolation())
                    ddtUrel = ddtUrelInterpolator_.interpolate(position,cellI);
                else
                    ddtUrel = ddtUrel_[cellI];

                //********* set/get reference point *********//

                scalar rs = particleCloud_.radius(index);
                scalar Re = 2.*rs*mag(Urel)/nufField[cellI];
                scalar gH = (0.75 + 0.105*Re)/(Re+SMALL); // Eq. 3.3
                scalar r;

                if (gH0_[index][0]!=NOTONCPU)
                {
                    r = pow(gH0_[index][0]/gH,1.5); // Eq. 3.4

                    if (r<0.25 || r>2.0)
                    {
                        gH0_[index][0]            = NOTONCPU; //reset reference
                        ddtUrelHist_[0][index][0] = NOTONCPU; //reset ddtU history (only component used for checking nKnown)
                    }

                }

                if (gH0_[index][0]==NOTONCPU)
                {
                    // calculate new reference values
                    scalar tRef = (rs*rs) / nufField[cellI] * (gH*gH) * 4.3354; // (256/pi)^(1/3) = 4.3354 (Eq 3.1)
                    scalar Vs   = rs*rs*rs*M_PI*4/3;
                    scalar mRef = Vs*rhoField[cellI] * gH * 5.2863;             // 9/(2*sqrt(pi))*(256/pi)^(1/6) = 5.2863 (Eq. 3.2)

                    gH0_[index][0] = gH;
                    tRef_[index][0] = tRef;
                    mRef_[index][0] = mRef;
                    lRef_[index][0] = rs;
                    r = 1;
                }

                scalar mps    = lRef_[index][0]/tRef_[index][0]; // m/s
                scalar mpss   = mps/tRef_[index][0]; // m/s^2
                scalar newton = mpss*mRef_[index][0]; // kg*m/s^2
                scalar dt     = U_.mesh().time().deltaT().value() / tRef_[index][0];  // dim.less
                scalar t0     = nInt_*dt; // dim.less

                //********* update histories *********//

                // non-dimensionlise
                Urel    /= mps;
                ddtUrel /= mpss;

                // update ddtUrel  and r history
                update_ddtUrelHist(ddtUrel,index); // add current dim.less ddtUrel to history
                update_rHist(r,index); // add current r to history

                // warning and reset for too small t0
                if (t0<t0min)
                {
                    Pout << "ParmarBassetForce WARNING: t0 = " << t0 << " at ID = " << index <<  endl;
                    gH0_[index][0]            = NOTONCPU; //reset reference
                    ddtUrelHist_[0][index][0] = NOTONCPU; //reset ddtU history (only component used for checking nKnown)
                }

                // check length of known history
                int nKnown = 0;
                for (int j=0; j<nHist_; j++) // loop over past times
                {
                    if (ddtUrelHist_[j][index][0] == NOTONCPU)
                        break;
                    else
                        nKnown++;
                }

                //********* short term force computing (full integral) *********//

                // number of values to be used for short term force computing
                int nShort = min(nKnown,nInt_+1);

                // int_0^dt K(r,xi) dxi * ddtU(t) dxi (singularity treated by assuming constant acceleration)
                if (nShort>0)
                {
                    for (int i=0; i<3; i++) // loop over dimensions
                        Fshort[i] = -calculateK0(r,dt) * ddtUrelHist_[0][index][i];
                }

                // int_dt^t0 K(r,xi) * ddtU(t-xi) dxi (trapezoid rule)
                if (nShort>2)
                {
                    for (int j=1; j<nShort; j++)
                    {
                        scalar xi = j*dt;
                        scalar K  = pow((pow(xi,.25) + rHist_[j][index][0]*xi),-2.); // Eq. 3.4

                        for (int i=0; i<3; i++) // loop over dimensions
                            Fshort[i] -= trapWeight(j,nShort) * K * ddtUrelHist_[j][index][i] * dt;
                    }
                }

                //********* long term force computing (differential form) *********//

                // update F1, F2 history
                update_FHist(vector::zero,vector::zero,index);

                // initialise ddtUrel(t0) and Flong(:) as 0 and r(t0) as 1
                if (nKnown == nInt_)
                {
                    for (int j=nInt_; j<nHist_; j++) // loop over past times
                    {
                        rHist_[j][index][0] = 1.;
                        for (int i=0; i<3; i++) // loop over dimensions
                            ddtUrelHist_[j][index][i] = 0.0;
                    }

                    for (int j=0; j<2*discOrder_+1; j++) // loop over past times
                        for (int i=0; i<3; i++) // loop over dimensions
                            for (int k=0; k<2; k++) // loop over F1, F2
                                FHist_[k][j][index][i] = 0.0;
                    nKnown = nHist_;
                }

                // solve ODEs
                if (nKnown == nHist_)
                {
                    for (int k=0; k<2; k++) // loop over F1, F2
                    {
                        //calculate coefficients
                        double C[4];
                        calculateCoeffs(k,t0,rHist_[nInt_][index][0],c,chi,C);

                        // solve Eq. 3.20
                        solveFlongODE(k,C,dt,index);
                    }
                }

                //********* total force *********//

                // sum and convert to N
                for (int i=0; i<3; i++) // loop over dimensions
                {
                    ParmarBassetForce[i] = Fshort[i];
                    for (int k=0; k<2; k++) // loop over F1, F2
                        ParmarBassetForce[i] += FHist_[k][0][index][i];
                }
                ParmarBassetForce *= newton;

                // Set value fields and write the probe
                if(probeIt_)
                {
                    scalar ReRef = 0.75/(gH0_[index][0]-0.105);

                    vector Flong1(0,0,0);
                    vector Flong2(0,0,0);

                    for (int i=0; i<3; i++) // loop over dimensions
                    {
                        Flong1[i] = FHist_[0][0][index][i];
                        Flong2[i] = FHist_[1][0][index][i];
                    }

                    //
                    // relative velocity (m/s)
                    vector UrelNoSmooth;
                    vector ddtUrelNoSmooth;

                    if(forceSubM(0).interpolation())
                        UrelNoSmooth = UrelNoSmoothInterpolator_.interpolate(position,cellI);
                    else
                        UrelNoSmooth = UrelNoSmooth_[cellI];

                    // acceleration (m/s2)
                    if(forceSubM(0).interpolation())
                        ddtUrelNoSmooth = ddtUrelNoSmoothInterpolator_.interpolate(position,cellI);
                    else
                        ddtUrelNoSmooth = ddtUrelNoSmooth_[cellI];

                    UrelNoSmooth /= mps;
                    ddtUrelNoSmooth /= mpss;
                    //

                    #include "setupProbeModelfields.H"
                    vValues.append(ParmarBassetForce);           //first entry must the be the force
                    vValues.append(Urel);
                    vValues.append(ddtUrel);
                    //
                    vValues.append(UrelNoSmooth);
                    vValues.append(ddtUrelNoSmooth);
                    //
                    vValues.append(Fshort);
                    vValues.append(Flong1);
                    vValues.append(Flong2);
                    sValues.append(ReRef);
                    sValues.append(tRef_[index][0]);
                    sValues.append(mRef_[index][0]);
                    sValues.append(lRef_[index][0]);
                    sValues.append(r);
                    sValues.append(t0);
                    sValues.append(nKnown);
                    particleCloud_.probeM().writeProbe(index, sValues, vValues);
                }
            }
            else
            {
                // not on CPU
                gH0_[index][0]            = NOTONCPU; //reset reference
                ddtUrelHist_[0][index][0] = NOTONCPU; //reset ddtU history (only component used for checking nKnown)
            }

            // write particle based data to global array
            forceSubM(0).partToArray(index,ParmarBassetForce,vector::zero);
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
void Foam::ParmarBassetForce::reAllocArrays() const
{
    Pout << "ParmarBassetForce::reAllocArrays..." << endl;

    particleCloud_.dataExchangeM().allocateArray(gH0_,     NOTONCPU,1);
    particleCloud_.dataExchangeM().allocateArray(tRef_,    NOTONCPU,1);
    particleCloud_.dataExchangeM().allocateArray(mRef_,    NOTONCPU,1);
    particleCloud_.dataExchangeM().allocateArray(lRef_,    NOTONCPU,1);

    for (int i=0; i<nHist_; i++)
    {
        particleCloud_.dataExchangeM().allocateArray(ddtUrelHist_[i],NOTONCPU,3);
        particleCloud_.dataExchangeM().allocateArray(rHist_      [i],NOTONCPU,1);
    }

    for (int i=0; i<2*discOrder_+1; i++)
        for (int k=0; k<2; k++)
            particleCloud_.dataExchangeM().allocateArray(FHist_[k][i],0.0,3);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
scalar Foam::ParmarBassetForce::calculateK0(scalar r, scalar dt) const
{
    scalar gamma = pow(r,0.333)*pow(dt,0.25);

    /*
    scalar K0 = 2./(9.*pow(r,0.666)) *
            (
                    6.*gamma*gamma/(gamma*gamma*gamma+1)
                    + 2.*sqrt(3.)*atan(sqrt(3.)*gamma/(2.-gamma))
                    + log((gamma*gamma-gamma+1.)/(gamma+1.)/(gamma+1.))
            ); // int_0^dt K(r,xi) dxi
    */
    scalar K0 = 2./(9.*pow(r,0.666)) *
            (
              -  0.37224 * gamma
              + 12.1587  * gamma*gamma
              -  6.4884  * gamma*gamma*gamma
            ); // third order approximation to the eq. above

    return K0;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
scalar Foam::ParmarBassetForce::trapWeight(int i, int n) const
{
    if ( (i==1) || (i==(n-1)) )
        return 0.5;
    else
        return 1.0;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
void Foam::ParmarBassetForce::update_ddtUrelHist(const vector& ddtUrel, int index) const
{
    for (int i=0; i<3; i++) // loop over dimensions
    {
        for (int j=nHist_-1; j>0; j--) // loop over past times
            ddtUrelHist_[j][index][i] = ddtUrelHist_[j-1][index][i];

        ddtUrelHist_[0][index][i] = ddtUrel[i];
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
void Foam::ParmarBassetForce::update_rHist(scalar r, int index) const
{
    for (int j=nHist_-1; j>0; j--) // loop over past times
        rHist_[j][index][0] = rHist_[j-1][index][0];

    rHist_[0][index][0] = r;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
void Foam::ParmarBassetForce::update_FHist(const vector& F1, const vector& F2, int index) const
{
    for (int i=0; i<3; i++) // loop over dimensions
    {
        for (int k=0; k<2; k++) // loop over F1, F2
            for (int j=2*discOrder_; j>0; j--) // loop over past times
                FHist_[k][j][index][i] = FHist_[k][j-1][index][i];

        FHist_[0][0][index][i] = F1[i];
        FHist_[1][0][index][i] = F2[i];
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
void Foam::ParmarBassetForce::calculateCoeffs(int k, scalar t0, scalar r, double c[2][4][3], double chi[2][4][2], double C[4]) const
{
    for (int j=0; j<4; j++) // loop over terms
    {
        C[j]  = (c[k][j][0] * pow((t0 + c[k][j][1]),c[k][j][2]));    //t0 dependency
        C[j] *= (1 + chi[k][j][0]*(r-1) + chi[k][j][1]*(r-1)*(r-1)); // r dependency
    }
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
void Foam::ParmarBassetForce::solveFlongODE(int k, double C[4], scalar dt, int index) const
{
    if (discOrder_==1)
    {
        for (int i=0; i<3; i++) // loop over dimensions
        {
            FHist_[k][0][index][i] =
                    (
                        (     C[1]/dt+2/(dt*dt)) * FHist_[k][1][index][i]
                      - (             1/(dt*dt)) * FHist_[k][2][index][i]
                      - (C[2]+C[3]/dt        ) * ddtUrelHist_[nInt_  ][index][i]
                      + (     C[3]/dt        ) * ddtUrelHist_[nInt_+1][index][i]
                    ) / (C[0]+C[1]/dt+1/(dt*dt)); // Eq. 3.20 using first order temporal discretisation
        }
    }
    if (discOrder_==2)
    {
        for (int i=0; i<3; i++) // loop over dimensions
        {
            FHist_[k][0][index][i] =
                    (
                        (       4*C[1]/(2*dt) + 24/(4*dt*dt)) * FHist_[k][1][index][i]
                      - (         C[1]/(2*dt) + 22/(4*dt*dt)) * FHist_[k][2][index][i]
                      + (                        8/(4*dt*dt)) * FHist_[k][3][index][i]
                      - (                        1/(4*dt*dt)) * FHist_[k][4][index][i]

                      - (C[2] + 3*C[3]/(2*dt)               ) * ddtUrelHist_[nInt_  ][index][i]
                      + (       4*C[3]/(2*dt)               ) * ddtUrelHist_[nInt_+1][index][i]
                      - (         C[3]/(2*dt)               ) * ddtUrelHist_[nInt_+2][index][i]
                    ) / (C[0] + 3*C[1]/(2*dt) +  9/(4*dt*dt)); // Eq. 3.20 using second order temporal discretisation
        }
    }
}

void Foam::ParmarBassetForce::rescaleHist(scalar tScale, scalar mScale, scalar lScale, scalar rScale, int index) const
{
    for (int i=0; i<3; i++) // loop over dimensions
    {
        // rescale ddtU history
        for (int j=0; j<nHist_; j++) // loop over past times
            if (ddtUrelHist_[j][index][i] != NOTONCPU)
                ddtUrelHist_[j][index][i] *= lScale/(tScale*tScale);

        // rescale F1, F2 history
        for (int k=0; k<2; k++) // loop over F1, F2
            for (int j=0; j<(2*discOrder_+1); j++) // loop over past times
                FHist_[k][j][index][i] *= mScale*lScale/(tScale*tScale);
    }
    // rescale r history
    for (int j=0; j<nHist_; j++) // loop over past times
        if (rHist_[j][index][0] != NOTONCPU)
            rHist_[j][index][0] /= pow(rScale,1.5);
}

} // End namespace Foam

// ************************************************************************* //

