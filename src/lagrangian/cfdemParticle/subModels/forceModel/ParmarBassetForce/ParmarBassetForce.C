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
    nInt_(propsDict_.lookupOrDefault<int>("nIntegral", -1)),
    discOrder_(propsDict_.lookupOrDefault<int>("discretisationOrder", 1)),
    ddtUrelHistSize_(nInt_+discOrder_),
    rHistSize_(nInt_),
    FHistSize_(2*discOrder_),
    ddtUrelHistRegName_(typeName + "ddtUrelHist"),  // indexed as: ddtUrelHist_[particle ID][iDim*ddtUrelHistSize_+iHist]
    rHistRegName_(typeName + "rHist"),              // indexed as: rHist_[particle ID][iHist]
    FHistRegName_(typeName + "FHist"),              // indexed as: ddtUrelHist_[particle ID][iDim*FHistSize_*2+iHist*2+k]
    gH0RegName_(typeName + "gH0"),
    tRefRegName_(typeName + "tRef"),
    mRefRegName_(typeName + "mRef"),
    lRefRegName_(typeName + "lRef"),
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
    forceSubM(0).setSwitchesList(SW_VERBOSE,true); // activate search for verbose switch
    forceSubM(0).readSwitches();

    //Extra switches/settings
    particleCloud_.checkCG(true);

    if (discOrder_ < 1 || discOrder_ > 2)
        FatalError << "Parmar Basset Force: Discretisation order > 2 not implemented!" << abort(FatalError);

    if (nInt_ < 1)
        FatalError << "Parmar Basset Force: nIntegral missing or invalid, must be > 1" << abort(FatalError);

    // allocate particle properties
    particleCloud_.registerParticleProperty<double**>(ddtUrelHistRegName_,3*ddtUrelHistSize_,NOTONCPU,false);
    particleCloud_.registerParticleProperty<double**>(      rHistRegName_,        rHistSize_,NOTONCPU,false);
    particleCloud_.registerParticleProperty<double**>(      FHistRegName_,      6*FHistSize_,NOTONCPU,false);
    particleCloud_.registerParticleProperty<double**>(        gH0RegName_,                1 ,NOTONCPU,false);
    particleCloud_.registerParticleProperty<double**>(       tRefRegName_,                1 ,NOTONCPU,false);
    particleCloud_.registerParticleProperty<double**>(       mRefRegName_,                1 ,NOTONCPU,false);
    particleCloud_.registerParticleProperty<double**>(       lRefRegName_,                1 ,NOTONCPU,false);

    //Append the field names to be probed
    particleCloud_.probeM().initialize(typeName, typeName+".logDat");
    particleCloud_.probeM().vectorFields_.append("ParmarBassetForce"); //first entry must the be the force
    particleCloud_.probeM().vectorFields_.append("Urel");
    particleCloud_.probeM().vectorFields_.append("ddtUrel");
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

}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void ParmarBassetForce::setForce() const
{
    double**& ddtUrelHist_ = particleCloud_.getParticlePropertyRef<double**>(ddtUrelHistRegName_);
    double**& rHist_       = particleCloud_.getParticlePropertyRef<double**>(      rHistRegName_);
    double**& FHist_       = particleCloud_.getParticlePropertyRef<double**>(      FHistRegName_);
    double**& gH0_         = particleCloud_.getParticlePropertyRef<double**>(        gH0RegName_);
    double**& tRef_        = particleCloud_.getParticlePropertyRef<double**>(       tRefRegName_);
    double**& mRef_        = particleCloud_.getParticlePropertyRef<double**>(       mRefRegName_);
    double**& lRef_        = particleCloud_.getParticlePropertyRef<double**>(       lRefRegName_);

    vector position(0,0,0);
    vector Urel(0,0,0);
    vector ddtUrel(0,0,0);

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

    smoothingM().smoothen(Urel_);
    smoothingM().smoothen(ddtUrel_);

    interpolationCellPoint<vector> UrelInterpolator_(Urel_);
    interpolationCellPoint<vector> ddtUrelInterpolator_(ddtUrel_);

    for(int index = 0;index <  particleCloud_.numberOfParticles(); index++)
    {
            vector ParmarBassetForce(0,0,0);
            vector Fshort(0,0,0);
            vector Flong[2]={vector::zero, vector::zero};
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
                    scalar gHratio = gH0_[index][0]/gH;
                    r = gHratio*sqrt(gHratio); // gHratio^1.5, Eq. 3.4

                    if (r<0.25 || r>2.0)
                    {
                        gH0_[index][0]         = NOTONCPU; //reset reference
                        ddtUrelHist_[index][0] = NOTONCPU; //reset ddtU history (only component used for checking nKnown)
                    }

                }

                if (gH0_[index][0]==NOTONCPU)
                {
                    // calculate new reference values
                    scalar tRef = (rs*rs) / nufField[cellI] * (gH*gH) * 4.3354; // (256/pi)^(1/3) = 4.3354 (Eq 3.1)
                    scalar Vs   = rs*rs*rs*M_PI*4/3;
                    scalar mRef = Vs*rhoField[cellI] * gH * 5.2863;             // 9/(2*sqrt(pi))*(256/pi)^(1/6) = 5.2863 (Eq. 3.2)

                    gH0_[index][0]  = gH;
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

                // non-dimensionlise
                Urel    /= mps;
                ddtUrel /= mpss;

                // check length of known history
                int nKnown = 1; // we always know the current step
                for (int j=0; j<ddtUrelHistSize_; j++) // loop over past times
                {
                    if (ddtUrelHist_[index][j] == NOTONCPU)
                        break;
                    else
                        nKnown++;
                }

                //********* short term force computing (full integral) *********//

                // number of values to be used for short term force computing
                int nShort = min(nKnown,nInt_+1);

                // int_0^dt K(r,xi) dxi * ddtU(t) dxi (singularity treated by assuming constant acceleration)
                for (int i=0; i<3; i++) // loop over dimensions
                    Fshort[i] = -calculateK0(r,dt) * ddtUrel[i];

                // int_dt^t0 K(r,xi) * ddtU(t-xi) dxi (trapezoid rule)
                if (nShort>2)
                {
                    for (int j=0; j<(nShort-1); j++) // we don't use the current step here, hence nShort-1
                    {
                        scalar xi = (j+1)*dt;
                        scalar invsqrtK = sqrt(sqrt(xi)) + rHist_[index][j]*xi; // K^-0.5
                        scalar K  = 1./(invsqrtK*invsqrtK); // Eq. 3.4

                        for (int i=0; i<3; i++) // loop over dimensions
                            Fshort[i] -= trapWeight(j,nShort-1) * K * ddtUrelHist_[index][i*ddtUrelHistSize_+j] * dt;
                    }
                }

                //********* long term force computing (differential form) *********//

                // initialise ddtUrel(t0) and Flong(:) as 0 and r(t0) as 1
                if (nKnown == nInt_)
                {
                    // initialise the histories beyond nInt
                    for (int j=nInt_-1; j<rHistSize_; j++) // loop over past times
                        rHist_[index][j] = 1.;

                    for (int j=nInt_-1; j<ddtUrelHistSize_; j++) // loop over past times
                        for (int i=0; i<3; i++) // loop over dimensions
                            ddtUrelHist_[index][i*ddtUrelHistSize_+j] = 0.0;

                    for (int k=0; k<2; k++) // loop over F1, F2
                        for (int j=0; j<FHistSize_; j++) // loop over past times
                            for (int i=0; i<3; i++) // loop over dimensions
                                FHist_[index][i*FHistSize_*2+j*2+k] = 0.0;
                    nKnown = ddtUrelHistSize_+1;
                }

                // solve ODEs
                if (nKnown == ddtUrelHistSize_+1)
                {
                    for (int k=0; k<2; k++) // loop over F1, F2
                    {
                        //calculate coefficients
                        double C[4];
                        calculateCoeffs(k,t0,rHist_[index][nInt_-1],c,chi,C);

                        // solve Eq. 3.20
                        Flong[k] = solveFlongODE(FHist_,ddtUrelHist_,k,C,dt,index);
                    }
                }

                //********* update histories *********//
                update_ddtUrelHist(ddtUrelHist_,ddtUrel,index); // add current dim.less ddtUrel to history
                update_rHist(rHist_,r,index); // add current r to history
                update_FHist(FHist_,Flong[0],Flong[1],index);

                //********* total force *********//

                // sum and convert to N
                ParmarBassetForce = Fshort;
                for (int k=0; k<2; k++) // loop over F1, F2
                    ParmarBassetForce += Flong[k];
                ParmarBassetForce *= newton;

                if (forceSubM(0).verbose() && index >= 0 && index < 2)
                {
                    Pout << "cellI = " << cellI << endl;
                    Pout << "index = " << index << endl;
                    Pout << "Fshort = " << Fshort*newton << endl;
                    Pout << "Flong1 = " << Flong[0]*newton << endl;
                    Pout << "Flong2 = " << Flong[1]*newton << endl;
                    Pout << "Ftotal = " << ParmarBassetForce << endl;
                }

                // Set value fields and write the probe
                if(probeIt_)
                {
                    scalar ReRef = 0.75/(gH0_[index][0]-0.105);

                    #include "setupProbeModelfields.H"
                    vValues.append(ParmarBassetForce);           //first entry must the be the force
                    vValues.append(Urel);
                    vValues.append(ddtUrel);
                    vValues.append(Fshort);
                    vValues.append(Flong[0]);
                    vValues.append(Flong[1]);
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
                gH0_[index][0]         = NOTONCPU; //reset reference
                ddtUrelHist_[index][0] = NOTONCPU; //reset ddtU history (only component used for checking nKnown)
            }

            // write particle based data to global array
            forceSubM(0).partToArray(index,ParmarBassetForce,vector::zero);
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
scalar Foam::ParmarBassetForce::calculateK0(scalar r, scalar dt) const
{
    scalar cbrtr = cbrt(r); // cube root of r
    scalar gamma = cbrtr*sqrt(sqrt(dt));

    /*
    scalar K0 = 2./(9.*pow(r,0.666)) *
            (
                    6.*gamma*gamma/(gamma*gamma*gamma+1)
                    + 2.*sqrt(3.)*atan(sqrt(3.)*gamma/(2.-gamma))
                    + log((gamma*gamma-gamma+1.)/(gamma+1.)/(gamma+1.))
            ); // int_0^dt K(r,xi) dxi
    */  
    scalar K0 = 2./(9.*cbrtr*cbrtr) *
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
    if ( (i==0) || (i==(n-1)) )
        return 0.5;
    else
        return 1.0;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
void Foam::ParmarBassetForce::update_ddtUrelHist(double**& ddtUrelHist_, const vector& ddtUrel, int index) const
{
    for (int i=0; i<3; i++) // loop over dimensions
    {
        for (int j=ddtUrelHistSize_-1; j>0; j--) // loop over past times
            ddtUrelHist_[index][i*ddtUrelHistSize_+j] = ddtUrelHist_[index][i*ddtUrelHistSize_+j-1];

        ddtUrelHist_[index][i*ddtUrelHistSize_] = ddtUrel[i];
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
void Foam::ParmarBassetForce::update_rHist(double**& rHist_, scalar r, int index) const
{
    for (int j=rHistSize_-1; j>0; j--) // loop over past times
        rHist_[index][j] = rHist_[index][j-1];

    rHist_[index][0] = r;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
void Foam::ParmarBassetForce::update_FHist(double**& FHist_, const vector& F1, const vector& F2, int index) const
{
    for (int i=0; i<3; i++) // loop over dimensions
    {
        for (int k=0; k<2; k++) // loop over F1, F2
            for (int j=FHistSize_-1; j>0; j--) // loop over past times
                FHist_[index][i*FHistSize_*2+j*2+k] = FHist_[index][i*FHistSize_*2+(j-1)*2+k];

        FHist_[index][i*FHistSize_*2  ] = F1[i];
        FHist_[index][i*FHistSize_*2+1] = F2[i];
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
vector Foam::ParmarBassetForce::solveFlongODE(double**& FHist_, double**& ddtUrelHist_, int k, double C[4], scalar dt, int index) const
{
    vector Flong = vector::zero;

    if (discOrder_==1)
    {
        for (int i=0; i<3; i++) // loop over dimensions
        {
            Flong[i] = 
                    (
                        (     C[1]/dt+2/(dt*dt)) * FHist_[index][i*FHistSize_*2  +k]
                      - (             1/(dt*dt)) * FHist_[index][i*FHistSize_*2+2+k]
                      - (C[2]+C[3]/dt        ) * ddtUrelHist_[index][i*ddtUrelHistSize_+nInt_-1]
                      + (     C[3]/dt        ) * ddtUrelHist_[index][i*ddtUrelHistSize_+nInt_  ]
                    ) / (C[0]+C[1]/dt+1/(dt*dt)); // Eq. 3.20 using first order temporal discretisation
        }
    }
    if (discOrder_==2)
    {
        for (int i=0; i<3; i++) // loop over dimensions
        {
            Flong[i] = 
                    (
                        (       4*C[1]/(2*dt) + 24/(4*dt*dt)) * FHist_[index][i*FHistSize_*2  +k]
                      - (         C[1]/(2*dt) + 22/(4*dt*dt)) * FHist_[index][i*FHistSize_*2+2+k]
                      + (                        8/(4*dt*dt)) * FHist_[index][i*FHistSize_*2+4+k]
                      - (                        1/(4*dt*dt)) * FHist_[index][i*FHistSize_*2+6+k]

                      - (C[2] + 3*C[3]/(2*dt)               ) * ddtUrelHist_[index][i*ddtUrelHistSize_+nInt_-1]
                      + (       4*C[3]/(2*dt)               ) * ddtUrelHist_[index][i*ddtUrelHistSize_+nInt_  ]
                      - (         C[3]/(2*dt)               ) * ddtUrelHist_[index][i*ddtUrelHistSize_+nInt_+1]
                    ) / (C[0] + 3*C[1]/(2*dt) +  9/(4*dt*dt)); // Eq. 3.20 using second order temporal discretisation
        }
    }
    return Flong;
}

} // End namespace Foam

// ************************************************************************* //

