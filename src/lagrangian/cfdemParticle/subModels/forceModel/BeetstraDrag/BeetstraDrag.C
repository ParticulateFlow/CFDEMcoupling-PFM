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
    Contributing author: Paul Kieckhefen, TU Hamburg, Germany

\*---------------------------------------------------------------------------*/

#include "error.H"

#include "BeetstraDrag.H"
#include "addToRunTimeSelectionTable.H"
#include "averagingModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(BeetstraDrag, 0);

addToRunTimeSelectionTable
(
    forceModel,
    BeetstraDrag,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
BeetstraDrag::BeetstraDrag
(
    const dictionary& dict,
    cfdemCloud& sm
)
:
    forceModel(dict,sm),
    propsDict_(dict.subDict(typeName + "Props")),
    multiTypes_(false),
    velFieldName_(propsDict_.lookup("velFieldName")),
    U_(sm.mesh().lookupObject<volVectorField> (velFieldName_)),
    voidfractionFieldName_(propsDict_.lookup("voidfractionFieldName")),
    voidfraction_(sm.mesh().lookupObject<volScalarField> (voidfractionFieldName_)),
    minVoidfraction_(propsDict_.lookupOrDefault<scalar>("minVoidfraction",0.1)),
    UsFieldName_(propsDict_.lookup("granVelFieldName")),
    UsField_(sm.mesh().lookupObject<volVectorField> (UsFieldName_)),
    scaleDia_(1.),
    typeCG_(propsDict_.lookupOrDefault<scalarList>("coarseGrainingFactors",scalarList(1,1.0))),
    scaleDrag_(1.),
    rhoP_(0.),
    rho_(0.),
    Lc2_(0.),
    dPrim_(0.),
    nuf_(0.),
    g_(9.81),
    k_(0.05),
    useGC_(false),
    usePC_(false)
{
    //Append the field names to be probed
    particleCloud_.probeM().initialize(typeName, typeName+".logDat");
    particleCloud_.probeM().vectorFields_.append("dragForce"); //first entry must  be the force
    particleCloud_.probeM().vectorFields_.append("Urel");
    particleCloud_.probeM().scalarFields_.append("Rep");
    particleCloud_.probeM().scalarFields_.append("betaP");
    particleCloud_.probeM().scalarFields_.append("voidfraction");
    particleCloud_.probeM().writeHeader();

    // init force sub model
    setForceSubModels(propsDict_);
    // define switches which can be read from dict
    forceSubM(0).setSwitchesList(SW_TREAT_FORCE_EXPLICIT,true); // activate treatExplicit switch
    forceSubM(0).setSwitchesList(SW_IMPL_FORCE_DEM,true); // activate implDEM switch
    forceSubM(0).setSwitchesList(SW_VERBOSE,true); // activate search for verbose switch
    forceSubM(0).setSwitchesList(SW_INTERPOLATION,true); // activate search for interpolate switch
    forceSubM(0).setSwitchesList(SW_SCALAR_VISCOSITY,true); // activate scalarViscosity switch
    forceSubM(0).readSwitches();

    particleCloud_.checkCG(true);
    if (propsDict_.found("scale") && typeCG_.size()==1)
    {
        // if "scale" is specified and there's only one single type, use "scale"
        scaleDia_=scalar(readScalar(propsDict_.lookup("scale")));
        typeCG_[0] = scaleDia_;
    }
    if (propsDict_.found("scaleDrag"))
    {
        scaleDrag_=scalar(readScalar(propsDict_.lookup("scaleDrag")));
    }

    if (typeCG_.size()>1) multiTypes_ = true;

    if (propsDict_.found("useFilteredDragModel"))
    {
        useGC_ = true;
        g_=propsDict_.lookupOrDefault<scalar>("g", 9.81);
        dPrim_=scalar(readScalar(propsDict_.lookup("dPrim")));
        rhoP_=scalar(readScalar(propsDict_.lookup("rhoP")));
        rho_=scalar(readScalar(propsDict_.lookup("rho")));
        nuf_=scalar(readScalar(propsDict_.lookup("nuf")));
	    scalar ut = terminalVelocity(1., dPrim_, nuf_, rho_, rhoP_, g_);
	    scalar Frp = ut*ut/g_/dPrim_;
        Lc2_ = ut*ut/g_*pow(Frp, -.6666667); // n is hardcoded as -2/3
        Info << "using grid coarsening correction with Lc2 = " << Lc2_ << " and ut = " << ut << " and Frp = " << Frp<< endl;

        if (propsDict_.found("useParcelSizeDependentFilteredDrag"))
        {
            usePC_ = true;
            if (propsDict_.found("k"))
                k_=scalar(readScalar(propsDict_.lookup("k")));
            Info << "using particle coarsening correction with k = " << k_ << endl;
        }
    }



}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

BeetstraDrag::~BeetstraDrag()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void BeetstraDrag::setForce() const
{
    if (typeCG_.size()>1 || typeCG_[0] > 1)
    {
        Info << "Beetstra using scale = " << typeCG_ << endl;
    }
    else if (particleCloud_.cg() > 1)
    {
        scaleDia_=particleCloud_.cg();
        Info << "Beetstra using scale from liggghts cg = " << scaleDia_ << endl;
    }

    const volScalarField& nufField = forceSubM(0).nuField();
    const volScalarField& rhoField = forceSubM(0).rhoField();

    vector position(0,0,0);
    scalar voidfraction(1);
    vector Ufluid(0,0,0);
    vector drag(0,0,0);
    label cellI=0;

    vector Us(0,0,0);
    vector Ur(0,0,0);
    scalar ds(0);
    scalar ds_scaled(0);
    scalar dSauter(0);
    scalar scaleDia3 = typeCG_[0]*typeCG_[0]*typeCG_[0];
    scalar nuf(0);
    scalar rho(0);
    scalar magUr(0);
    scalar Rep(0);
    scalar localPhiP(0);
    scalar GCcorr(1.);
    scalar PCcorr(1.);

    scalar cg = typeCG_[0];
    label partType = 1;

    vector dragExplicit(0,0,0);
    scalar dragCoefficient(0);

    interpolationCellPoint<scalar> voidfractionInterpolator_(voidfraction_);
    interpolationCellPoint<vector> UInterpolator_(U_);

    #include "setupProbeModel.H"

    for(int index = 0; index < particleCloud_.numberOfParticles(); ++index)
    {
            cellI = particleCloud_.cellIDs()[index][0];
            drag = vector::zero;
            dragExplicit = vector::zero;
            Ufluid = vector::zero;
            voidfraction=0;
            dragCoefficient = 0;

            if (cellI > -1) // particle found
            {

                if ( forceSubM(0).interpolation() )
                {
                    position     = particleCloud_.position(index);
                    voidfraction = voidfractionInterpolator_.interpolate(position,cellI);
                    Ufluid       = UInterpolator_.interpolate(position,cellI);
                    //Ensure interpolated void fraction to be meaningful
                    // Info << " --> voidfraction: " << voidfraction << endl;
                    if (voidfraction > 1.00) voidfraction = 1.0;
                    if (voidfraction < minVoidfraction_) voidfraction = minVoidfraction_;
                }
                else
                {
                    voidfraction = voidfraction_[cellI];
                    Ufluid = U_[cellI];
                }
                // in case a fines phase is present, void fraction needs to be adapted
                adaptVoidfraction(voidfraction, cellI);

                if (multiTypes_)
                {
                    partType = particleCloud_.particleType(index);
                    cg = typeCG_[partType - 1];
                    scaleDia3 = cg*cg*cg;
                }

                Us = particleCloud_.velocity(index);
                Ur = Ufluid-Us;
                magUr = mag(Ur);
                ds = 2*particleCloud_.radius(index);
                ds_scaled = ds/cg;
                dSauter = meanSauterDiameter(ds_scaled, cellI);

                rho = rhoField[cellI];
                nuf = nufField[cellI];

                localPhiP = 1.0f-voidfraction+SMALL;

                // calc particle's drag coefficient (i.e., Force per unit slip velocity and Stokes drag)

                Rep=dSauter*voidfraction*magUr/nuf + SMALL;

                dragCoefficient = F(voidfraction, Rep)
                                   *3*M_PI*nuf*rho*voidfraction
                                   *effDiameter(ds_scaled, cellI, index)
                                   *scaleDia3*scaleDrag_;

                // calculate filtering corrections
                if (useGC_)
                {
                    GCcorr = 1.-h(localPhiP)
                                /(  a(localPhiP)
                                    *Lc2_
                                    /std::cbrt(U_.mesh().V()[cellI])
                                  + 1.
                                 );
                    if (usePC_)
                    {
                        PCcorr = Foam::exp(k_*(1.-cg));
                    }
                }

                // apply filtering corrections
                dragCoefficient *= GCcorr*PCcorr;

                if (modelType_=="B")
                    dragCoefficient /= voidfraction;

                drag = dragCoefficient * Ur;

                // explicitCorr
                forceSubM(0).explicitCorr(drag,dragExplicit,dragCoefficient,Ufluid,U_[cellI],Us,UsField_[cellI],forceSubM(0).verbose());

                if(forceSubM(0).verbose() && index >=0 && index <2)
                {
                    Pout << "cellI = " << cellI << endl;
                    Pout << "index = " << index << endl;
                    Pout << "Us = " << Us << endl;
                    Pout << "Ur = " << Ur << endl;
                    Pout << "ds = " << ds << endl;
                    Pout << "ds/scale = " << ds/cg << endl;
                    Pout << "rho = " << rho << endl;
                    Pout << "nuf = " << nuf << endl;
                    Pout << "voidfraction = " << voidfraction << endl;
                    Pout << "Rep = " << Rep << endl;
                    Pout << "GCcorr = " << GCcorr << endl;
		            Pout << "PCcorr = " << PCcorr << endl;
                    Pout << "drag = " << drag << endl;
                }

                //Set value fields and write the probe
                if(probeIt_)
                {
                    #include "setupProbeModelfields.H"
                    vValues.append(drag);   //first entry must the be the force
                    vValues.append(Ur);
                    sValues.append(Rep);
                    sValues.append(voidfraction);
                    particleCloud_.probeM().writeProbe(index, sValues, vValues);
                }
            }

            // write particle based data to global array
            forceSubM(0).partToArray(index,drag,dragExplicit,Ufluid,dragCoefficient);
    }
}

/*********************************************************
 * "Drag Force of Intermediate Reynolds Number Flow Past *
 *  Mono- and Bidisperse Arrays of Spheres", eq. 16      *
 *  R Beetstra, M. A. van der Hoef, JAM Kuipers          *
 *  AIChE Journal 53(2) (2007)                           *
 *********************************************************/
double BeetstraDrag::F(double voidfraction, double Rep) const
{
    double localPhiP = max(SMALL,min(1.-SMALL,1.-voidfraction));
    return   10.0*localPhiP/(voidfraction*voidfraction)
           + voidfraction*voidfraction*(1.0+1.5*Foam::sqrt(localPhiP))
           + 0.413*Rep/(24*voidfraction*voidfraction)
                       *(1.0/voidfraction
                          +3*voidfraction*localPhiP
                          +8.4*Foam::pow(Rep,-0.343)
                        )
                        /(1+Foam::pow(10,3*localPhiP)
                            *Foam::pow(Rep,-0.5*(1+4*localPhiP))
                         );

}
 

/*********************************************************
 * "A drag model for filtered Euler-Lagange simulations  *
 *  of clustered gas-particle suspension",               *
 *  S. Radl, S. Sundaresan,                              *
 *  Chemical Engineering Science 117 (2014).             *
 *********************************************************/
double BeetstraDrag::a(double phiP) const
{
    double a0m = 0.;
    double a1m = 0.;
    double a2m = 0.;
    double a3m = 0.;
    double phipam = 0.;
    if      (phiP < 0.016)
    {
      a0m = 21.51;
    }
    else if (phiP < 0.100)
    {
      a0m =  1.96; a1m =  29.40; a2m =  164.91; a3m = -1923.;
    }
    else if (phiP < 0.180)
    {
      a0m =  4.63; a1m =   4.68; a2m = -412.04; a3m =  2254.; phipam = 0.10;
    }
    else if (phiP < 0.250)
    {
      a0m =  3.52; a1m = -17.99; a2m =  128.80; a3m =  -603.; phipam = 0.18;
    }
    else if (phiP < 0.400)
    {
      a0m =  2.68; a1m =  -8.20; a2m =    2.18; a3m = 112.33; phipam = 0.25;
    }
    else
    {
      a0m =  1.79;
    }
    const double phiP_minus_phipam = phiP-phipam;
    return a0m + a1m*phiP_minus_phipam + a2m*phiP_minus_phipam*phiP_minus_phipam
           + a3m*phiP_minus_phipam*phiP_minus_phipam*phiP_minus_phipam;
}

double BeetstraDrag::h(double phiP) const
{
  double h0m = 0.;
  double h1m = 0.;
  double h2m = 0.;
  double h3m = 0.;
  double phiphm = 0.;

  if      (phiP < 0.03)
  {
                 h1m = 7.97;
  }
  else if (phiP < 0.08)
  {
    h0m = 0.239; h1m =  4.640; h2m =   -4.410; h3m =  253.630; phiphm = 0.03;
  }
  else if (phiP < 0.12)
  {
    h0m = 0.492; h1m =  6.100; h2m =   33.630; h3m = -789.600; phiphm = 0.08;
  }
  else if (phiP < 0.18)
  {
    h0m = 0.739; h1m =  5.010; h2m =  -61.100; h3m =  310.800; phiphm = 0.12;
  }
  else if (phiP < 0.34)
  {
    h0m = 0.887; h1m =  1.030; h2m =   -5.170; h3m =    5.990; phiphm = 0.18;
  }
  else if (phiP < 0.48)
  {
    h0m = 0.943; h1m = -0.170; h2m =   -2.290; h3m =   -9.120; phiphm = 0.34;
  }
  else if (phiP < 0.55)
  {
    h0m = 0.850; h1m = -1.350; h2m =   -6.130; h3m = -132.600; phiphm = 0.48;
  }
  else
  {
    h0m = 0.680; h1m = -2.340; h2m = -225.200;                 phiphm = 0.55;
  }
  const double phiP_minus_phiphm = phiP-phiphm;
  return h0m + h1m*phiP_minus_phiphm + h2m*phiP_minus_phiphm*phiP_minus_phiphm
         + h3m*phiP_minus_phiphm*phiP_minus_phiphm*phiP_minus_phiphm;
}

double BeetstraDrag::terminalVelocity(double voidfraction, double dp, double nuf, double rhof, double rhop, double g) const
{
    scalar u0(dp*dp*fabs(rhof-rhop)*g/18./rhof/nuf);
    scalar Re(u0*dp/nuf);
    scalar res(1.);
    scalar u(u0);
    scalar Fi(0);
    scalar CdSt(0);
    Info << "uo: " << u0<<endl;
    int i = 0;
    while ((res > 1.e-6) && (i<100))
    {
        Info << "Iteration " << i;
        u0 = u;
        Info << ", u0 = " << u0;
        CdSt = 24/Re;
        Info << ", CdSt = " << CdSt;
        Fi = F(voidfraction, Re);
        Info << ", F = ";
        u = sqrt(1.333333333*fabs(rhof-rhop)*g*dp
                  /(CdSt*voidfraction*Fi*rhof)
            );
        Info << ", u = " << u;
        Re = fabs(u)*dp/nuf*voidfraction;
        res = fabs((u-u0)/u);
        Info << "Res: " << res << endl;
        i++;
    }
    if (res >1.e-6)
        FatalError << "Terminal velocity calculation diverged!" << endl;

    return u;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
