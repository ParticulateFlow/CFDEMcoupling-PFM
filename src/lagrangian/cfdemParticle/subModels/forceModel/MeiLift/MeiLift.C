/*---------------------------------------------------------------------------*\
    CFDEMcoupling - Open Source CFD-DEM coupling

    CFDEMcoupling is part of the CFDEMproject
    www.cfdem.com
                                Christoph Goniva, christoph.goniva@cfdem.com
                                Copyright 2009-2012 JKU Linz
                                Copyright 2012-2015 DCS Computing GmbH, Linz
                                Copyright 2015-     JKU Linz
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

    This function is based on the derivation in R. Mei,
    An approximate expression for shear lift force on a spherical  particle at a
    finite Reynolds number,
    Int. J. Multiph. Flow 18 (1992) 145–147

    The data for this functions is based on J.B. Mclaughlin,
    Inertial migration of a small sphere in linear shear flows,
    Journal of Fluid Mechanics. 224 (1991) 261-274.

    The second order terms are based on E. Loth and A. J. Dorgan,
    An equation of motion for particles of finite Reynolds number and size,
    Environ. Fluid Mech. 9 (2009) 187–206
    and can be added to the lift coefficient if desired
\*---------------------------------------------------------------------------*/

#include "error.H"

#include "MeiLift.H"
#include "addToRunTimeSelectionTable.H"

//#include <mpi.h>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(MeiLift, 0);

addToRunTimeSelectionTable
(
    forceModel,
    MeiLift,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
MeiLift::MeiLift
(
    const dictionary& dict,
    cfdemCloud& sm
)
:
    forceModel(dict,sm),
    propsDict_(dict.subDict(typeName + "Props")),
    velFieldName_(propsDict_.lookup("velFieldName")),
    U_(sm.mesh().lookupObject<volVectorField> (velFieldName_)),
	useShearInduced_(propsDict_.lookupOrDefault<bool>("useShearInduced",true)),
	useSpinInduced_(propsDict_.lookupOrDefault<bool>("useSpinInduced",false)),
	combineShearSpin_(propsDict_.lookupOrDefault<bool>("combineShearSpin",false))
{
	// read switches

	if(useShearInduced_)
	    Info << "Lift model: including shear-induced term.\n";

	if(useSpinInduced_)
    {
	    Info << "Lift model: including spin-induced term.\n";
	    Info << "Make sure to use a rolling friction model in LIGGGHTS!\n";
	    if(!dict.lookupOrDefault<bool>("getParticleAngVels",false))
	    	FatalError << "Lift model: useSpinInduced=true requires getParticleAngVels=true in couplingProperties" << abort(FatalError);
	}

	if(combineShearSpin_)
	{
	    Info << "Lift model: combining shear- and spin-terms by assuming equilibrium spin-rate.\n";
	    if(!useShearInduced_ || !useSpinInduced_)
	        FatalError << "Shear- and spin-induced lift must be activated in order to combine." << abort(FatalError);
	}

    // init force sub model
    setForceSubModels(propsDict_);
    // define switches which can be read from dict
    forceSubM(0).setSwitchesList(SW_TREAT_FORCE_EXPLICIT,true); // activate treatExplicit switch
    forceSubM(0).setSwitchesList(SW_VERBOSE,true); // activate search for verbose switch
    forceSubM(0).setSwitchesList(SW_INTERPOLATION,true); // activate search for interpolate switch
    forceSubM(0).setSwitchesList(SW_SCALAR_VISCOSITY,true); // activate scalarViscosity switch
    forceSubM(0).readSwitches();

    particleCloud_.checkCG(false);

    //Append the field names to be probed
    particleCloud_.probeM().initialize(typeName, typeName+".logDat");
    particleCloud_.probeM().vectorFields_.append("liftForce"); //first entry must the be the force
    particleCloud_.probeM().vectorFields_.append("Urel");      //other are debug
    particleCloud_.probeM().vectorFields_.append("vorticity");
    particleCloud_.probeM().vectorFields_.append("Ang_velocity");
    particleCloud_.probeM().scalarFields_.append("Rep");
    particleCloud_.probeM().scalarFields_.append("Rew");
    particleCloud_.probeM().scalarFields_.append("J*");
    particleCloud_.probeM().scalarFields_.append("Cl(shear)");
    particleCloud_.probeM().scalarFields_.append("Cl*(spin)");
    particleCloud_.probeM().scalarFields_.append("Omega_eq");
    particleCloud_.probeM().scalarFields_.append("Cl(combined)");
    particleCloud_.probeM().writeHeader();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

MeiLift::~MeiLift()
{
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void MeiLift::setForce() const
{
    const volScalarField& nufField = forceSubM(0).nuField();
    const volScalarField& rhoField = forceSubM(0).rhoField();

    // vectors
    vector position(0,0,0);
    vector lift(0,0,0);
    vector Us(0,0,0);
    vector Ur(0,0,0);
    vector Omega(0,0,0);
    vector vorticity(0,0,0);

    // properties
    scalar magUr(0);
    scalar magVorticity(0);
    scalar magOmega(0);
    scalar ds(0);
    scalar nuf(0);
    scalar rho(0);
    scalar voidfraction(1);

    // dimensionless groups
    scalar Rep(0);
    scalar Rew(0);

    // shear induced
    scalar omega_star(0);
    scalar Clshear(0);
    scalar J_star(0);
    scalar alphaStar(0);
    scalar epsilonSqr(0);
    scalar epsilon(0);

    // spin induced
    scalar Omega_star(0);
    scalar Clspin_star(0);

    // shear-spin combination
    scalar Omega_eq(0);
    scalar Clcombined(0);


    volVectorField vorticityField = fvc::curl(U_);

    interpolationCellPoint<vector> UInterpolator_(U_);
    interpolationCellPoint<vector> VorticityInterpolator_(vorticityField);

    #include "setupProbeModel.H"

    for (int index = 0; index < particleCloud_.numberOfParticles(); ++index)
    {
            lift = vector::zero;
            label cellI = particleCloud_.cellIDs()[index][0];

            if (cellI > -1) // particle Found
            {
                // properties
            	Us = particleCloud_.velocity(index);
            	Omega = particleCloud_.particleAngVel(index);

                if (forceSubM(0).interpolation())
                {
                    position  = particleCloud_.position(index);
                    Ur        = UInterpolator_.interpolate(position,cellI) - Us;
                    vorticity = VorticityInterpolator_.interpolate(position,cellI);
                }
                else
                {
                    Ur        = U_[cellI] - Us;
                    vorticity = vorticityField[cellI];
                }

                ds  	= 2. * particleCloud_.radius(index);
                nuf 	= nufField[cellI];
                rho 	= rhoField[cellI];

                magUr 	= mag(Ur);
                Rep 	= ds*magUr/nuf;

                // shear-induced lift
                if (useShearInduced_)
                {
                	magVorticity 	= mag(vorticity);
                	Rew 			= magVorticity*ds*ds/nuf;
                	omega_star 		= magVorticity * ds / magUr;
                    alphaStar 		= 0.5 * omega_star;
                    epsilonSqr 		= omega_star / Rep;
                    epsilon 		= sqrt(epsilonSqr);


                    //Basic model for the correction to the Saffman lift
                    //McLaughlin (1991), Mei (1992), Loth and Dorgan (2009)
                    //J_star = 0.443 * J
                    if (epsilon < 0.1) //epsilon << 1
                    {
                        //McLaughlin (1991), Eq (3.27): J = 32 * pi^2 * epsilon^5 * ln(1 / epsilon^2)
                        J_star = -140.0 * epsilonSqr * epsilonSqr * epsilon * log(1. / (epsilonSqr+SMALL));
                    }
                    else if (epsilon > 20.0) //epsilon >> 1
                    {
                        //McLaughlin (1991), Eq (3.26): J = 2.255 - 0.6463 / epsilon^2
                        J_star = 1.0 - 0.287 / epsilonSqr;
                    }
                    else
                    {
                        //Mei (1992), Eq (10)
                        //Loth and Dorgan (2009), Eq (32)
                        J_star = 0.3
                                 * (1.0   + tanh(2.5 * (log10(epsilon) + 0.191)))
                                 * (0.667 + tanh(6.0 * (      epsilon  - 0.32 )));
                    }

                    //Loth and Dorgan (2009), Eq (31), Eq (32)
                    Clshear = J_star * 4.11 * epsilon; //multiply correction to the basic Saffman model
                }

                if (useSpinInduced_)
                {
                	magOmega 		= mag(Omega);
                	Omega_star 		= magOmega * ds / magUr;

                	//Loth and Dorgan (2009), Eq (34)
                	Clspin_star		= 1.0 - (0.675 + 0.15 * (1.0 + tanh(0.28 * (Omega_star - 2.0)))) * tanh(0.18 * sqrt(Rep));
                }

                if (combineShearSpin_)
                {
                    //Loth and Dorgan (2009), Eq (38)
                    Omega_eq 		= alphaStar * (1.0 - 0.0075 * Rew) * (1.0 - 0.062 * sqrt(Rep) - 0.001 * Rep);
                    //Loth and Dorgan (2009), Eq (39)
                    Clcombined 		= Clshear + Clspin_star * Omega_eq;

                    if (magUr>0.0 && magVorticity>0.0)
                    {
                    	//Loth and Dorgan (2009), Eq (27)
                    	lift = 0.125 * constant::mathematical::pi
                    			* rho
								* Clcombined
								* magUr * magUr
								* (vorticity ^ Ur) / mag(vorticity ^ Ur) // force direction
								* ds * ds;
                    }
                }
                else
                {
                	//Loth and Dorgan (2009), Eq (36)
                	if (useShearInduced_)
                	{
                        if (magUr>0.0 && magVorticity>0.0)
                        {
                        	//Loth and Dorgan (2009), Eq (27)
                        	lift += 0.125 * constant::mathematical::pi
                        			* rho
									* Clshear
									* magUr * magUr
									* (vorticity ^ Ur) / mag(vorticity ^ Ur) // force direction
									* ds * ds;
                        }
                	}
                	if (useSpinInduced_)
                	{
                        if (magUr>0.0 && magOmega>0.0)
                        {
                        	//Loth and Dorgan (2009), Eq (33)
                        	lift += 0.125 * constant::mathematical::pi
                        			* rho
									* Clspin_star
									* (Omega ^ Ur)
									* ds * ds * ds;
                        }
                	}
                }

                if (modelType_ == "B")
                {
                    voidfraction = particleCloud_.voidfraction(index);
                    lift /= voidfraction;
                }


                //**********************************
                //SAMPLING AND VERBOSE OUTOUT
                if ( forceSubM(0).verbose() )
                {
                    Pout << "index = " << index << endl;
                    Pout << "Us = " << Us << endl;
                    Pout << "Ur = " << Ur << endl;
                    Pout << "vorticity = " << vorticity << endl;
                    Pout << "Omega = " << Omega << endl;
                    Pout << "ds = " << ds << endl;
                    Pout << "rho = " << rho << endl;
                    Pout << "nuf = " << nuf << endl;
                    Pout << "Rep = " << Rep << endl;
                    Pout << "Rew = " << Rew << endl;
                    Pout << "alphaStar = " << alphaStar << endl;
                    Pout << "epsilon = " << epsilon << endl;
                    Pout << "J_star = " << J_star << endl;
                    Pout << "Omega_eq = " << Omega_eq << endl;
                    Pout << "Clshear = " <<  Clshear<< endl;
                    Pout << "Clspin_star = " << Clspin_star << endl;
                    Pout << "Clcombined = " << Clcombined << endl;
                    Pout << "lift = " << lift << endl;
                }

                //Set value fields and write the probe
                if (probeIt_)
                {
                    #include "setupProbeModelfields.H"
                    vValues.append(lift);   //first entry must the be the force
                    vValues.append(Ur);
                    vValues.append(vorticity);
                    vValues.append(Omega);
                    sValues.append(Rep);
                    sValues.append(Rew);
                    sValues.append(J_star);
                    sValues.append(Clshear);
                    sValues.append(Clspin_star);
                    sValues.append(Omega_eq);
                    sValues.append(Clcombined);
                    particleCloud_.probeM().writeProbe(index, sValues, vValues);
                }
                // END OF SAMPLING AND VERBOSE OUTOUT
                //**********************************

            }
            // write particle based data to global array
            forceSubM(0).partToArray(index,lift,vector::zero);
    }

}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
