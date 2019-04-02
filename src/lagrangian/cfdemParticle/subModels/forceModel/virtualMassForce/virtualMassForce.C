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

#include "virtualMassForce.H"
#include "addToRunTimeSelectionTable.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

#define NOTONCPU 9999

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(virtualMassForce, 0);

addToRunTimeSelectionTable
(
    forceModel,
    virtualMassForce,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
virtualMassForce::virtualMassForce
(
    const dictionary& dict,
    cfdemCloud& sm
)
:
    forceModel(dict,sm),
    propsDict_(dict.subDict(typeName + "Props")),
    velFieldName_(propsDict_.lookup("velFieldName")),
    U_(sm.mesh().lookupObject<volVectorField> (velFieldName_)),
	voidfractionFieldName_(propsDict_.lookup("voidfractionFieldName")),
	voidfraction_(sm.mesh().lookupObject<volScalarField> (voidfractionFieldName_)),
	UsFieldName_(propsDict_.lookup("granVelFieldName")),
	Us_(sm.mesh().lookupObject<volVectorField> (UsFieldName_)),
    phiFieldName_(propsDict_.lookup("phiFieldName")),
    phi_(sm.mesh().lookupObject<surfaceScalarField> (phiFieldName_)),
    UrelOld_(NULL),
	splitUrelCalculation_(propsDict_.lookupOrDefault<bool>("splitUrelCalculation",false)),
	useUs_(propsDict_.lookupOrDefault<bool>("useUs",false)),
	useFelderhof_(propsDict_.lookupOrDefault<bool>("useFelderhof",false)),
    Cadd_(0.5),
    DDtUrel_
    (   IOobject
        (
            "DDtUrel",
            sm.mesh().time().timeName(),
            sm.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        sm.mesh(),
        dimensionedVector("zero", dimensionSet(0,1,-2,0,0,0,0), vector(0,0,0))
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

    if (particleCloud_.dataExchangeM().maxNumberOfParticles() > 0)
    {
        // get memory for 2d array
        particleCloud_.dataExchangeM().allocateArray(UrelOld_,NOTONCPU,3);
    }

    // init force sub model
    setForceSubModels(propsDict_);
    // define switches which can be read from dict
    forceSubM(0).setSwitchesList(SW_TREAT_FORCE_EXPLICIT,true); // activate treatExplicit switch
    forceSubM(0).setSwitchesList(SW_INTERPOLATION,true); // activate search for interpolate switch
    forceSubM(0).readSwitches();

    //Extra switches/settings
    if(splitUrelCalculation_)
    {
        Info << "Virtual mass model: will split the Urel calculation\n";
        Info << "WARNING: be sure that LIGGGHTS integration takes ddtv_p implicitly into account! \n";
    }
    if(propsDict_.found("Cadd"))
    {
        Cadd_ = readScalar(propsDict_.lookup("Cadd"));
        Info << "Virtual mass model: using non-standard Cadd = " << Cadd_ << endl;
    }
    if(useUs_)
    {
        Info << "Virtual mass model: using averaged Us \n";
        Info << "WARNING: ignoring virtual mass of relative particle motion/collisions \n";
        
        if(splitUrelCalculation_)
        {
            FatalError << "Virtual mass model: useUs=true requires splitUrelCalculation_=false" << abort(FatalError);
        }
    }
    if(useFelderhof_)
    {
        Info << "Virtual mass model: using Cadd correlation by Felderhof \n";
        Info << "WARNING: ignoring user-set Cadd \n";

	    if(!dict.lookupOrDefault<bool>("getParticleDensities",false))
	    	FatalError << "Virtual mass model: useFelderhof=true requires getParticleDensities=true in couplingProperties" << abort(FatalError);
    }

    particleCloud_.checkCG(true);

    //Append the field names to be probed
    particleCloud_.probeM().initialize(typeName, typeName+".logDat");
    particleCloud_.probeM().vectorFields_.append("virtualMassForce"); //first entry must the be the force
    particleCloud_.probeM().vectorFields_.append("acceleration");
    //
    particleCloud_.probeM().vectorFields_.append("accelerationNoSmooth");
    //
    particleCloud_.probeM().scalarFields_.append("Cadd");
    particleCloud_.probeM().scalarFields_.append("Vs");
    particleCloud_.probeM().scalarFields_.append("rho");
    particleCloud_.probeM().scalarFields_.append("voidfraction");
    particleCloud_.probeM().scalarFields_.append("specificGravity");
    particleCloud_.probeM().writeHeader();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

virtualMassForce::~virtualMassForce()
{
    particleCloud_.dataExchangeM().destroy(UrelOld_,3);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void virtualMassForce::setForce() const
{
	if (!useUs_ || !splitUrelCalculation_ )
		reAllocArrays();

    scalar dt = U_.mesh().time().deltaT().value();

    //acceleration field
    if(splitUrelCalculation_)
        DDtUrel_ = fvc::ddt(U_) + fvc::div(phi_, U_); //Total Derivative of fluid velocity
    else if(useUs_)
        DDtUrel_ = fvc::ddt(U_) + fvc::div(phi_, U_) - fvc::ddt(Us_); //Total Derivative of fluid velocity minus average particle velocity

    //
        volVectorField DDtUrelNoSmooth_ = DDtUrel_;
        interpolationCellPoint<vector> DDtUrelNoSmoothInterpolator_(DDtUrelNoSmooth_);
    //
    
    // smoothen
    smoothingM().smoothen(DDtUrel_);

    interpolationCellPoint<vector> UInterpolator_(U_);
    interpolationCellPoint<vector> DDtUrelInterpolator_(DDtUrel_);
    interpolationCellPoint<scalar> voidfractionInterpolator_(voidfraction_);

    #include "setupProbeModel.H"

    bool haveUrelOld_(false);

    for(int index = 0;index <  particleCloud_.numberOfParticles(); index++)
    {
            vector virtualMassForce(0,0,0);
            vector position(0,0,0);
            vector DDtUrel(0,0,0);

            scalar voidfraction(1);
            scalar sg(1);

            label cellI = particleCloud_.cellIDs()[index][0];

            if (cellI > -1) // particle Found
            {
                // particle position
                if(forceSubM(0).interpolation())
                    position = particleCloud_.position(index);

                //********* acceleration value *********//
                if(splitUrelCalculation_ || useUs_ )
                {
                    // DDtUrel from acceleration field 
                    if(forceSubM(0).interpolation())
                		DDtUrel = DDtUrelInterpolator_.interpolate(position,cellI);
                	else
                		DDtUrel = DDtUrel_[cellI];
                }
                else
                {
                    // DDtUrel from UrelOld
                    vector Ufluid(0,0,0);

                    // relative velocity
                    if(forceSubM(0).interpolation())
                        Ufluid      = UInterpolator_.interpolate(position,cellI);
                    else
                        Ufluid = U_[cellI];

                    vector Us   = particleCloud_.velocity(index);
                    vector Urel = Ufluid - Us;

                    //Check of particle was on this CPU the last step
                	if(UrelOld_[index][0]==NOTONCPU) //use 1. element to indicate that particle was on this CPU the last time step
                		haveUrelOld_ = false;
                	else
                		haveUrelOld_ = true;

                	vector UrelOld(0,0,0);

                	for(int j=0;j<3;j++)
                	{
                		UrelOld[j]         = UrelOld_[index][j];
                		UrelOld_[index][j] = Urel[j];
                	}
                	if(haveUrelOld_ ) //only compute force if we have old (relative) velocity
                		DDtUrel = (Urel-UrelOld)/dt;
                }

                //********* Cadd value *********//
                scalar rho  = forceSubM(0).rhoField()[cellI];
                scalar ds = 2*particleCloud_.radius(index);
                scalar Vs = ds*ds*ds*M_PI/6;
                
                scalar Cadd;
                
                if (useFelderhof_)
                {                    
                    scalar epsilons(0);
                    scalar logsg(0);

                	if(forceSubM(0).interpolation())
		                voidfraction = voidfractionInterpolator_.interpolate(position,cellI);
		            else
		                voidfraction = voidfraction_[cellI];

                    sg       = particleCloud_.particleDensity(index) / rho;
                    logsg    = log(sg);
                    epsilons = 1-voidfraction;

                    Cadd = 0.5
                    		+ ( 0.047*logsg + 0.13)*epsilons
							+ (-0.066*logsg - 0.58)*epsilons*epsilons
							+ (               1.42)*epsilons*epsilons*epsilons;
                }
                else
                {
                    // use predefined value
                    Cadd = Cadd_;
                }

                //********* calculate force *********//
                virtualMassForce = Cadd_ * rho * Vs * DDtUrel;

                if( forceSubM(0).verbose() ) //&& index>100 && index < 105)
                {
                    Pout << "index / cellI = " << index << tab << cellI << endl;
                    Pout << "position = " << particleCloud_.position(index) << endl;
                }

                // Set value fields and write the probe
                if(probeIt_)
                {
                    //
                    vector DDtUrelNoSmooth(0,0,0);
                    if(forceSubM(0).interpolation())
                		DDtUrelNoSmooth = DDtUrelInterpolator_.interpolate(position,cellI);
                	else
                		DDtUrelNoSmooth = DDtUrelNoSmooth_[cellI];
                    //
                    #include "setupProbeModelfields.H"
                    vValues.append(virtualMassForce);           //first entry must the be the force
                    vValues.append(DDtUrel);
                    //
                    vValues.append(DDtUrelNoSmooth);
                    //
                    sValues.append(Cadd);
                    sValues.append(Vs);
                    sValues.append(rho);
                    sValues.append(voidfraction);
                    sValues.append(sg);
                    particleCloud_.probeM().writeProbe(index, sValues, vValues);
                }
            }
            else    //particle not on this CPU
                if (!useUs_ || !splitUrelCalculation_ )
                    UrelOld_[index][0] = NOTONCPU;

            // write particle based data to global array
            forceSubM(0).partToArray(index,virtualMassForce,vector::zero);
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
void Foam::virtualMassForce::reAllocArrays() const
{
    if(particleCloud_.numberOfParticlesChanged())
    {
        Pout << "virtualMassForce::reAllocArrays..." << endl;
        particleCloud_.dataExchangeM().allocateArray(UrelOld_,NOTONCPU,3);
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
