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

#include "DiFeliceDragMS.H"
#include "addToRunTimeSelectionTable.H"

//#include "mpi.h"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(DiFeliceDragMS, 0);

addToRunTimeSelectionTable
(
    forceModelMS,
    DiFeliceDragMS,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
DiFeliceDragMS::DiFeliceDragMS
(
    const dictionary& dict,
    cfdemCloudMS& sm
)
:
    forceModelMS(dict,sm),
    propsDict_(dict.subDict(typeName + "Props")),
    verbose_(false),
    velFieldName_(propsDict_.lookup("velFieldName")),
    U_(sm.mesh().lookupObject<volVectorField> (velFieldName_)),
    densityFieldName_(propsDict_.lookup("densityFieldName")),
    rho_(sm.mesh().lookupObject<volScalarField> (densityFieldName_)),
    voidfractionFieldName_(propsDict_.lookup("voidfractionFieldName")),
    voidfraction_(sm.mesh().lookupObject<volScalarField> (voidfractionFieldName_)),
    interpolation_(false),
    //sphereToClump_(readScalar(propsDict_.lookup("sphereToClump")))
    dH_(readScalar(propsDict_.lookup("hydraulicDiameter")))
{
    //Append the field names to be probed
    particleCloud_.probeM().initialize(typeName, "diFeliceDrag.logDat");
    particleCloud_.probeM().vectorFields_.append("dragForce"); //first entry must the be the force
    particleCloud_.probeM().vectorFields_.append("Urel");        //other are debug
    particleCloud_.probeM().scalarFields_.append("Rep");          //other are debug
    particleCloud_.probeM().scalarFields_.append("Cd");                 //other are debug
    particleCloud_.probeM().scalarFields_.append("voidfraction");       //other are debug
    particleCloud_.probeM().writeHeader();

    if (propsDict_.found("verbose")) verbose_=true;
    if (propsDict_.found("treatExplicit")) treatExplicit_=true;
    if (propsDict_.found("interpolation"))
    {
        Info << "using interpolated value of U." << endl;
        interpolation_=true;
    }
    particleCloud_.checkCG(false);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

DiFeliceDragMS::~DiFeliceDragMS()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void DiFeliceDragMS::setForce() const
{
    // get viscosity field
    #ifdef comp
        const volScalarField& nufField = cloudRefMS().turbulence().mu() / rho_;
    #else
        const volScalarField& nufField = cloudRefMS().turbulence().nu();
    #endif

    vector position(0,0,0);
    scalar voidfraction(1);
    vector Ufluid(0,0,0);
    vector drag(0,0,0);
    label cellI=0;
    vector Us(0,0,0);
    vector Ur(0,0,0);
    scalar ds(0);
    scalar nuf(0);
    scalar rho(0);
    scalar magUr(0);
    scalar Rep(0);
    scalar Cd(0);

    interpolationCellPoint<scalar> voidfractionInterpolator(voidfraction_);
    interpolationCellPoint<vector> UInterpolator(U_);

    #include "setupProbeModel.H"

    for(int index = 0;index <  cloudRefMS().numberOfClumps(); index++)
    {

        //if(mask[index][0])  // would have to be transformed from body ID to particle ID
        //{
            cellI = cloudRefMS().cellIDCM(index);
			drag = vector(0,0,0);

            if (cellI > -1) // particle Found
            {
                if(interpolation_)
                {
                    position = cloudRefMS().positionCM(index);
                    Ufluid = UInterpolator.interpolate(position,cellI);
                    voidfraction = voidfractionInterpolator.interpolate(position,cellI);
                }else
                {
                    Ufluid = U_[cellI];
                    voidfraction = voidfraction_[cellI]; //particleCloud_.voidfraction(index); ???//
                }

                Us = cloudRefMS().velocityCM(index);
                Ur = Ufluid-Us;
                //ds = 2*cloudRefMS().radius(index)/sphereToClump_; // scale from particle diameter
                ds = dH_;                                           // use dict defined diameter

                nuf = nufField[cellI];
                rho = rho_[cellI];
                magUr = mag(Ur);
                Rep = 0;
                Cd = 0;
                scalar phi(0);
                scalar phiN(0);

                if (magUr > 0)
                {
                    // calc particle Re Nr
                    Rep = ds*voidfraction*magUr/(nuf+SMALL);

                    // calc fluid drag Coeff
//                    phi=1; //AsurfAequi/Asurf;
//                    phiN=1; //AcrosssecAequi/Acrosssec;
                    // paper uses different Re definition!?
//                    Cd=8/(Rep*sqrt(phiN))+16/(Rep*sqrt(phi))+3/(sqrt(Rep)*pow(phi,0.75))+0.42*pow(10,(-0.4*pow(log(phi),0.2)))/phiN;

                    // calc fluid drag Coeff
                    Cd = sqr(0.63 + 4.8/sqrt(Rep));

                    // calc model coefficient Xi
                    scalar Xi = 3.7 - 0.65 * exp(-sqr(1.5-log10(Rep))/2);

                    // calc particle's drag
                    drag = 0.125*Cd*rho*M_PI*ds*ds*pow(voidfraction,(2-Xi))*magUr*Ur;

                    if (modelType_=="B")
                        drag /= voidfraction;
                }

                if(verbose_ && index >=0 && index <10)
                {
                    Info << "index = " << index << endl;
                    Info << "Us = " << Us << endl;
                    Info << "Ur = " << Ur << endl;
                    Info << "ds = " << ds << endl;
                    Info << "rho = " << rho << endl;
                    Info << "nuf = " << nuf << endl;
                    Info << "voidfraction = " << voidfraction << endl;
                    Info << "Rep = " << Rep << endl;
                    Info << "Cd = " << Cd << endl;
                    Info << "drag = " << drag << endl;
                }
            }
            // set force on bodies
            if(treatExplicit_) for(int j=0;j<3;j++) cloudRefMS().expForcesCM()[index][j] += drag[j];
            else  for(int j=0;j<3;j++) cloudRefMS().impForcesCM()[index][j] += drag[j];
            for(int j=0;j<3;j++) cloudRefMS().DEMForcesCM()[index][j] += drag[j];
        //}
    }

    // set force on particles
    setForcesOnParticle();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
