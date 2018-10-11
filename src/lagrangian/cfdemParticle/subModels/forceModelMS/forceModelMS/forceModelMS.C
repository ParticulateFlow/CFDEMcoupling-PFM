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
#include "forceModelMS.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(forceModelMS, 0);

defineRunTimeSelectionTable(forceModelMS, dictionary);

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
forceModelMS::forceModelMS
(
    const dictionary& dict,
    cfdemCloudMS& sm
)
:
    forceModel(dict,sm),
    particleCloudMS_(sm)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

forceModelMS::~forceModelMS()
{}

// * * * * * * * * * * * * * * * * Member Fct  * * * * * * * * * * * * * * * //
void forceModelMS::setForcesOnParticle() const
{
    if(!cloudRefMS().useforcePerClump())
    {
        int nrigidC(-1);
        label ind(-1);
        for(int index = 0;index <  cloudRefMS().numberOfParticles(); index++)
        {

            if (particleCloud_.cellIDs()[index][0] > -1) // particle Found
            {

                ind=cloudRefMS().body(index);
                if (ind < 0)
                {
                    Warning <<"clump was deleted??? ind = "<< ind << endl;
                }
                else
                {
                    nrigidC=cloudRefMS().nrigid(ind);

                    if (nrigidC <= 0)
                    {
                        Warning <<"A BUG occurred in GidaspowDragMS::setForce!!! nrigidC = " 
                                << nrigidC <<", ind = " << ind <<", index=" << index <<"\n" << endl;
                        nrigidC = 1000;
                    }
                    if(forceSubM(0).switches()[0]) for(int j=0;j<3;j++) particleCloud_.expForces()[index][j] += cloudRefMS().expForcesCM()[ind][j] / nrigidC;
                    else{
                        for(int j=0;j<3;j++){
                        particleCloud_.impForces()[index][j] += cloudRefMS().impForcesCM()[ind][j] / nrigidC;
                        }
                    }
                }
            }
        }
    }
}

cfdemCloudMS& forceModelMS::cloudRefMS() const
{
    return particleCloudMS_;
}

void forceModelMS::readDhbyV(dictionary& dict)
{
    if (dict.found("manDHdev"))
    {
        cloudRefMS().setManDHdev(Switch(dict.lookup("manDHdev")));
        cloudRefMS().setDHbyV(scalarList(dict.lookup("dHbyV")));
    }
    Warning << "You defined:" << cloudRefMS().dHbyV().size() 
            << " diameters manually - is this the number of clump types defined in *.in file?" << endl;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
