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

Class
    Foam::cfdemCloudMS derived from cfdemCloud

Description
    cloud class managing DEM data for CFD-DEM coupling including multi-sphere
    functionality

SourceFiles
    cfdemCloudMS.C


\*---------------------------------------------------------------------------*/

#ifndef cfdemCloudMS_H
#define cfdemCloudMS_H

#include "cfdemCloud.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class cfdemCloudMS Declaration
\*---------------------------------------------------------------------------*/

class cfdemCloudMS
:
    public cfdemCloud
{
private:

    mutable double **positionsCM_;
    mutable double **velocitiesCM_;
    mutable double **cellIDsCM_;
    mutable int **bodies_;
    mutable int **nrigids_;

    mutable double **exCM_;
    mutable double **eyCM_;
    mutable double **ezCM_;

    mutable double **typeCM_;
    mutable double **typeVolCM_;
    mutable double **VclumpCM_;
    mutable double **particleWeightsCM_;

    mutable double **dHCM_;

    //mutable double **SclumpCM_;
    //mutable double **scalingCM_;
    //mutable double **Cclump_ex_;
    //mutable double **Cclump_ey_;

    mutable double **impForcesCM_;
    mutable double **expForcesCM_;
    mutable double **DEMForcesCM_;

    int numberOfClumps_;
    scalar overlapCorr_;        // ratio between clump vol and sum(particle vol)
    Switch monoMS_;
    bool numberOfClumpsChanged_;

    bool useforcePerClump_;

    const wordList forceModels_;
    autoPtr<class forceModelMS>* forceModel_;

// Private member functions
    void getDEMdata();
    void giveDEMdata();
    bool reAllocArrays() const;
    void setNumberOfParticles(int);
    void findCells();
    void setForces();
    void setParticleForceField();
    void setClumpVolume(); // can be replaced once volume is communicated!!!
    void setdHCM();

public:

    // Constructors

        //- Construct from components
        cfdemCloudMS
        (
            const fvMesh& mesh
        );

    // Destructor

        virtual ~cfdemCloudMS();


    // Member Functions
      // Access

       vector positionCM(int);

       vector velocityCM(int);

       label cellIDCM(int);

       label body(int);

       label nrigid(int);

       inline double **& positionsCM() const;

       inline double **& velocitiesCM() const;

       inline double **& typeCM() const;

       inline double **& typeVolCM() const;

       inline double **& VclumpCM() const;

       inline double **& dHCM() const;

       inline double **& impForcesCM() const;

       inline double **& expForcesCM() const;

       inline double **& DEMForcesCM() const;

       inline int numberOfClumps() const;

       inline bool numberOfClumpsChanged() const;

       inline bool useforcePerClump() const;

       const forceModel& forceM(int);

       int nrForceModels();

       bool evolve(volScalarField&,volVectorField&,volVectorField&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "cfdemCloudMSI.H"

#endif

// ************************************************************************* //
