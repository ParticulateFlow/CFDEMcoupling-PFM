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

    cloud class managing DEM data for CFD-DEM coupling

Class
    Foam::cfdemCloud

SourceFiles
    cfdemCloud.C
    cfdemCloudIO.C

\*---------------------------------------------------------------------------*/

#ifndef cfdemCloud_H
#define cfdemCloud_H

// choose version
#include "OFversion.H"
#include <vector>

#include "fvCFD.H"
#include "IFstream.H"

#include "turbulenceModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// forward declarations
class forceModel;
class locateModel;
class voidFractionModel;
class dataExchangeModel;
class IOModel;
class probeModel;
class averagingModel;
class clockModel;
class smoothingModel;
class momCoupleModel;
class meshMotionModel;
class liggghtsCommandModel;
class otherForceModel;

/*---------------------------------------------------------------------------*\
                           Class cfdemCloud Declaration
\*---------------------------------------------------------------------------*/

class cfdemCloud
{

// protected data
protected:
    const fvMesh& mesh_;

    IOdictionary couplingProperties_;

    IOdictionary liggghtsCommandDict_;

    const bool solveFlow_;

    const bool verbose_;

    const bool ignore_;

    bool allowCFDsubTimestep_;

    const bool limitDEMForces_;

    const bool getParticleDensities_;

    const bool getParticleEffVolFactors_;

    const bool getParticleTypes_;

    scalar maxDEMForce_;

    const word modelType_;

    double **positions_;

    double **velocities_;

    double **fluidVel_;

    double **fAcc_;

    double **impForces_;

    double **expForces_;

    double **DEMForces_;

    double **Cds_;

    double **radii_;

    double **voidfractions_;

    int **cellIDs_;

    double **particleDensities_;

    double **particleEffVolFactors_;

    int **particleTypes_;

    double **particleWeights_;

    double **particleVolumes_;

    double **particleV_;

    int numberOfParticles_;

    scalar d32_;

    bool numberOfParticlesChanged_;

    bool arraysReallocated_;

    const wordList forceModels_;

    const wordList momCoupleModels_;

    const wordList liggghtsCommandModelList_;

    const wordList otherForceModels_;

    const word turbulenceModelType_;

    scalar cg_;

    bool cgOK_;

    bool impDEMdrag_;

    bool impDEMdragAcc_;

    const scalar imExSplitFactor_;

    const bool treatVoidCellsAsExplicitForce_; //will treat the coupling force in cells with no Us data explicitly

    const bool useDDTvoidfraction_;

    volScalarField ddtVoidfraction_;

    volScalarField particleDensityField_;

    Switch checkPeriodicCells_;

    const turbulenceModel& turbulence_;

    autoPtr<dataExchangeModel> dataExchangeModel_;

    PtrList<forceModel> forceModel_;

    autoPtr<locateModel> locateModel_;

    PtrList<momCoupleModel> momCoupleModel_;

    autoPtr<IOModel> IOModel_;

    autoPtr<probeModel> probeModel_;

    autoPtr<voidFractionModel> voidFractionModel_;

    autoPtr<averagingModel> averagingModel_;

    autoPtr<clockModel> clockModel_;

    autoPtr<smoothingModel> smoothingModel_;

    autoPtr<meshMotionModel> meshMotionModel_;

    PtrList<liggghtsCommandModel> liggghtsCommand_;

    PtrList<otherForceModel> otherForceModel_;

// Protected member functions
    virtual void getDEMdata();

    virtual void giveDEMdata();

    virtual void setNumberOfParticles(int);

    virtual void findCells();

    virtual void setForces();

    virtual void setParticleForceField();

    virtual void setScalarAverages();

    virtual void setVectorAverages();

public:

    friend class dataExchangeModel;
    friend class voidFractionModel;
    friend class forceModel;
    friend class forceSubModel;

// Constructors

    //- Construct from mesh and a list of particles
    cfdemCloud
    (
            const fvMesh& mesh
    );

    //- Destructor
    virtual ~cfdemCloud();

// public Member Functions

    // Access
       bool allowCFDsubTimestep() const { return allowCFDsubTimestep_; }

       void setAllowCFDsubTimestep(bool b) { allowCFDsubTimestep_ = b; }

       void checkCG(bool);

       void setPos(double **&);

       const word& modelType() const { return modelType_; }

       label particleCell(int) const;

       vector position(int) const;

       vector velocity(int) const;

       vector expForce(int) const;

       vector fluidVel(int) const;

       virtual const forceModel& forceM(int);

       virtual label nrForceModels() const;

       virtual label nrMomCoupleModels() const;

       scalar voidfraction(int) const;

       label liggghtsCommandModelIndex(const word&) const;

       inline void setCG(double);

       inline scalar cg() const;

       inline bool impDEMdrag() const;

       inline bool impDEMdragAcc() const;

       inline scalar imExSplitFactor() const;

       inline bool treatVoidCellsAsExplicitForce() const;

       inline bool ignore() const;

       inline const fvMesh& mesh() const;

       inline bool solveFlow() const;

       inline bool verbose() const;

       inline const IOdictionary& couplingProperties() const;

       inline double ** positions() const;

       inline double ** velocities() const;

       inline double ** fluidVels() const;

       inline double ** fAccs() const;

       inline double ** impForces() const;

       inline double ** expForces() const;

       inline double ** DEMForces() const;

       inline double ** Cds() const;

       inline double ** radii() const;

       inline double ** voidfractions() const;

       inline void get_radii(double**&) const;

       inline int ** cellIDs() const;

       inline void get_cellIDs(int**&) const;

       inline double ** particleWeights() const;

       virtual inline label body(int) const;

       virtual inline double particleVolume(int) const;

       inline scalar radius(int) const;

       virtual inline double d(int) const;

       inline scalar d32(bool recalc=true);
       virtual inline double dMin() {return -1;}
       virtual inline double dMax() {return -1;}
       virtual inline int minType() {return -1;}
       virtual inline int maxType() {return -1;}
       virtual inline bool multipleTypesDMax() {return false;}
       virtual inline bool multipleTypesDMin() {return false;}

       inline bool getParticleDensities() const;
       virtual inline double ** particleDensities() const;
       virtual inline double particleDensity(label index) const;

       inline bool getParticleEffVolFactors() const;
       virtual inline double particleEffVolFactor(label index) const;

       inline bool getParticleTypes() const;
       virtual inline int ** particleTypes() const;
       virtual inline label particleType(label index) const;

       //access to the particle's rotation and torque data
       virtual inline double ** DEMTorques() const {return NULL;}
       virtual inline double ** omegaArray() const {return NULL;}
       virtual vector omega(int) const {return vector::zero;}

       //access to the particles' orientation information
       virtual inline double ** exArray() const {return NULL;}
       virtual vector ex(int) const {return vector::zero;}

       //Detector if SRF module is enable or not
       virtual inline bool SRFOn(){return false;}

       inline int numberOfParticles() const;

       inline bool numberOfParticlesChanged() const;

       inline int numberOfClumps() const;

       inline bool arraysReallocated() const;

       inline const wordList& forceModels() const;

       inline voidFractionModel& voidFractionM();
       inline const voidFractionModel& voidFractionM() const;

       inline const locateModel& locateM() const;

       inline momCoupleModel& momCoupleM(int);
       inline const momCoupleModel& momCoupleM(int) const;

       inline dataExchangeModel& dataExchangeM();

       inline const IOModel& IOM() const;

       inline probeModel& probeM();

       inline averagingModel& averagingM();
       inline const averagingModel& averagingM() const;

       inline const clockModel& clockM() const;

       inline const smoothingModel& smoothingM() const;

       inline const meshMotionModel& meshMotionM() const;

       inline const wordList& liggghtsCommandModelList() const;

       inline liggghtsCommandModel& liggghtsCommand(const label);

       inline const turbulenceModel& turbulence() const;

    // Write

      // write cfdemCloud internal data
        virtual bool evolve(volScalarField&,volVectorField&,volVectorField&);

        virtual void postFlow() {}

        virtual bool reAllocArrays();

    // IO
        void writeScalarFieldToTerminal(double**&) const;

        void writeVectorFieldToTerminal(double**&) const;

    // functions
        tmp<fvVectorMatrix> divVoidfractionTau(volVectorField& ,volScalarField&) const;

        tmp<volScalarField> ddtVoidfraction() const;

        void calcDdtVoidfraction(volScalarField& voidfraction);

        //tmp<fvVectorMatrix> ddtVoidfractionU(volVectorField& ,volScalarField&) const;

        tmp<volScalarField> voidfractionNuEff(volScalarField&) const;

        void resetArray(double**&,int,int,double resetVal=0.);

        void otherForces(volVectorField&);

        bool checkPeriodicCells() const { return checkPeriodicCells_; }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "cfdemCloudI.H"

#endif

// ************************************************************************* //
