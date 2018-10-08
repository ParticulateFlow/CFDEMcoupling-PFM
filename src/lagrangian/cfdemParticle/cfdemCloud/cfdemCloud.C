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

#include "fileName.H"
#include "cfdemCloud.H"
#include "global.H"
#include "forceModel.H"
#include "locateModel.H"
#include "momCoupleModel.H"
#include "meshMotionModel.H"
#include "voidFractionModel.H"
#include "dataExchangeModel.H"
#include "IOModel.H"
#include "probeModel.H"
#include "averagingModel.H"
#include "clockModel.H"
#include "smoothingModel.H"
#include "liggghtsCommandModel.H"
#include "otherForceModel.H"

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
cfdemCloud::cfdemCloud
(
    const fvMesh& mesh
)
:
    mesh_(mesh),
    couplingProperties_
    (
        IOobject
        (
            "couplingProperties",
            mesh_.time().constant(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    liggghtsCommandDict_
    (
        IOobject
        (
            "liggghtsCommands",
            mesh_.time().constant(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    solveFlow_(couplingProperties_.lookupOrDefault<bool>("solveFlow", true)),
    verbose_(couplingProperties_.found("verbose")),
    ignore_(couplingProperties_.found("ignore")),
    allowCFDsubTimestep_(true),
    limitDEMForces_(couplingProperties_.found("limitDEMForces")),
    getParticleDensities_(couplingProperties_.lookupOrDefault<bool>("getParticleDensities",false)),
    getParticleEffVolFactors_(couplingProperties_.lookupOrDefault<bool>("getParticleEffVolFactors",false)),
    getParticleTypes_(couplingProperties_.lookupOrDefault<bool>("getParticleTypes",false)),
    maxDEMForce_(0.),
    modelType_(couplingProperties_.lookup("modelType")),
    positions_(NULL),
    velocities_(NULL),
    fluidVel_(NULL),
    fAcc_(NULL),
    impForces_(NULL),
    expForces_(NULL),
    DEMForces_(NULL),
    Cds_(NULL),
    radii_(NULL),
    voidfractions_(NULL),
    cellIDs_(NULL),
    particleDensities_(NULL),
    particleEffVolFactors_(NULL),
    particleTypes_(NULL),
    particleWeights_(NULL),
    particleVolumes_(NULL),
    particleV_(NULL),
    numberOfParticles_(0),
    d32_(-1),
    numberOfParticlesChanged_(false),
    arraysReallocated_(false),
    forceModels_(couplingProperties_.lookup("forceModels")),
    momCoupleModels_(couplingProperties_.lookup("momCoupleModels")),
    liggghtsCommandModelList_(liggghtsCommandDict_.lookup("liggghtsCommandModels")),
    otherForceModels_(couplingProperties_.lookupOrDefault<wordList>("otherForceModels",wordList(0))),
    turbulenceModelType_(couplingProperties_.lookup("turbulenceModelType")),
    cg_(1.),
    cgOK_(true),
    impDEMdrag_(false),
    impDEMdragAcc_(false),
    imExSplitFactor_(couplingProperties_.lookupOrDefault<scalar>("imExSplitFactor",1.0)),
    treatVoidCellsAsExplicitForce_(couplingProperties_.lookupOrDefault<bool>("treatVoidCellsAsExplicitForce",false)),
    useDDTvoidfraction_(couplingProperties_.found("useDDTvoidfraction")),
    ddtVoidfraction_
    (
        IOobject
        (
            "ddtVoidfraction",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimensionSet(0,0,-1,0,0), 0)  // 1/s
    ),
    particleDensityField_
    (
        IOobject
        (
            "partRho",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimensionSet(1,-3,0,0,0), 0.0)
    ),
    checkPeriodicCells_(false),
    turbulence_
    (
        mesh.lookupObject<turbulenceModel>
        (
            turbulenceModelType_
        )
    ),
    dataExchangeModel_
    (
        dataExchangeModel::New
        (
            couplingProperties_,
            *this
        )
    ),
    forceModel_(nrForceModels()),
    locateModel_
    (
        locateModel::New
        (
            couplingProperties_,
            *this
        )
    ),
    momCoupleModel_(nrMomCoupleModels()),
    IOModel_
    (
        IOModel::New
        (
            couplingProperties_,
            *this
        )
    ),
    probeModel_
    (
        probeModel::New
        (
            couplingProperties_,
            *this,
            "none",
            "none"
        )
    ),
    voidFractionModel_
    (
        voidFractionModel::New
        (
            couplingProperties_,
            *this
        )
    ),
    averagingModel_
    (
        averagingModel::New
        (
            couplingProperties_,
            *this
        )
    ),
    clockModel_
    (
        clockModel::New
        (
            couplingProperties_,
            mesh.time()
        )
    ),
    smoothingModel_
    (
        smoothingModel::New
        (
            couplingProperties_,
            *this
        )
    ),
    meshMotionModel_
    (
        meshMotionModel::New
        (
            couplingProperties_,
            *this
        )
    ),
    liggghtsCommand_(liggghtsCommandModelList_.size()),
    otherForceModel_(otherForceModels_.size())
{
    #include "versionInfo.H"
    global buildInfo(couplingProperties_,*this);
    buildInfo.info();

    Info << "If BC are important, please provide volScalarFields -imp/expParticleForces-" << endl;

    if(imExSplitFactor_ > 1.0)
        FatalError << "You have set imExSplitFactor > 1 in your couplingProperties. Must be <= 1."
                   << abort(FatalError);
    if(imExSplitFactor_ < 0.0)
        FatalError << "You have set imExSplitFactor < 0 in your couplingProperties. Must be >= 0."
                   << abort(FatalError);

    if (limitDEMForces_)
    {
        maxDEMForce_ = readScalar(couplingProperties_.lookup("limitDEMForces"));
    }

    if (turbulenceModelType_=="LESProperties")
    {
        Info << "WARNING - LES functionality not yet tested!" << endl;
    }

    if (!useDDTvoidfraction_)
    {
        Info << "ignoring ddt(voidfraction)" << endl;
    }

    const bool adjustTimeStep  = mesh_.time().controlDict().lookupOrDefault("adjustTimeStep", false);
    if (adjustTimeStep)
        FatalError << "CFDEMcoupling does not support adjustable time steps."
                   << abort(FatalError);

    forAll(momCoupleModel_, modeli)
    {
        momCoupleModel_.set
        (
            modeli,
            momCoupleModel::New
            (
                couplingProperties_,
                *this,
                momCoupleModels_[modeli]
            )
        );
    }

    forAll(forceModels_, modeli)
    {
        forceModel_.set
        (
            modeli,
            forceModel::New
            (
                couplingProperties_,
                *this,
                forceModels_[modeli]
            )
        );
    }

    // run liggghts commands from cfdem
    forAll(liggghtsCommand_, modeli)
    {
        liggghtsCommand_.set
        (
            modeli,
            liggghtsCommandModel::New
            (
                liggghtsCommandDict_,
                *this,
                liggghtsCommandModelList_[modeli],
                modeli
            )
        );
    }

    forAll(otherForceModel_, modeli)
    {
        otherForceModel_.set
        (
            modeli,
            otherForceModel::New
            (
                couplingProperties_,
                *this,
                otherForceModels_[modeli]
            )
        );
    }

    setCG(dataExchangeM().getCG());
    Switch cgWarnOnly(couplingProperties_.lookupOrDefault<Switch>("cgWarnOnly", true));

    if (!cgOK_ && cg_ > 1)
    {
        if (cgWarnOnly)
            Warning << "at least one of your models is not fit for cg !!!" << endl;
        else
            FatalError << "at least one of your models is not fit for cg !!!" << abort(FatalError);
    }

    // check if simulation is a fully periodic box
    const polyBoundaryMesh& patches = mesh_.boundaryMesh();
    int nPatchesCyclic(0);
    int nPatchesNonCyclic(0);
    forAll(patches, patchI)
    {
        const polyPatch& pp = patches[patchI];
        if (isA<cyclicPolyPatch>(pp) || isA<cyclicAMIPolyPatch>(pp))
            ++nPatchesCyclic;
        else if (!isA<processorPolyPatch>(pp))
            ++nPatchesNonCyclic;
    }

    if (nPatchesNonCyclic == 0)
    {
        checkPeriodicCells_ = true;
    }

    //hard set checkperiodic cells if wished
    if(this->couplingProperties().found("checkPeriodicCells"))
    {
        checkPeriodicCells_ = couplingProperties().lookupOrDefault<Switch>("checkPeriodicCells", checkPeriodicCells_);
    }

    if (nPatchesCyclic > 0 && nPatchesNonCyclic > 0)
    {
        if (verbose_) Info << "nPatchesNonCyclic=" << nPatchesNonCyclic << ", nPatchesCyclic=" << nPatchesCyclic << endl;
        Warning << "Periodic handing is disabled because the domain is not fully periodic!\n" << endl;
    }
}

// * * * * * * * * * * * * * * * * Destructors  * * * * * * * * * * * * * * //
cfdemCloud::~cfdemCloud()
{
    clockM().evalPar();
    clockM().normHist();
    dataExchangeM().destroy(positions_,3);
    dataExchangeM().destroy(velocities_,3);
    dataExchangeM().destroy(fluidVel_,3);
    dataExchangeM().destroy(fAcc_,3);
    dataExchangeM().destroy(impForces_,3);
    dataExchangeM().destroy(expForces_,3);
    dataExchangeM().destroy(DEMForces_,3);
    dataExchangeM().destroy(Cds_,1);
    dataExchangeM().destroy(radii_,1);
    dataExchangeM().destroy(voidfractions_,1);
    dataExchangeM().destroy(cellIDs_,1);
    dataExchangeM().destroy(particleWeights_,1);
    dataExchangeM().destroy(particleVolumes_,1);
    dataExchangeM().destroy(particleV_,1);
    if(getParticleDensities_) dataExchangeM().destroy(particleDensities_,1);
    if(getParticleEffVolFactors_) dataExchangeM().destroy(particleEffVolFactors_,1);
    if(getParticleTypes_) dataExchangeM().destroy(particleTypes_,1);
}

// * * * * * * * * * * * * * * * private Member Functions  * * * * * * * * * * * * * //
void cfdemCloud::getDEMdata()
{
    dataExchangeM().getData("radius","scalar-atom",radii_);
    dataExchangeM().getData("x","vector-atom",positions_);
    dataExchangeM().getData("v","vector-atom",velocities_);

    if(impDEMdragAcc_)
        dataExchangeM().getData("dragAcc","vector-atom",fAcc_); // array is used twice - might be necessary to clean it first

    if(getParticleDensities_) dataExchangeM().getData("density","scalar-atom",particleDensities_);
    if(getParticleEffVolFactors_) dataExchangeM().getData("effvolfactor","scalar-atom",particleEffVolFactors_);
    if(getParticleTypes_) dataExchangeM().getData("type","scalar-atom",particleTypes_);
}

void cfdemCloud::giveDEMdata()
{
    if(forceM(0).coupleForce())
    {
        dataExchangeM().giveData("dragforce","vector-atom",DEMForces_);

        if(impDEMdrag_)
        {
            if(verbose_) Info << "sending Ksl and uf" << endl;
            dataExchangeM().giveData("Ksl","scalar-atom",Cds_);
            dataExchangeM().giveData("uf","vector-atom",fluidVel_);
        }
    }
    if(verbose_) Info << "giveDEMdata done." << endl;
}

// * * *   write top level fields   * * * //

// * * * * * * * * * * * * * * * protected Member Functions  * * * * * * * * * * * * * //

void cfdemCloud::setNumberOfParticles(int nP)
{
    if(nP != numberOfParticles())
    {
        numberOfParticlesChanged_ = true;
        numberOfParticles_ = nP;
    }
}

void cfdemCloud::findCells()
{
    locateM().findCell(NULL,positions_,cellIDs_,numberOfParticles());
}

void cfdemCloud::setForces()
{
    resetArray(fluidVel_,numberOfParticles(),3);
    resetArray(impForces_,numberOfParticles(),3);
    resetArray(expForces_,numberOfParticles(),3);
    resetArray(DEMForces_,numberOfParticles(),3);
    resetArray(Cds_,numberOfParticles(),1);
    for (int i=0;i<cfdemCloud::nrForceModels();i++) cfdemCloud::forceM(i).setForce();

    if (limitDEMForces_)
    {
        scalar maxF = 0.0;
        for (int index = 0;index <  numberOfParticles(); ++index)
        {
            scalar F = mag(expForce(index));
            if (F > maxF) maxF = F;
            if (F > maxDEMForce_)
              for(int i=0;i<3;i++) DEMForces_[index][i] *= maxDEMForce_/F;
        }
        Info << "largest particle-fluid interaction on particle: " << maxF << endl;
    }
}

void cfdemCloud::setParticleForceField()
{
    averagingM().setVectorSum
    (
        forceM(0).impParticleForces(),
        impForces_,
        particleWeights_,
        NULL //mask
    );
    averagingM().setVectorSum
    (
        forceM(0).expParticleForces(),
        expForces_,
        particleWeights_,
        NULL //mask
    );
}

void cfdemCloud::setScalarAverages()
{
    if(!getParticleDensities_) return;
    if(verbose_) Info << "- setScalarAverage" << endl;

    particleDensityField_.primitiveFieldRef() = 0.0;
    averagingM().resetWeightFields();
    averagingM().setScalarAverage
    (
        particleDensityField_,
        particleDensities_,
        particleWeights_,
        averagingM().UsWeightField(),
        NULL
    );

    if(verbose_) Info << "setScalarAverage done." << endl;
}

void cfdemCloud::setVectorAverages()
{
    if(verbose_) Info << "- setVectorAverage(Us,velocities_,weights_)" << endl;
    averagingM().setVectorAverage
    (
        averagingM().UsNext(),
        velocities_,
        particleWeights_,
        averagingM().UsWeightField(),
        NULL //mask
    );
    if(verbose_) Info << "setVectorAverage done." << endl;
}

// * * * * * * * * * * * * * * * public Member Functions  * * * * * * * * * * * * * //

void cfdemCloud::checkCG(bool ok)
{
    if(!cgOK_) return;
    if(!ok) cgOK_ = ok;
}

void cfdemCloud::setPos(double**& pos)
{
    for(int index = 0; index <  numberOfParticles(); ++index)
    {
        for(int i=0; i<3; i++)
        {
            positions_[index][i] = pos[index][i];
        }
    }
}

// * * * * * * * * * * * * * * * ACCESS  * * * * * * * * * * * * * //

label cfdemCloud::particleCell(int index) const
{
    return cellIDs()[index][0];
}

vector cfdemCloud::position(int index) const
{
    return vector(positions()[index][0],positions()[index][1],positions()[index][2]);
}

vector cfdemCloud::velocity(int index) const
{
    return vector(velocities()[index][0],velocities()[index][1],velocities()[index][2]);
}

vector cfdemCloud::expForce(int index) const
{
    return vector(DEMForces()[index][0],DEMForces()[index][1],DEMForces()[index][2]);
}

vector cfdemCloud::fluidVel(int index) const
{
    return vector(fluidVels()[index][0],fluidVels()[index][1],fluidVels()[index][2]);
}

const forceModel& cfdemCloud::forceM(int i)
{
    return forceModel_[i];
}

label cfdemCloud::nrForceModels() const
{
    return forceModels_.size();
}

label cfdemCloud::nrMomCoupleModels() const
{
    return momCoupleModels_.size();
}

scalar cfdemCloud::voidfraction(int index) const
{
    return voidfractions()[index][0];
}

label cfdemCloud::liggghtsCommandModelIndex(const word& name) const
{
    return findIndex(liggghtsCommandModelList_, name);
}

// * * * * * * * * * * * * * * * WRITE  * * * * * * * * * * * * * //

// * * *   write cfdemCloud internal data   * * * //

bool cfdemCloud::evolve
(
    volScalarField& alpha,
    volVectorField& Us,
    volVectorField& U
)
{
    numberOfParticlesChanged_ = false;
    arraysReallocated_=false;
    bool doCouple=false;

    if(!ignore())
    {
        if (dataExchangeM().doCoupleNow())
        {
            Info << "\n Coupling..." << endl;
            dataExchangeM().couple(0);
            doCouple=true;

            // reset vol Fields
            clockM().start(16,"resetVolFields");
            if(verbose_)
            {
                Info << "couplingStep:" << dataExchangeM().couplingStep()
                     << "\n- resetVolFields()" << endl;
            }
            averagingM().resetVectorAverage(averagingM().UsPrev(),averagingM().UsNext(),false);
            voidFractionM().resetVoidFractions();
            averagingM().resetVectorAverage(forceM(0).impParticleForces(),forceM(0).impParticleForces(),true);
            averagingM().resetVectorAverage(forceM(0).expParticleForces(),forceM(0).expParticleForces(),true);
            averagingM().resetWeightFields();
            for (int i=0;i<momCoupleModels_.size(); i++)
                momCoupleM(i).resetMomSourceField();
            if(verbose_) Info << "resetVolFields done." << endl;
            clockM().stop("resetVolFields");

            if(verbose_) Info << "- getDEMdata()" << endl;
            clockM().start(17,"getDEMdata");
            getDEMdata();
            clockM().stop("getDEMdata");
            if(verbose_) Info << "- getDEMdata done." << endl;

            // search cellID of particles
            clockM().start(18,"findCell");
            if(verbose_) Info << "- findCell()" << endl;
            findCells();
            if(verbose_) Info << "findCell done." << endl;
            clockM().stop("findCell");

            // set void fraction field
            clockM().start(19,"setvoidFraction");
            if(verbose_) Info << "- setvoidFraction()" << endl;
            voidFractionM().setvoidFraction(NULL,voidfractions_,particleWeights_,particleVolumes_,particleV_);
            if(verbose_) Info << "setvoidFraction done." << endl;
            clockM().stop("setvoidFraction");

            // set average particles velocity field
            clockM().start(20,"setAverages");
            setScalarAverages();
            setVectorAverages();


            //Smoothen "next" fields
            smoothingM().dSmoothing();
            smoothingM().smoothen(voidFractionM().voidFractionNext());

            //only smoothen if we use implicit force coupling in cells void of particles
            //because we need unsmoothened Us field to detect cells for explicit
            //force coupling
            if(!treatVoidCellsAsExplicitForce())
                smoothingM().smoothenReferenceField(averagingM().UsNext());

            clockM().stop("setAverages");
        }

        //============================================
        //CHECK JUST TIME-INTERPOATE ALREADY SMOOTHENED VOIDFRACTIONNEXT AND UsNEXT FIELD
        //      IMPLICIT FORCE CONTRIBUTION AND SOLVER USE EXACTLY THE SAME AVERAGED
        //      QUANTITIES AT THE GRID!
        scalar timeStepFrac = dataExchangeM().timeStepFraction();
        Info << "\n timeStepFraction() = " << timeStepFrac << endl;
        if(timeStepFrac > 1.0000001)
        {
   //         FatalError << "cfdemCloud::dataExchangeM().timeStepFraction()>1: Do not do this, since dangerous. This might be due to the fact that you used a adjustable CFD time step. Please use a fixed CFD time step." << abort(FatalError);
              Warning << "cfdemCloud::dataExchangeM().timeStepFraction() = " << timeStepFrac << endl;
        }
        clockM().start(24,"interpolateEulerFields");

        // update voidFractionField
        alpha = voidFractionM().voidFractionInterp();
        if(dataExchangeM().couplingStep() < 2)
        {
            alpha.oldTime() = alpha; // supress volume src
            alpha.oldTime().correctBoundaryConditions();
        }
        alpha.correctBoundaryConditions();

        // calc ddt(voidfraction)
        calcDdtVoidfraction(alpha);

        // update mean particle velocity Field
        Us = averagingM().UsInterp();
        Us.correctBoundaryConditions();

        clockM().stop("interpolateEulerFields");
        //============================================

        if(doCouple)
        {
            // set particles forces
            clockM().start(21,"setForce");
            if(verbose_) Info << "- setForce(forces_)" << endl;
            setForces();
            if(verbose_) Info << "setForce done." << endl;
            clockM().stop("setForce");

            // get next force field
            clockM().start(22,"setParticleForceField");
            if(verbose_) Info << "- setParticleForceField()" << endl;
            setParticleForceField();
            if(verbose_) Info << "- setParticleForceField done." << endl;
            clockM().stop("setParticleForceField");

            // write DEM data
            if(verbose_) Info << " -giveDEMdata()" << endl;
            clockM().start(23,"giveDEMdata");
            giveDEMdata();
            clockM().stop("giveDEMdata");

            dataExchangeM().couple(1);
        }//end dataExchangeM().couple()


        if(verbose_){
            #include "debugInfo.H"
        }

        clockM().start(25,"dumpDEMdata");
        // do particle IO
        IOM().dumpDEMdata();
        clockM().stop("dumpDEMdata");

    }//end ignore
    return doCouple;
}

bool cfdemCloud::reAllocArrays()
{
    if(numberOfParticlesChanged_ && !arraysReallocated_)
    {
        // get arrays of new length
        dataExchangeM().allocateArray(positions_,0.,3);
        dataExchangeM().allocateArray(velocities_,0.,3);
        dataExchangeM().allocateArray(fluidVel_,0.,3);
        dataExchangeM().allocateArray(fAcc_,0.,3);
        dataExchangeM().allocateArray(impForces_,0.,3);
        dataExchangeM().allocateArray(expForces_,0.,3);
        dataExchangeM().allocateArray(DEMForces_,0.,3);
        dataExchangeM().allocateArray(Cds_,0.,1);
        dataExchangeM().allocateArray(radii_,0.,1);
        dataExchangeM().allocateArray(voidfractions_,1.,voidFractionM().maxCellsPerParticle());
        dataExchangeM().allocateArray(cellIDs_,-1,voidFractionM().maxCellsPerParticle());
        dataExchangeM().allocateArray(particleWeights_,0.,voidFractionM().maxCellsPerParticle());
        dataExchangeM().allocateArray(particleVolumes_,0.,voidFractionM().maxCellsPerParticle());
        dataExchangeM().allocateArray(particleV_,0.,1);
        if(getParticleDensities_) dataExchangeM().allocateArray(particleDensities_,0.,1);
        if(getParticleEffVolFactors_) dataExchangeM().allocateArray(particleEffVolFactors_,0.,1);
        if(getParticleTypes_) dataExchangeM().allocateArray(particleTypes_,0,1);
        arraysReallocated_ = true;
        return true;
    }
    return false;
}

tmp<fvVectorMatrix> cfdemCloud::divVoidfractionTau(volVectorField& U,volScalarField& voidfraction) const
{
    return
    (
      - fvm::laplacian(voidfractionNuEff(voidfraction), U)
      - fvc::div(voidfractionNuEff(voidfraction)*dev2(fvc::grad(U)().T()))
    );
}

tmp<volScalarField> cfdemCloud::ddtVoidfraction() const
{
    if (!useDDTvoidfraction_)
    {
        Info << "suppressing ddt(voidfraction)" << endl;
        return tmp<volScalarField> (ddtVoidfraction_ * 0.);
    }
    return tmp<volScalarField> (ddtVoidfraction_ * 1.) ;
}

void cfdemCloud::calcDdtVoidfraction(volScalarField& voidfraction)
{
    // version if ddt is calculated only at coupling time
    //Info << "calculating ddt(voidfraction) based on couplingTime" << endl;
    //scalar scale=mesh().time().deltaT().value()/dataExchangeM().couplingTime();
    //ddtVoidfraction_ = fvc::ddt(voidfraction) * scale;

    ddtVoidfraction_ = fvc::ddt(voidfraction);
}

/*tmp<fvVectorMatrix> cfdemCloud::ddtVoidfractionU(volVectorField& U,volScalarField& voidfraction) const
{
    if (dataExchangeM().couplingStep() <= 2) return fvm::ddt(U);

    return fvm::ddt(voidfraction,U);
}*/

tmp<volScalarField> cfdemCloud::voidfractionNuEff(volScalarField& voidfraction) const
{
    if (modelType_=="B" || modelType_=="Bfull")
    {
        return tmp<volScalarField>
        (
            #ifdef compre
                new volScalarField("viscousTerm", (turbulence_.mut() + turbulence_.mu()))
            #else
                new volScalarField("viscousTerm", (turbulence_.nut() + turbulence_.nu()))
            #endif
        );
    }
    else
    {
        return tmp<volScalarField>
        (
            #ifdef compre
                new volScalarField("viscousTerm", voidfraction*(turbulence_.mut() + turbulence_.mu()))
            #else
                new volScalarField("viscousTerm", voidfraction*(turbulence_.nut() + turbulence_.nu()))
            #endif
        );
    }
}

void cfdemCloud::resetArray(double**& array,int length,int width,double resetVal)
{
    for(int index = 0;index < length; ++index){
        for(int i=0;i<width;i++){
            array[index][i] = resetVal;
        }
    }
}

void cfdemCloud::otherForces(volVectorField& forcefield)
{
  forcefield.primitiveFieldRef() = vector::zero;
  forcefield.boundaryFieldRef() = vector::zero;
  for (int i=0;i<otherForceModels_.size();i++)
      forcefield += otherForceModel_[i].exportForceField();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
