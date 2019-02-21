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

Contributing authors
    Paul Kieckhefen (TUHH)  2018-

\*---------------------------------------------------------------------------*/


#include "twoWayOne2One.H"
#include "addToRunTimeSelectionTable.H"
#include "clockModel.H"
#include "pair.h"
#include "force.h"
#include "forceModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(twoWayOne2One, 0);

addToRunTimeSelectionTable
(
    dataExchangeModel,
    twoWayOne2One,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
twoWayOne2One::twoWayOne2One
(
    const dictionary& dict,
    cfdemCloud& sm
)
:
    dataExchangeModel(dict,sm),
    propsDict_(dict.subDict(typeName + "Props")),
    thisLigPartner_(0),
    thisFoamPartner_(0),
    lig2foam_(nullptr),
    foam2lig_(nullptr),
    lig2foam_mask_(nullptr),
    lig2foam_ids_(nullptr),
    foam2lig_ids_(nullptr),
    lig2foam_vec_tmp_(nullptr),
    lig2foam_scl_tmp_(nullptr),
    foam2lig_vec_tmp_(nullptr),
    foam2lig_scl_tmp_(nullptr),
    staticProcMap_(propsDict_.lookupOrDefault<Switch>("useStaticProcMap", false)),
    cellIdComm_(propsDict_.lookupOrDefault<Switch>("useCellIdComm", false)),
    my_prev_cell_ids_fix_(nullptr),
    verbose_(propsDict_.lookupOrDefault("verbose", false)),
    lmp(nullptr)
{
    Info<<"Starting up LIGGGHTS for first time execution"<<endl;

    comm_liggghts_ = MPI_COMM_WORLD;

    // read path from dictionary
    const fileName liggghtsPath(propsDict_.lookup("liggghtsPath"));

    // open LIGGGHTS input script
    Info<<"Executing input script '"<< liggghtsPath.c_str() <<"'"<<endl;
    lmp = new LAMMPS_NS::LAMMPS(0,nullptr,comm_liggghts_);
    lmp->input->file(liggghtsPath.c_str());

    // get DEM time step size
    DEMts_ = lmp->update->dt;
    checkTSsize();

    // calculate boundingBox of FOAM subdomain
    primitivePatch tmpBoundaryFaces
    (
        SubList<face>
        (
            sm.mesh().faces(),
            sm.mesh().nFaces() - sm.mesh().nInternalFaces(),
            sm.mesh().nInternalFaces()
        ),
        sm.mesh().points()
    );
    typedef PrimitivePatch<face, List, const pointField> bPatch;
    bPatch boundaryFaces
    (
        tmpBoundaryFaces.localFaces(),
        tmpBoundaryFaces.localPoints()
    );
    thisFoamBox_ = treeBoundBox(boundaryFaces.localPoints());
    if (staticProcMap_)
    {
        createProcMap();
    }

    if (cellIdComm_)
    {
        my_prev_cell_ids_fix_ = static_cast<LAMMPS_NS::FixPropertyAtom*>
        (   lmp->modify->find_fix_property
            (
                "prev_cell_ids",
                "property/atom",
                "scalar",
                0,
                0,
                "cfd coupling",
                true
            )
        );
    }
}

void twoWayOne2One::createProcMap()
{
    List<treeBoundBox> foamBoxes(Pstream::nProcs());
    foamBoxes[Pstream::myProcNo()] = thisFoamBox_;
    Pstream::gatherList(foamBoxes);
    Pstream::scatterList(foamBoxes);

    // calculate bounding box of LIG subdomain
    // this may have to move to couple when dynamic LB occurs
    List<boundBox> ligBoxes(Pstream::nProcs());
    double** ligbb = o2o_liggghts_get_boundingbox(lmp);
    boundBox thisLigBox
    (
        point(ligbb[0][0], ligbb[0][1], ligbb[0][2]),
        point(ligbb[1][0], ligbb[1][1], ligbb[1][2])
    );
    ligBoxes[Pstream::myProcNo()] = thisLigBox;
    Pstream::gatherList(ligBoxes);
    Pstream::scatterList(ligBoxes);

    thisLigPartner_.clear();
    thisFoamPartner_.clear();

    // detect LIG subdomains which this FOAM has to interact with
    forAll(ligBoxes, ligproci)
    {
        if (thisFoamBox_.overlaps(ligBoxes[ligproci]))
        {
            thisLigPartner_.append(ligproci);
        }
    }
    // detect FOAM subdomains this LIG has to interact with
    // TODO: refactor to invert this list here
    forAll(foamBoxes, foamproci)
    {
        if (thisLigBox.overlaps(foamBoxes[foamproci]))
        {
            thisFoamPartner_.append(foamproci);
        }
    }

    if (verbose_)
    {
        Pout<< "FOAM bounding box: " << thisFoamBox_
            << " LIG bounding box: " << thisLigBox
            << nl
            << "FOAM comm partners: " << thisFoamPartner_
            << " LIG comm partners: " << thisLigPartner_
            << endl;
    }
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

twoWayOne2One::~twoWayOne2One()
{
    delete foam2lig_;
    delete lig2foam_;

    destroy(lig2foam_ids_);
    destroy(foam2lig_ids_);

    destroy(lig2foam_vec_tmp_);
    destroy(lig2foam_scl_tmp_);
    destroy(foam2lig_vec_tmp_);
    destroy(foam2lig_scl_tmp_);

    delete lmp;
}

// * * * * * * * * * * * * * * * public Member Functions  * * * * * * * * * * * * * //
void twoWayOne2One::getData
(
    word name,
    word type,
    double ** const& field,
    label /*step*/
) const
{
    if (name == "x") // the location is transferred by couple()
    {
        return;
    }
    if (type == "vector-atom")
    {
        double **tmp_= static_cast<double **>(lammps_extract_atom(lmp,name.c_str()));
        if (!tmp_)
        {
            LAMMPS_NS::Fix *fix = nullptr;
            fix = lmp->modify->find_fix_property
            (
                name.c_str(),
                "property/atom",
                "vector",
                0,
                0,
                "cfd coupling",
                false
            );
            if (fix)
            {
                tmp_ = static_cast<LAMMPS_NS::FixPropertyAtom*>(fix)->array_atom;
            }
            else
            {
                Warning<< "coupling fix not found!" << endl;
            }

            if (!tmp_)
            {
                FatalError<< "find_fix_property " << name
                          << " array_atom not found."
                          << abort(FatalError);
            }
        }

        lig2foam_->exchange<double>
        (
            tmp_,
            lig2foam_vec_tmp_,
            3
        );
        extractCollected<double>
        (
            lig2foam_vec_tmp_,
            const_cast<double**&>(field),
            3
        );
    }
    else if (type == "scalar-atom")
    {
        double *tmp_ = static_cast<double *>(lammps_extract_atom(lmp,name.c_str()));
        if (!tmp_)
        {
            LAMMPS_NS::Fix *fix = nullptr;
            fix = lmp->modify->find_fix_property
            (
                name.c_str(),
                "property/atom",
                "scalar",
                0,
                0,
                "cfd coupling",
                true
            );

            if (fix)
            {
                tmp_ = static_cast<LAMMPS_NS::FixPropertyAtom*>(fix)->vector_atom;
            }
            else
            {
                FatalError<< "coupling fix not found!" << abort(FatalError);
            }

            if (!tmp_)
            {
                FatalError<< "find_fix_property " << name
                          << " vector_atom not found."
                          << abort(FatalError);
            }
        }
        lig2foam_->exchange<double>
        (
            tmp_,
            lig2foam_scl_tmp_
        );
        extractCollected<double>
        (
            lig2foam_scl_tmp_,
            const_cast<double**&>(field)
        );
    }
    else
    {
        FatalError << "requesting type " << type << " and name " << name << abort(FatalError);
    }
}

void twoWayOne2One::getData
(
    word name,
    word type,
    int ** const& field,
    label /*step*/
) const
{
    FatalError << "do not use this getData!!!" << abort(FatalError);
/*
    o2o_data_liggghts_to_of
    (
        name.c_str(),
        type.c_str(),
        lmp,
        (void*&) field,
        "int",
        comm_liggghts_
    );
*/
}

void twoWayOne2One::giveData
(
    word name,
    word type,
    double ** const& field,
    const char* datatype
) const
{
    if (type == "vector-atom")
    {
        foam2lig_->exchange
        (
            const_cast<double**&>(field),
            foam2lig_vec_tmp_,
            3
        );
        o2o_data_of_to_liggghts
        (
            name.c_str(),
            type.c_str(),
            lmp,
            foam2lig_vec_tmp_,
            datatype,
            foam2lig_ids_,
            foam2lig_->ncollected_
        );
    }
    else if (type == "scalar-atom")
    {
        foam2lig_->exchange
        (
            const_cast<double**&>(field),
            foam2lig_scl_tmp_,
            1
        );
        o2o_data_of_to_liggghts
        (
            name.c_str(),
            type.c_str(),
            lmp,
            foam2lig_scl_tmp_,
            datatype,
            foam2lig_ids_,
            foam2lig_->ncollected_
        );
    }
    else
    {
        FatalError<< "twoWayMany2Many::giveData requested type " << type
                  << " not implemented!"
                  << abort(FatalError);
    }
}

//============
// double **
void twoWayOne2One::allocateArray
(
    double**& array,
    double initVal,
    int width,
    int length
) const
{
    int len = max(length,1);
    lmp->memory->grow(array, len, width, "o2o:dbl**");
    for (int i = 0; i < len; i++)
        for (int j = 0; j < width; j++)
            array[i][j] = initVal;
}

void twoWayOne2One::allocateArray
(
    double**& array,
    double initVal,
    int width,
    const char* length
) const
{
    int len = max(particleCloud_.numberOfParticles(),1);
    lmp->memory->grow(array, len, width, "o2o:dbl**:autolen");
    for (int i = 0; i < len; i++)
        for (int j = 0; j < width; j++)
            array[i][j] = initVal;
}

void inline twoWayOne2One::destroy(double** array,int len) const
{
    lmp->memory->destroy(array);
}

//============
// int **
void twoWayOne2One::allocateArray
(
    int**& array,
    int initVal,
    int width,
    int length
) const
{
    int len = max(length,1);
    lmp->memory->grow(array, len, width, "o2o:int**");
    for (int i = 0; i < len; i++)
        for (int j = 0; j < width; j++)
            array[i][j] = initVal;
}

void twoWayOne2One::allocateArray
(
    int**& array,
    int initVal,
    int width,
    const char* length
) const
{
    int len = max(particleCloud_.numberOfParticles(),1);
    lmp->memory->grow(array, len, width, "o2o:int**:autolen");
    for (int i = 0; i < len; i++)
        for (int j = 0; j < width; j++)
            array[i][j] = initVal;
}

void inline twoWayOne2One::destroy(int** array,int len) const
{
    lmp->memory->destroy(array);
}

//============
// double *
void twoWayOne2One::allocateArray(double*& array, double initVal, int length) const
{
    int len = max(length,1);
    lmp->memory->grow(array, len, "o2o:dbl*");
    for (int i = 0; i < len; i++)
        array[i] = initVal;
}

void inline twoWayOne2One::destroy(double* array) const
{
    lmp->memory->destroy(array);
}

//==============
// int *
void twoWayOne2One::allocateArray(int*& array, int initVal, int length) const
{
    int len = max(length,1);
    lmp->memory->grow(array, len, "o2o:int*");
    for (int i = 0; i < len; i++)
        array[i] = initVal;
}

void inline twoWayOne2One::destroy(int* array) const
{
    lmp->memory->destroy(array);
}
//==============

bool twoWayOne2One::couple(int i)
{
    bool coupleNow = false;
    if (i==0)
    {
        couplingStep_++;
        coupleNow = true;


        // run commands from liggghtsCommands dict
        Info<< "Starting up LIGGGHTS" << endl;
        particleCloud_.clockM().start(3,"LIGGGHTS");

        // check if liggghtsCommandModels with exaxt timing are being run
        bool exactTiming(false);
        int runComNr = -10;
        DynamicList<scalar> interruptTimes(0);
        DynamicList<int> DEMstepsToInterrupt(0);
        DynamicList<int> lcModel(0);

        forAll(particleCloud_.liggghtsCommandModelList(),i)
        {
            // Check if exact timing is needed
            // get time for execution
            // store time for execution in list
            if(particleCloud_.liggghtsCommand(i).exactTiming())
            {
                exactTiming = true;
                DynamicList<scalar> h
                  = particleCloud_.liggghtsCommand(i).executionsWithinPeriod
                (
                    TSstart(),
                    TSend()
                );

                forAll(h,j)
                {
                    // save interrupt times (is this necessary)
                    interruptTimes.append(h[j]);

                    // calc stepsToInterrupt
                    DEMstepsToInterrupt.append(DEMstepsTillT(h[j]));

                    // remember which liggghtsCommandModel to run
                    lcModel.append(i);
                }

                // make cumulative
                label len = DEMstepsToInterrupt.size();
                label ind(0);
                forAll(DEMstepsToInterrupt,i)
                {
                    ind = len - i - 1;
                    if(ind > 0)
                    {
                        DEMstepsToInterrupt[ind] -= DEMstepsToInterrupt[ind-1];
                    }
                }

                Info<< "Foam::twoWayOne2One::couple(i): interruptTimes="      << interruptTimes      << nl
                    << "Foam::twoWayOne2One::couple(i): DEMstepsToInterrupt=" << DEMstepsToInterrupt << nl
                    << "Foam::twoWayOne2One::couple(i): lcModel="             << lcModel
                    << endl;
            }

            if(particleCloud_.liggghtsCommand(i).type() == "runLiggghts")
            {
                runComNr = i;
            }
        }

        // models with exact timing exists
        label commandLines(0);
        if(exactTiming)
        {
            // extension for more liggghtsCommands active the same time:
            //    sort interrupt list within this run period
            //    keep track of corresponding liggghtsCommand
            int DEMstepsRun(0);

            forAll(interruptTimes,j)
            {
                // set run command till interrupt
                DEMstepsRun += DEMstepsToInterrupt[j];
                particleCloud_.liggghtsCommand(runComNr).set(DEMstepsToInterrupt[j]);
                const char* command = particleCloud_.liggghtsCommand(runComNr).command(0);
                Info<< "Executing run command: '"<< command <<"'"<< endl;
                lmp->input->one(command);

                // run liggghts command with exact timing
                command = particleCloud_.liggghtsCommand(lcModel[j]).command(0);
                Info << "Executing command: '"<< command <<"'"<< endl;
                lmp->input->one(command);
            }

            // do the run
            if(particleCloud_.liggghtsCommand(runComNr).runCommand(couplingStep()))
            {
                particleCloud_.liggghtsCommand(runComNr).set(couplingInterval() - DEMstepsRun);
                const char* command = particleCloud_.liggghtsCommand(runComNr).command(0);
                Info<< "Executing run command: '"<< command <<"'"<< endl;
                lmp->input->one(command);
            }

            // do the other non exact timing models
            forAll(particleCloud_.liggghtsCommandModelList(),i)
            {
                if
                (
                  ! particleCloud_.liggghtsCommand(i).exactTiming() &&
                    particleCloud_.liggghtsCommand(i).runCommand(couplingStep())
                )
                {
                    commandLines=particleCloud_.liggghtsCommand(i).commandLines();
                    for(int j=0;j<commandLines;j++)
                    {
                        const char* command = particleCloud_.liggghtsCommand(i).command(j);
                        Info << "Executing command: '"<< command <<"'"<< endl;
                        lmp->input->one(command);
                    }
                }
            }
        }
        // no model with exact timing exists
        else
        {
            forAll(particleCloud_.liggghtsCommandModelList(),i)
            {
                if(particleCloud_.liggghtsCommand(i).runCommand(couplingStep()))
                {
                    commandLines=particleCloud_.liggghtsCommand(i).commandLines();
                    for(int j=0;j<commandLines;j++)
                    {
                        const char* command = particleCloud_.liggghtsCommand(i).command(j);
                        Info << "Executing command: '"<< command <<"'"<< endl;
                        lmp->input->one(command);
                    }
                }
            }
        }

        particleCloud_.clockM().stop("LIGGGHTS");
        Info<< "LIGGGHTS finished" << endl;

        if (!staticProcMap_)
        {
            createProcMap();
        }

        setupLig2FoamCommunication();

        locateParticles();

        setupFoam2LigCommunication();

        if (verbose_)
        {
            Pout<< "FOAM owns " << getNumberOfParticles()
                << " LIG owns " << lmp->atom->nlocal
                << nl
                << "FOAM collects " << lig2foam_->ncollected_
                << " LIG collects " << foam2lig_->ncollected_
                << endl;
        }
    }

    return coupleNow;
}

void twoWayOne2One::setupLig2FoamCommunication()
{
    int* src_procs = new int[thisLigPartner_.size()];
    for (int proci = 0; proci < thisLigPartner_.size(); proci++)
    {
        src_procs[proci] = thisLigPartner_[proci];
    }
    int* dst_procs = new int[thisFoamPartner_.size()];
    for (int proci = 0; proci < thisFoamPartner_.size(); proci++)
    {
        dst_procs[proci] = thisFoamPartner_[proci];
    }

    delete lig2foam_;
    lig2foam_ = new One2One(comm_liggghts_);
    lig2foam_->setup
    (
      thisLigPartner_.size(),
      src_procs,
      thisFoamPartner_.size(),
      dst_procs,
      lmp->atom->nlocal
    );
    allocateArray
    (
        lig2foam_vec_tmp_,
        0.,
        3 * lig2foam_->ncollected_
    );
    allocateArray
    (
        lig2foam_scl_tmp_,
        0.,
        lig2foam_->ncollected_
    );
}


void twoWayOne2One::locateParticles()
{
    // get positions for locate
    double** my_positions = static_cast<double**>(lmp->atom->x);
    double*  my_flattened_positions = nullptr;
    allocateArray(my_flattened_positions, 0., 3*lmp->atom->nlocal);
    for (int atomi = 0; atomi < lmp->atom->nlocal; atomi++)
    {
        for (int coordi = 0; coordi < 3; coordi++)
        {
            my_flattened_positions[atomi*3+coordi] = my_positions[atomi][coordi];
        }
    }

    double* collected_flattened_positions = nullptr;
    allocateArray(collected_flattened_positions, 0., 3*lig2foam_->ncollected_);

    lig2foam_->exchange(my_flattened_positions, collected_flattened_positions, 3);
    destroy(my_flattened_positions);

    double* my_prev_cell_ids = nullptr;
    double* prev_cell_ids = nullptr;
    if (cellIdComm_)
    {
        my_prev_cell_ids = my_prev_cell_ids_fix_->vector_atom;
        allocateArray(prev_cell_ids, -1, lig2foam_->ncollected_);
        lig2foam_->exchange(my_prev_cell_ids, prev_cell_ids);
    }

    if (lig2foam_mask_)
    {
        delete [] lig2foam_mask_;
    }
    lig2foam_mask_ = new bool[lig2foam_->ncollected_];

    DynamicList<label> cellIds;
    cellIds.setCapacity(lig2foam_->ncollected_);
    label n_located(0);
    label roundedCelli(-1);
    const label nCells(particleCloud_.mesh().cells().size());
    for (int atomi = 0; atomi < lig2foam_->ncollected_; atomi++)
    {
        const vector position = vector
        (
            collected_flattened_positions[3*atomi+0],
            collected_flattened_positions[3*atomi+1],
            collected_flattened_positions[3*atomi+2]
        );
        if (!thisFoamBox_.contains(position))
        {
            lig2foam_mask_[atomi] = false;
            continue;
        }
        const label cellI = particleCloud_.locateM().findSingleCell
        (
            position,
            cellIdComm_
            ?   // don't know whether using round is efficient
                (roundedCelli = round(prev_cell_ids[atomi])) < nCells
                ?
                roundedCelli
                :
                -1
            :
            -1
        );

        lig2foam_mask_[atomi] = false;
        if (cellI >= 0) // in domain
        {
            lig2foam_mask_[atomi] = true;
            n_located++;
            cellIds.append(cellI);
        }
    }
    if (cellIdComm_)
    {
        destroy(prev_cell_ids);
    }

    setNumberOfParticles(n_located);
    particleCloud_.reAllocArrays();

    reduce(n_located, sumOp<label>());
    if (verbose_ || n_located != returnReduce(lmp->atom->nlocal, sumOp<label>()))
    {
        Warning << "Have located " << n_located
                << " ouf of " << returnReduce(lmp->atom->nlocal, sumOp<label>())
                << " particles in FOAM. "
                << endl;
    }

    // copy positions/cellids/ids of located particles into arrays
    allocateArray(lig2foam_ids_, 0, getNumberOfParticles());
    int* collected_ids = nullptr;
    allocateArray(collected_ids, 0, lig2foam_->ncollected_);
    lig2foam_->exchange<int>(lmp->atom->tag, collected_ids);
    extractCollected<int>(collected_ids, lig2foam_ids_);
    destroy(collected_ids);

    double* extracted_flattened_positions = new double[getNumberOfParticles()*3];
    extractCollected<double>
    (
        collected_flattened_positions,
        extracted_flattened_positions,
        3
    );
    setPositions(getNumberOfParticles(), extracted_flattened_positions);
    destroy(extracted_flattened_positions);
    destroy(collected_flattened_positions);

    setCellIDs(cellIds);
}

void twoWayOne2One::setupFoam2LigCommunication()
{
    int* src_procs = new int[thisFoamPartner_.size()];
    for (int proci = 0; proci < thisFoamPartner_.size(); proci++)
    {
        src_procs[proci] = thisFoamPartner_[proci];
    }

    int* dst_procs = new int[thisLigPartner_.size()];
    for (int proci = 0; proci < thisLigPartner_.size(); proci++)
    {
        dst_procs[proci] = thisLigPartner_[proci];
    }

    delete foam2lig_;
    foam2lig_ = new One2One(comm_liggghts_);

    foam2lig_->setup
    (
        thisFoamPartner_.size(),
        src_procs,
        thisLigPartner_.size(),
        dst_procs,
        getNumberOfParticles()
    );
    allocateArray
    (
        foam2lig_ids_,
        0,
        foam2lig_->ncollected_
    );
    foam2lig_->exchange<int>(lig2foam_ids_, foam2lig_ids_);

    allocateArray
    (
        foam2lig_vec_tmp_,
        0.,
        3 * foam2lig_->ncollected_
    );
    allocateArray
    (
        foam2lig_scl_tmp_,
        0.,
        foam2lig_->ncollected_
    );

    if (cellIdComm_)
    {
        double** dbl_cell_ids = new double*[getNumberOfParticles()];
        for (int atomi = 0; atomi < getNumberOfParticles(); atomi++)
        {   // TEMPORARY: if this persists after 19.07.2018, call me.
            dbl_cell_ids[atomi] = new double[1];
            dbl_cell_ids[atomi][0] = particleCloud_.cellIDs()[atomi][0];
        }
        giveData("prev_cell_ids", "scalar-atom", dbl_cell_ids, "double");
        delete [] dbl_cell_ids;
    }
}

template <typename T>
void twoWayOne2One::extractCollected(T**& src, T**& dst, int width) const
{
    int locali = 0;

    for (int atomi = 0; atomi < lig2foam_->ncollected_; atomi++)
    {
        if (!lig2foam_mask_[atomi]) continue;

        for (int coordi = 0; coordi < width; coordi++)
        {
            dst[locali][coordi] = src[atomi][coordi];
        }
        locali++;
    }
}

template <typename T>
void twoWayOne2One::extractCollected(T*& src, T*& dst, int width) const
{
    int locali = 0;

    for (int atomi = 0; atomi < lig2foam_->ncollected_; atomi++)
    {
        if (!lig2foam_mask_[atomi]) continue;

        for (int coordi = 0; coordi < width; coordi++)
        {
            dst[locali] = src[atomi*width+coordi];
            locali++;
        }
    }
}

template <typename T>
void twoWayOne2One::extractCollected(T*& src, T**& dst, int width) const
{
    int locali = 0;

    for (int atomi = 0; atomi < lig2foam_->ncollected_; atomi++)
    {
        if (!lig2foam_mask_[atomi]) continue;

        for (int coordi = 0; coordi < width; coordi++)
        {
            dst[locali][coordi] = src[atomi*width+coordi];
        }
        locali++;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
