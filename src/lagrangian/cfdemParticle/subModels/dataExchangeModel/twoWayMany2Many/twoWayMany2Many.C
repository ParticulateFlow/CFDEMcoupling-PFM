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
#include "twoWayMany2Many.H"
#include "addToRunTimeSelectionTable.H"
#include "clockModel.H"
#include "memory.h"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(twoWayMany2Many, 0);

addToRunTimeSelectionTable
(
    dataExchangeModel,
    twoWayMany2Many,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
twoWayMany2Many::twoWayMany2Many
(
    const dictionary& dict,
    cfdemCloud& sm
)
:
    dataExchangeModel(dict,sm),
    propsDict_(dict.subDict(typeName + "Props")),
    pbm_(sm.mesh().boundaryMesh()),
    pData_(sm.mesh().globalData()),
    procPatches_(pData_.processorPatches()),
    procPatchIndices_(pData_.processorPatchIndices()),
    neighbourProcs_(pData_[Pstream::myProcNo()]),
    neighbourProcIndices_(Pstream::nProcs(), -1)
{
    allowDiagComm_=true;
    if (propsDict_.found("allowDiagComm"))
        allowDiagComm_=Switch(propsDict_.lookup("allowDiagComm"));
    if(!allowDiagComm_)
        Warning << "Make sure you decompose only in one direction as allowDiagComm flag is false!" << endl;

    forAll(neighbourProcs_, i) neighbourProcIndices_[neighbourProcs_[i]] = i;

    Info<<"Starting up LIGGGHTS for first time execution"<<endl;

    MPI_Comm_rank(MPI_COMM_WORLD,&me);
    MPI_Comm_size(MPI_COMM_WORLD,&nprocs);

    if (me < nprocs) liggghts = 1;
    else liggghts = MPI_UNDEFINED;

    MPI_Comm_split(MPI_COMM_WORLD,liggghts,0,&comm_liggghts);

    // open LIGGGHTS input script
    char *liggghtsPathChar = new char[256];
    int n = 0;
    if (me == 0)
    {
      // read path from dictionary
      const fileName liggghtsPath(propsDict_.lookup("liggghtsPath"));
      strcpy(liggghtsPathChar, liggghtsPath.c_str());
      n = strlen(liggghtsPathChar) + 1;

      Info<<"Executing input script '"<< liggghtsPath.c_str() <<"'"<<endl;
    }

    if (liggghts == 1) lmp = new LAMMPS_NS::LAMMPS(0,NULL,comm_liggghts);

    MPI_Bcast(&n,1,MPI_INT,0,MPI_COMM_WORLD);
    if (n > 0)
    {
        MPI_Bcast(liggghtsPathChar,n,MPI_CHAR,0,MPI_COMM_WORLD);
        if (liggghts == 1) lmp->input->file(liggghtsPathChar);
    }

    delete [] liggghtsPathChar;

    // get DEM time step size
    DEMts_ = lmp->update->dt;
    checkTSsize();

    // m2m stuff
    firstRun_ = true;
    particleLost_ = false;
    Npart_ = -1;
    lmp2foam_ = NULL;
    lmp2foam_vec_ = NULL;
    foam2lmp_vec_ = NULL;
    foam2lmp_ = NULL;
    nlocal_lammps_ = -1;
    id_lammps_ = NULL;
    id_lammpsVec_ = NULL;
    nlocal_foam_ = -1;
    id_foam_ = NULL;
    id_foamVec_ = NULL;
    tmp_ = NULL;
    tmpI_ = NULL;
    pos_lammps_ = NULL;
    nlocal_foam_lost_ = -1;
    id_foamLost_ = NULL;
    id_foamLostAll = NULL;
    lost_pos_ = NULL;
    lost_posAll = NULL;
    cellID_foam_ = NULL;
    pos_foam_ = NULL;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

twoWayMany2Many::~twoWayMany2Many()
{
    free(id_lammps_);
    free(id_lammpsVec_);
    free(id_foamVec_);
    free(id_foam_);
    free(id_foamLost_);
    free(pos_foam_);
    free(cellID_foam_);
    delete[] lost_posAll;
    free(lost_pos_);
    delete[] id_foamLostAll;
    destroy(tmpI_);
    destroy(tmp_);
    destroy(pos_lammps_,3);
    delete lmp2foam_;
    delete lmp2foam_vec_;
    delete foam2lmp_vec_;
    delete foam2lmp_;
    delete lmp;
}


// * * * * * * * * * * * * * * * private Member Functions  * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * public Member Functions  * * * * * * * * * * * * * //
void twoWayMany2Many::getData
(
    word name,
    word type,
    double ** const& field,
    label /*step*/
) const
{
    if (name != "x")
    {
        if (type == "vector-atom")
        {
            double **tmp_ = static_cast<double **>(lammps_extract_atom(lmp,name.c_str()));
            if (!tmp_)
            {
                LAMMPS_NS::Fix *fix = NULL;
                fix = lmp->modify->find_fix_property(name.c_str(),"property/atom","vector",0,0,"cfd coupling",false);
                if (fix)
                    tmp_ = static_cast<LAMMPS_NS::FixPropertyAtom*>(fix)->array_atom;
                else
                    Warning << "coupling fix not found!" << endl;

                if (!tmp_)
                    FatalError << "find_fix_property " << name << " array_atom not found." << abort(FatalError);
            }

            lmp2foam_vec_->exchange(tmp_ ? tmp_[0] : NULL, field[0]);
        }
        else if (type == "scalar-atom")
        {
            double *tmp_ = static_cast<double *>(lammps_extract_atom(lmp,name.c_str()));
            if (!tmp_)
            {
                LAMMPS_NS::Fix *fix = NULL;
                fix = lmp->modify->find_fix_property(name.c_str(),"property/atom","scalar",0,0,"cfd coupling",true);
                if (fix)
                    tmp_ = static_cast<LAMMPS_NS::FixPropertyAtom*>(fix)->vector_atom;
                else
                    FatalError << "coupling fix not found!" << abort(FatalError);

                if (!tmp_)
                    FatalError << "find_fix_property " << name << " vector_atom not found." << abort(FatalError);
            }

            lmp2foam_->exchange(tmp_, field[0]);
        }
        else
        {
            FatalError << "requesting type " << type << " and name " << name << abort(FatalError);
        }
    }
}

void twoWayMany2Many::getData
(
    word name,
    word type,
    int ** const& field,
    label /*step*/
) const
{
    data_liggghts_to_of(name.c_str(), type.c_str(), lmp, (void*&)field, "int");
}

void twoWayMany2Many::giveData
(
    word name,
    word type,
    double ** const& field,
    const char* /*datatype*/
) const
{
    if (type == "vector-atom")
    {
        double **tmp_ = NULL;
        LAMMPS_NS::Fix *fix = NULL;
        fix = lmp->modify->find_fix_property(name.c_str(),"property/atom","vector",0,0,"cfd coupling",false);
        if (fix)
            tmp_ = static_cast<LAMMPS_NS::FixPropertyAtom*>(fix)->array_atom;
        else
            Warning << "coupling fix not found!"<<endl;

        if (!tmp_)
            FatalError << "find_fix_property " << name << " array_atom not found." << abort(FatalError);

        foam2lmp_vec_->exchange(field[0],tmp_ ? tmp_[0] : NULL);
    }
    else if (type == "scalar-atom")
    {
        double *tmp_ = NULL;
        LAMMPS_NS::Fix *fix = NULL;
        fix = lmp->modify->find_fix_property(name.c_str(),"property/atom","scalar",0,0,"cfd coupling",false);
        if (fix)
            tmp_ = static_cast<LAMMPS_NS::FixPropertyAtom*>(fix)->vector_atom;
        else
            FatalError << "coupling fix not found!" << abort(FatalError);

        if (!tmp_)
            FatalError << "find_fix_property " << name << " vector_atom not found." << abort(FatalError);

        foam2lmp_->exchange(field[0],tmp_); // for double *
    }
    else
    {
        FatalError << "twoWayMany2Many::giveData requested type " << type << " not implemented! \n" << abort(FatalError);
    }
}

//============
// double **
void twoWayMany2Many::allocateArray
(
    double**& array,
    double initVal,
    int width,
    int length
) const
{
    int len = max(length,1);
    lmp->memory->grow(array, len, width, "m2m:dbl**");
    for (int i = 0; i < len; i++)
        for (int j = 0; j < width; j++)
            array[i][j] = initVal;
}

void twoWayMany2Many::allocateArray
(
    double**& array,
    double initVal,
    int width,
    const char* length
) const
{
    int len = max(particleCloud_.numberOfParticles(),1);
    lmp->memory->grow(array, len, width, "m2m:dbl**:autolen");
    for (int i = 0; i < len; i++)
        for (int j = 0; j < width; j++)
            array[i][j] = initVal;
}

void inline twoWayMany2Many::destroy(double** array,int len) const
{
    lmp->memory->destroy(array);
}

//============
// int **
void twoWayMany2Many::allocateArray
(
    int**& array,
    int initVal,
    int width,
    int length
) const
{
    int len = max(length,1);
    lmp->memory->grow(array, len, width, "m2m:int**");
    for (int i = 0; i < len; i++)
        for (int j = 0; j < width; j++)
            array[i][j] = initVal;
}

void twoWayMany2Many::allocateArray
(
    int**& array,
    int initVal,
    int width,
    const char* length
) const
{
    int len = max(particleCloud_.numberOfParticles(),1);
    lmp->memory->grow(array, len, width, "m2m:int**:autolen");
    for (int i = 0; i < len; i++)
        for (int j = 0; j < width; j++)
            array[i][j] = initVal;
}

void inline twoWayMany2Many::destroy(int** array,int len) const
{
    lmp->memory->destroy(array);
}

//============
// double *
void twoWayMany2Many::allocateArray(double*& array, double initVal, int length) const
{
    int len = max(length,1);
    lmp->memory->grow(array, len, "m2m:dbl*");
    for (int i = 0; i < len; i++)
        array[i] = initVal;
}

void inline twoWayMany2Many::destroy(double* array) const
{
    lmp->memory->destroy(array);
}

//==============
// int *
void twoWayMany2Many::allocateArray(int*& array, int initVal, int length) const
{
    int len = max(length,1);
    lmp->memory->grow(array, len, "m2m:int*");
    for (int i = 0; i < len; i++)
        array[i] = initVal;
}

void inline twoWayMany2Many::destroy(int* array) const
{
    lmp->memory->destroy(array);
}
//==============


bool twoWayMany2Many::couple(int i)
{
    bool coupleNow = false;
    if (i==0)
    {
        couplingStep_++;
        coupleNow = true;

        // start liggghts
        if (liggghts == 1)
        {
            // run commands from liggghtsCommands dict
            Info<<"Starting up LIGGGHTS" << endl;
            particleCloud_.clockM().start(3,"LIGGGHTS");
            forAll(particleCloud_.liggghtsCommandModelList(),i)
            {

                if (particleCloud_.liggghtsCommand(i).runCommand(couplingStep()))
                {
                    label commandLines = particleCloud_.liggghtsCommand(i).commandLines();
                    for (int j=0; j<commandLines; j++)
                    {
                        const char* command = particleCloud_.liggghtsCommand(i).command(j);
                        Info << "Executing command: '" << command << "'" << endl;
                        lmp->input->one(command);
                    }
                }
            }
            particleCloud_.clockM().stop("LIGGGHTS");
            Info << "LIGGGHTS finished" << endl;
        }

        double newNpart = liggghts_get_maxtag(lmp);
        setNumberOfParticles(newNpart);

        if (Npart_ != newNpart)
        {
            Npart_ = newNpart;
            firstRun_ = true;
        }

        particleCloud_.clockM().start(4,"CoupleSyncIDs()");
        syncIDs();
        firstRun_=false;
        particleCloud_.clockM().stop("CoupleSyncIDs()");

        setNumberOfParticles(nlocal_foam_);

        // re-allocate arrays of cloud
        particleCloud_.reAllocArrays();

        setPositions(nlocal_foam_,pos_foam_);
        setCellIDs(nlocal_foam_,cellID_foam_);

        Info <<"Foam::twoWayMany2Many::couple(i) done." << endl;
    }
    return coupleNow;
}

int twoWayMany2Many::getNumberOfParticles() const
{
    return liggghts_get_maxtag(lmp);
}

int twoWayMany2Many::getNumberOfClumps() const
{
    Warning << "Foam::twoWayMany2Many::getNumberOfClumps() - changes necessary here" << endl;
    //return liggghts_get_maxtag_ms(lmp);
    return 1;
}

void twoWayMany2Many::syncIDs() const
{
    particleCloud_.clockM().start(5,"recv_DEM_ids");

    delete lmp2foam_;
    delete lmp2foam_vec_;
    delete foam2lmp_vec_;
    delete foam2lmp_;
    lmp2foam_ = new Many2Many(MPI_COMM_WORLD);
    lmp2foam_vec_ = new Many2Many(MPI_COMM_WORLD);
    foam2lmp_vec_ = new Many2Many(MPI_COMM_WORLD);
    foam2lmp_ = new Many2Many(MPI_COMM_WORLD);

    nlocal_lammps_ = *(static_cast<int *>(lammps_extract_global(lmp,"nlocal")));

    int*  id_lammpsSync = NULL;
    double** pos_lammpsSync = NULL;
    bool pos_lammps_alloc_flag = false;
    bool id_lammps_alloc_flag = false;

    if (firstRun_ || particleLost_)
    {
        id_lammps_ = NULL;
        id_lammps_ = static_cast<int *>(lammps_extract_atom(lmp,"id"));

        allocateArray(id_lammpsVec_,0,nlocal_lammps_*3);
        allocateArray(id_lammpsSync,0,nlocal_lammps_);
        for (int i = 0; i < nlocal_lammps_; i++)
        {
            id_lammpsSync[i]=id_lammps_[i];

            for (int j=0;j<3;j++)
                id_lammpsVec_[i*3+j] = id_lammps_[i]*3+j;
        }
        destroy(pos_lammps_,0);
        pos_lammps_ = NULL;
        pos_lammps_ = static_cast<double **>(lammps_extract_atom(lmp,"x"));
        pos_lammps_alloc_flag = false;
        id_lammps_alloc_flag = false;
    }
    else
    {
        id_lammps_ = static_cast<int *>(lammps_extract_atom(lmp,"id"));
        allocateArray(id_lammpsSync,0,nlocal_lammps_);
        for (int i = 0; i < nlocal_lammps_; i++)
            id_lammpsSync[i] = id_lammps_[i];

        allocateArray(id_lammpsVec_,0,nlocal_lammps_*3);
        for (int i = 0; i < nlocal_lammps_; i++)
            for (int j=0;j<3;j++)
                id_lammpsVec_[i*3+j] = id_lammpsSync[i]*3+j;

        lmp2foam_->setup(nlocal_lammps_,id_lammpsSync,nlocal_foam_,id_foam_);
        lmp2foam_vec_->setup(nlocal_lammps_*3,id_lammpsVec_,nlocal_foam_*3,id_foamVec_);
        foam2lmp_vec_->setup(nlocal_foam_*3,id_foamVec_,nlocal_lammps_*3,id_lammpsVec_);
        foam2lmp_->setup(nlocal_foam_,id_foam_,nlocal_lammps_,id_lammpsSync);

        id_lammps_ = NULL;
        allocateArray(id_lammps_,-1.,nlocal_foam_);
        id_lammps_alloc_flag = true;

        allocateArray(tmpI_,-1.,nlocal_foam_);
        lmp2foam_->exchange(id_lammpsSync, tmpI_);
        for (int i=0;i<nlocal_foam_;i++)
            id_lammps_[i] = tmpI_[i];

        pos_lammpsSync = static_cast<double **>(lammps_extract_atom(lmp,"x"));

        allocateArray(tmp_,-1.,nlocal_foam_*3);
        lmp2foam_vec_->exchange(pos_lammpsSync ? pos_lammpsSync[0] : NULL, tmp_);

        allocateArray(pos_lammps_,0,3,nlocal_foam_);
        pos_lammps_alloc_flag = true;
        for (int i=0;i<nlocal_foam_;i++)
            for (int j=0;j<3;j++)
                pos_lammps_[i][j] = tmp_[i*3+j];

    }
    particleCloud_.clockM().stop("recv_DEM_ids");

    particleCloud_.clockM().start(6,"locateParticle()");
    locateParticle(id_lammpsSync, id_lammps_alloc_flag);
    id_lammps_alloc_flag = true;
    particleCloud_.clockM().stop("locateParticle()");

    particleCloud_.clockM().start(11,"setup_Comm");

    delete lmp2foam_;
    delete lmp2foam_vec_;
    delete foam2lmp_vec_;
    delete foam2lmp_;
    lmp2foam_ = new Many2Many(MPI_COMM_WORLD);
    lmp2foam_vec_ = new Many2Many(MPI_COMM_WORLD);
    foam2lmp_vec_ = new Many2Many(MPI_COMM_WORLD);
    foam2lmp_ = new Many2Many(MPI_COMM_WORLD);

    if (firstRun_ || particleLost_)
    {
        lmp2foam_->setup(nlocal_lammps_,id_lammps_,nlocal_foam_,id_foam_);
        lmp2foam_vec_->setup(nlocal_lammps_*3,id_lammpsVec_,nlocal_foam_*3,id_foamVec_);
        foam2lmp_vec_->setup(nlocal_foam_*3,id_foamVec_,nlocal_lammps_*3,id_lammpsVec_);
        foam2lmp_->setup(nlocal_foam_,id_foam_,nlocal_lammps_,id_lammps_);
    }
    else
    {
        lmp2foam_->setup(nlocal_lammps_,id_lammpsSync,nlocal_foam_,id_foam_);
        lmp2foam_vec_->setup(nlocal_lammps_*3,id_lammpsVec_,nlocal_foam_*3,id_foamVec_);
        foam2lmp_vec_->setup(nlocal_foam_*3,id_foamVec_,nlocal_lammps_*3,id_lammpsVec_);
        foam2lmp_->setup(nlocal_foam_,id_foam_,nlocal_lammps_,id_lammpsSync);
    }
    if (id_lammps_alloc_flag) destroy(id_lammps_);
    id_lammps_ = NULL;    // free pointer from LIG
    destroy(id_lammpsSync);
    id_lammpsSync = NULL; // free pointer from LIG
    if (pos_lammps_alloc_flag) destroy(pos_lammps_,0);
    pos_lammps_ = NULL; // free pointer from LIG

    particleCloud_.clockM().stop("setup_Comm");
}

void twoWayMany2Many::locateParticle(int* id_lammpsSync, bool id_lammps_alloc_flag) const
{
#if defined(version21)

    int nop = particleCloud_.numberOfParticles();

    allocateArray(id_foamLost_,0,nop);
    allocateArray(lost_pos_,0.,nop*3);
    allocateArray(id_foam_,0,nop);
    allocateArray(id_foamVec_,0,nop*3);
    allocateArray(cellID_foam_,0,nop);
    if(firstRun_)
        allocateArray(pos_foam_,0,nop*3);

    particleCloud_.clockM().start(7,"locate_Stage1");
    int iterate;
    if (firstRun_ || particleLost_) iterate=nlocal_lammps_;
    else iterate = nlocal_foam_;

    nlocal_foam_ = 0;
    nlocal_foam_lost_ = 0;
    vector pos;
    label cellID = 0;
    label searchCellID = -1;
    List< DynamicList<int> > particleTransferID(neighbourProcs_.size());
    List< DynamicList<vector> > particleTransferPos(neighbourProcs_.size());

    for (int i = 0; i < iterate; i++)
    {
        pos = vector(pos_lammps_[i][0],pos_lammps_[i][1],pos_lammps_[i][2]);
        cellID = particleCloud_.locateM().findSingleCell(pos,searchCellID);
        point oldPos(pos_foam_[nlocal_foam_*3+0],pos_foam_[nlocal_foam_*3+1],pos_foam_[nlocal_foam_*3+2]);

        if (cellID >= 0)
        {
            id_foam_[nlocal_foam_] = id_lammps_[i];

            for (int j=0;j<3;j++)
            {
                id_foamVec_[nlocal_foam_*3+j] = id_lammps_[i]*3+j;
                pos_foam_[nlocal_foam_*3+j] = pos[j];
            }
            cellID_foam_[nlocal_foam_] = cellID;

            nlocal_foam_ += 1;
        }
        else
        {
            bool commPart=false;
            point newPos=pos;
            label nearestFace = particleCloud_.locateM().intersection(oldPos,newPos);

            if (nearestFace >= particleCloud_.mesh().nInternalFaces())
            {
                label patchI = pbm_.whichPatch(nearestFace);

                label n(-1);
                if (procPatchIndices_[patchI] != -1)
                {
                    n = neighbourProcIndices_
                    [
                        refCast<const processorPolyPatch>
                        (
                            pbm_[patchI]
                        ).neighbProcNo()
                    ];

                    if (n==Pstream::myProcNo())
                    {
                        //Pout << couplingStep_ << "st communicating particle " << id_lammps_[i]
                        //     << "communication fails as particle travels diagonal or jumps over proc" << endl;
                    }
                    else
                    {
                        particleTransferID[n].append(id_lammps_[i]);
                        particleTransferPos[n].append(pos);
                        commPart=true;
                       //Pout << couplingStep_ << "st communicating particle " << id_lammps_[i] << ", to proc# " << n << endl;
                    }
                }
            }
            if (!commPart)
            {
                id_foamLost_[nlocal_foam_lost_] = id_lammps_[i];

                for (int j=0; j<3; j++)
                    lost_pos_[nlocal_foam_lost_*3+j] = pos[j];

                nlocal_foam_lost_ += 1;
            }
        }
    }
    particleCloud_.clockM().stop("locate_Stage1");
    particleCloud_.clockM().start(8,"locate_Stage2");

    PstreamBuffers pBufs(Pstream::nonBlocking);

    forAll(particleTransferID, i)
    {
        if (particleTransferID[i].size())
        {
            UOPstream particleStream
            (
                neighbourProcs_[i],
                pBufs
            );
            particleStream << particleTransferID[i]<<particleTransferPos[i];
        }
    }

    labelListList allNTrans(Pstream::nProcs());

    pBufs.finishedSends(allNTrans);

    forAll(allNTrans, i)
    {
        forAll(allNTrans[i], j)
        {
            if (allNTrans[i][j])
            {
                break;
            }
        }
    }

    label neighbProci;
    label nRec;
    forAll(neighbourProcs_, i)
    {
        neighbProci = neighbourProcs_[i];
        nRec = allNTrans[neighbProci][Pstream::myProcNo()];

        if (nRec)
        {
            UIPstream particleStream(neighbProci, pBufs);

            labelList recvParticleTransferID(particleStream);
            List<vector> recvParticleTransferPos(particleStream);

            forAll(recvParticleTransferID,i)
            {
                pos = recvParticleTransferPos[i];
                cellID = particleCloud_.locateM().findSingleCell(pos,searchCellID);

                if (cellID >= 0)
                {
                    id_foam_[nlocal_foam_] = recvParticleTransferID[i];

                    for (int j=0;j<3;j++)
                    {
                        id_foamVec_[nlocal_foam_*3+j] = recvParticleTransferID[i]*3+j;
                        pos_foam_[nlocal_foam_*3+j] = pos[j];
                    }
                    cellID_foam_[nlocal_foam_] = cellID;

                    nlocal_foam_ += 1;
                }
                else
                {
                    id_foamLost_[nlocal_foam_lost_] = recvParticleTransferID[i];

                    for (int j=0; j<3; j++)
                        lost_pos_[nlocal_foam_lost_*3+j] = pos[j];

                    nlocal_foam_lost_ += 1;
                }

            }
        }
    }
    particleCloud_.clockM().stop("locate_Stage2");

    particleCloud_.clockM().start(9,"locate_Stage3");

    int nlocal_foam_lostAll(-1);
    if (firstRun_ || allowDiagComm_)
    {
        particleCloud_.clockM().start(10,"locate_Stage3_1");
        MPI_Allreduce(&nlocal_foam_lost_, &nlocal_foam_lostAll, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        particleCloud_.clockM().stop("locate_Stage3_1");
    }

    if (nlocal_foam_lostAll > 0)
    {
        Info << "all-to-all necessary: nlocal_foam_lostAll=" << nlocal_foam_lostAll << endl;
        if (lost_posAll)
        {
           delete[] lost_posAll;
           lost_posAll = NULL;
        }
        if (id_foamLostAll)
        {
           delete[] id_foamLostAll;
           id_foamLostAll = NULL;
        }
        int nlocal_foam_lostAll = LAMMPS_NS::MPI_Allgather_Vector(lost_pos_, nlocal_foam_lost_*3, lost_posAll, MPI_COMM_WORLD)/3;
        LAMMPS_NS::MPI_Allgather_Vector(id_foamLost_, nlocal_foam_lost_, id_foamLostAll, MPI_COMM_WORLD);
        Info << couplingStep_ << "st nlocal_foam_lostAll=" << nlocal_foam_lostAll << endl;

        for (int i = 0; i < nlocal_foam_lostAll; i++)
        {
            pos = vector(lost_posAll[i*3+0],lost_posAll[i*3+1],lost_posAll[i*3+2]);
            cellID = particleCloud_.locateM().findSingleCell(pos,searchCellID);

            if (cellID >= 0)
            {
                id_foam_[nlocal_foam_] = id_foamLostAll[i];

                for (int j=0;j<3;j++)
                {
                    id_foamVec_[nlocal_foam_*3+j] = id_foamLostAll[i]*3+j;
                    pos_foam_[nlocal_foam_*3+j] = pos[j];
                }
                cellID_foam_[nlocal_foam_] = cellID;

                id_foamLostAll[i]=-1;

                nlocal_foam_ += 1;
            }
        }
    }
    particleCloud_.clockM().stop("locate_Stage3");

    particleLost_ = false;
    if (firstRun_)
    {
        int* id_foam_nowhere_all;
        dataExchangeModel::allocateArray(id_foam_nowhere_all,1,nlocal_foam_lostAll);
        MPI_Allreduce(id_foamLostAll, id_foam_nowhere_all, nlocal_foam_lostAll, MPI_INT, MPI_MIN, MPI_COMM_WORLD);

        int i = 0;
        while (i < nlocal_foam_lostAll)
        {
            if (id_foam_nowhere_all[i] > 0)
            {
                particleLost_=true;

                for (int j=0;j<nlocal_lammps_;j++)
                {
                    if (id_lammpsSync[j]==id_foam_nowhere_all[i])
                    {
                        id_lammpsSync[j] = id_lammpsSync[nlocal_lammps_-1];

                        for (int k=0;k<3;k++)
                            id_lammpsVec_[j*3+k] = id_lammpsVec_[(nlocal_lammps_-1)*3+k];

                        nlocal_lammps_ -= 1;
                        break;
                    }
                }
            }
            i++;
        }
        dataExchangeModel::destroy(id_foam_nowhere_all);
        id_foam_nowhere_all = NULL;
        if (id_lammps_alloc_flag) destroy(id_lammps_);
        id_lammps_ = NULL;
        allocateArray(id_lammps_,-1.,nlocal_lammps_);
        for (int i = 0; i < nlocal_lammps_; i++)
            id_lammps_[i]=id_lammpsSync[i];
    }
#elif defined(version16ext)
    Info << "M2M does not work with 1.6.x" << endl;
#endif
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
