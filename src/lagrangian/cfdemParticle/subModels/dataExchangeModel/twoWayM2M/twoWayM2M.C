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
#include "twoWayM2M.H"
#include "addToRunTimeSelectionTable.H"
#include "clockModel.H"
#include "memory.h"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(twoWayM2M, 0);

addToRunTimeSelectionTable
(
    dataExchangeModel,
    twoWayM2M,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
twoWayM2M::twoWayM2M
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
    forAll(neighbourProcs_, i) neighbourProcIndices_[neighbourProcs_[i]] = i;

    Info<<"Starting up LIGGGHTS for first time execution"<<endl;

    MPI_Comm_rank(MPI_COMM_WORLD,&me);
    MPI_Comm_size(MPI_COMM_WORLD,&nprocs);

    if (me < nprocs) liggghts = 1;
    else liggghts = MPI_UNDEFINED;

    MPI_Comm_split(MPI_COMM_WORLD,liggghts,0,&comm_liggghts);

    // open LIGGGHTS input script
    FILE *fp=NULL;
    if (me == 0)
    {
      // read path from dictionary
      const fileName liggghtsPath(propsDict_.lookup("liggghtsPath"));
      char * liggghtsPathChar = (char*)liggghtsPath.c_str();

      Info<<"Executing input script '"<< liggghtsPath.c_str() <<"'"<<endl;

      fp = fopen(liggghtsPathChar,"r");

      if (fp == NULL) {
        printf("ERROR: Could not open LIGGGHTS input script\n");
        MPI_Abort(MPI_COMM_WORLD,1);
      }
    }

    if (liggghts == 1) lmp = new LAMMPS_NS::LAMMPS(0,NULL,comm_liggghts);

    int n;
    char line[1024];
    while (1) {
      if (me == 0) {
        if (fgets(line,1024,fp) == NULL) n = 0;
        else n = strlen(line) + 1;
        if (n == 0) fclose(fp);
      }
      MPI_Bcast(&n,1,MPI_INT,0,MPI_COMM_WORLD);
      if (n == 0) break;
      MPI_Bcast(line,n,MPI_CHAR,0,MPI_COMM_WORLD);
      if (liggghts == 1) lmp->input->one(line);
    }

    // get DEM time step size
    DEMts_ = lmp->update->dt;
    checkTSsize();

    // m2m stuff
    firstRun_=true;
	safeRun_=false;
    lmp2foam_ = NULL;
    lmp2foam_vec_ = NULL;
    foam2lmp_vec_ = NULL;
    nlocal_lammps_ = -1;
    id_lammps_ = NULL;
    id_lammpsVec_ = NULL;
    nlocal_foam_ = -1;
    id_foam_ = NULL;
    id_foamVec_ = NULL;
    tmp_ = NULL;
    tmpI_ = NULL;
    pos_lammps_=NULL;
    nlocal_foam_lost_ = -1;
    id_foamLost_ = NULL;
    id_foamLostAll = NULL;
    lost_pos_ = NULL;
    lost_posAll = NULL;
    cellID_foam_ = NULL;
    pos_foam_ = NULL;
	if (propsDict_.found("safeRun")) safeRun_=true;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

twoWayM2M::~twoWayM2M()
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
    destroy(pos_lammps_);
    delete lmp2foam_;
    delete lmp2foam_vec_;
    delete foam2lmp_vec_;
    delete lmp;
}


// * * * * * * * * * * * * * * * private Member Functions  * * * * * * * * * * * * * //
char* twoWayM2M::wordToChar(word& inWord) const
{
    string HH = string(inWord);
    return const_cast<char*>(HH.c_str());
}


// * * * * * * * * * * * * * * * public Member Functions  * * * * * * * * * * * * * //
void twoWayM2M::getData
(
    word name,
    word type,
    double ** const& field,
    label step
) const
{
    char* charName = wordToChar(name);
    if ( type == "vector-atom" && name != "x")
    {
        double **tmp_ = (double **) lammps_extract_atom(lmp,charName);
        //for (int i = 0; i < nlocal_lammps_; i++)
        //    for(int j=0;j<3; j++)
        //       Pout << couplingStep_ << "st tmp_[" << i << "][j]=" << tmp_[i][j] << "  - name="<< name <<endl;

        lmp2foam_vec_->exchange(tmp_ ? tmp_[0] : NULL, field[0]);

        //for (int i = 0; i < nlocal_foam_; i++)
        //    for(int j=0;j<3; j++)
        //       Pout << couplingStep_ << "st field[" << i << "][j]=" << field[i][j] << "  - name="<< name <<endl;

    }else if (name != "x")
    {
        double *tmp_ = (double *) lammps_extract_atom(lmp,charName);

        //for (int i = 0; i < nlocal_lammps_; i++)
        //    Pout << couplingStep_ << "st tmp_[" << i << "]=" << tmp_[i] << "  -  name=" << name <<endl;

        lmp2foam_->exchange(tmp_, field[0]);

        //for (int i = 0; i < nlocal_foam_; i++)
        //    Pout << couplingStep_ << "st field[0][" << i << "]=" << field[0][i] << "  -  name=" << name <<endl;
    }
    //Info << "getData done for :" << name << endl;
}

void twoWayM2M::getData
(
    word name,
    word type,
    int ** const& field,
    label step
) const
{
    char* charName = wordToChar(name);
    char* charType = wordToChar(type);
    data_liggghts_to_of(charName,charType, lmp, (void*&) field,"int");
}

void twoWayM2M::giveData
(
    word name,
    word type,
    double ** const& field,
    const char* datatype
) const
{
    char* charName = wordToChar(name);
    if ( type == "vector-atom")
    {
        double **tmp_=NULL;
        LAMMPS_NS::Fix *fix = NULL;
        fix = lmp->modify->find_fix_property(charName,"property/atom","vector",0,0,"cfd coupling",false);
        if(fix)
            tmp_ = (double **) static_cast<LAMMPS_NS::FixPropertyAtom*>(fix)->array_atom;
        else
            Warning << "coupling fix not found!"<<endl;

        foam2lmp_vec_->exchange(field[0],tmp_ ? tmp_[0] : NULL);

        //==================
        //for(int index = 0;index <  nlocal_lammps_; ++index){
        //    vector forceField(field[index][0],field[index][1],field[index][2]); 
        //    vector tataField(tata_[index][0],tata_[index][1],tata_[index][2]); 
        //    Pout << "particle=" << index << " ,forceField=" << forceField<< " ,tataField=" << tataField <<  endl;
        //}
        //==================
    }else{
        Warning << "not implemented!"<<endl;
    }
}

//============
// double **
void Foam::twoWayM2M::allocateArray
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

void Foam::twoWayM2M::allocateArray
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
void Foam::twoWayM2M::destroy(double** array) const
{
    lmp->memory->destroy(array);
}
//============
// int **
void Foam::twoWayM2M::allocateArray
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

void Foam::twoWayM2M::allocateArray
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
void Foam::twoWayM2M::destroy(int** array) const
{
    lmp->memory->destroy(array);
}
//============
// double *
void Foam::twoWayM2M::allocateArray(double*& array, double initVal, int length) const
{
    int len = max(length,1);
    lmp->memory->grow(array, len, "m2m:dbl*");
    for (int i = 0; i < len; i++)
        array[i] = initVal;
}
void Foam::twoWayM2M::destroy(double* array) const
{
    lmp->memory->destroy(array);
}
//==============
// int *
void Foam::twoWayM2M::allocateArray(int*& array, int initVal, int length) const
{
    int len = max(length,1);
    lmp->memory->grow(array, len, "m2m:int*");
    for (int i = 0; i < len; i++)
        array[i] = initVal;
}
void Foam::twoWayM2M::destroy(int* array) const
{
    lmp->memory->destroy(array);
}
//==============


bool Foam::twoWayM2M::couple() const
{
    bool coupleNow = false;
    if (doCoupleNow())
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

                if(particleCloud_.liggghtsCommand()[i]().runCommand(couplingStep()))
                {
                    const char* command = particleCloud_.liggghtsCommand()[i]().command();
                    Info << "Executing command: '"<< command <<"'"<< endl;
                    lmp->input->one(command);
                }
            }
            particleCloud_.clockM().stop("LIGGGHTS");
            Info<<"LIGGGHTS finished"<<endl;
        }

        // give nr of particles to cloud
        double newNpart = liggghts_get_maxtag(lmp);
        setNumberOfParticles(newNpart);

        // m2m stuff
        particleCloud_.clockM().start(4,"CoupleSyncIDs()");
        syncIDs();
        firstRun_=false;
        particleCloud_.clockM().stop("CoupleSyncIDs()");

        // give nr of particles to cloud
        setNumberOfParticles(nlocal_foam_);

        // re-allocate arrays of cloud      
        particleCloud_.reAllocArrays();

        // give existing position and cellID data to cloud
        setPositions(nlocal_foam_,pos_foam_);
        setCellIDs(nlocal_foam_,cellID_foam_);
        //free(cellID_foam_); // cannot free here?
        //free(pos_foam_); // cannot free here?
    }
    return coupleNow;
}

int Foam::twoWayM2M::getNumberOfParticles() const
{
    return liggghts_get_maxtag(lmp);
}

int Foam::twoWayM2M::getNumberOfClumps() const
{
    Warning << "Foam::twoWayM2M::getNumberOfClumps() - changes necessary here" << endl;
    //return liggghts_get_maxtag_ms(lmp);
    return 1;
}

void Foam::twoWayM2M::syncIDs() const
{
    particleCloud_.clockM().start(5,"recv_DEM_ids");

    // update communication 
    delete lmp2foam_;
    delete lmp2foam_vec_;
    delete foam2lmp_vec_;
    lmp2foam_ = new Many2Many(MPI_COMM_WORLD);
    lmp2foam_vec_ = new Many2Many(MPI_COMM_WORLD);
    foam2lmp_vec_ = new Many2Many(MPI_COMM_WORLD);

    // get data from lammps
    nlocal_lammps_ = *((int *) lammps_extract_global(lmp,"nlocal"));
    int*  id_lammpsSync=NULL;
    double** pos_lammpsSync=NULL;
    if(firstRun_)
    {
        // get access to id array
        id_lammps_ = (int *) lammps_extract_atom(lmp,"id");

        // genereate vector IDs
        allocateArray(id_lammpsVec_,0,nlocal_lammps_*3);
        for (int i = 0; i < nlocal_lammps_; i++)
            for (int j=0;j<3;j++)
                id_lammpsVec_[i*3+j] = id_lammps_[i]*3+j;
        
        // get access to "x"
        pos_lammps_ = (double **) lammps_extract_atom(lmp,"x");
    }
    else
    {
        // re-arrange data using map
       
        // get access to id array
        id_lammpsSync = (int *) lammps_extract_atom(lmp,"id");

        // genereate vector IDs
        allocateArray(id_lammpsVec_,0,nlocal_lammps_*3);
        for (int i = 0; i < nlocal_lammps_; i++)
            for (int j=0;j<3;j++)
                id_lammpsVec_[i*3+j] = id_lammpsSync[i]*3+j;

        // make setup of m2m
        lmp2foam_->setup(nlocal_lammps_,id_lammpsSync,nlocal_foam_,id_foam_);
        lmp2foam_vec_->setup(nlocal_lammps_*3,id_lammpsVec_,nlocal_foam_*3,id_foamVec_);
        foam2lmp_vec_->setup(nlocal_foam_*3,id_foamVec_,nlocal_lammps_*3,id_lammpsVec_);

        // map data according to last TS
        allocateArray(id_lammps_,-1.,nlocal_foam_);
        allocateArray(tmpI_,-1.,nlocal_foam_);
        lmp2foam_->exchange(id_lammpsSync, tmpI_);
        for(int i=0;i<nlocal_foam_;i++)
            id_lammps_[i]=tmpI_[i];

        // get access to "x"
        pos_lammpsSync = (double **) lammps_extract_atom(lmp,"x");
        allocateArray(tmp_,-1.,nlocal_foam_*3);
        lmp2foam_vec_->exchange(pos_lammpsSync ? pos_lammpsSync[0] : NULL, tmp_);

        allocateArray(pos_lammps_,0,3,nlocal_foam_);
        for(int i=0;i<nlocal_foam_;i++)
            for(int j=0;j<3;j++)
                pos_lammps_[i][j]=tmp_[i*3+j];

        //for (int i = 0; i < nlocal_foam_; i++)
        //    Pout << couplingStep_ << "st pos exchanged:" <<"=" << pos_lammps_[i][0]<<","<<pos_lammps_[i][1]<<","<<pos_lammps_[i][2] <<endl;
    }
    particleCloud_.clockM().stop("recv_DEM_ids");

    particleCloud_.clockM().start(6,"locateParticle()");
    locateParticle();
    particleCloud_.clockM().stop("locateParticle()");

//MPI_Barrier(MPI_COMM_WORLD);
//Pout << couplingStep_ << "st == syncIDs  " << endl;
//if(couplingStep_==30){
//FatalError<<"stop!!!"<< abort(FatalError);
//}

        // output
        /*Info << "LAMMPS " << endl;
        for (int i = 0; i < nlocal_lammps_; i++)
        {
            if(firstRun_)
            {
                Pout << couplingStep_ << "st id_lammps_[" << i << "]=" << id_lammps_[i] << "  -  "<<endl;
            }else{
                Pout << couplingStep_ << "st id_lammpsSync[" << i << "]=" << id_lammpsSync[i] << "  -  "<<endl;
            }
        }
        for (int i = 0; i < nlocal_lammps_*3; i++)
        {
            Pout << couplingStep_ << "st id_lammpsVec_[" << i << "]=" << id_lammpsVec_[i] << "  -  "<<endl;
        }
        Info << "FOAM "<< endl;
        for (int i = 0; i < nlocal_foam_; i++)
        {
            Pout << couplingStep_ << "st id_foam_[" << i << "]=" << id_foam_[i] << "  -  "<<endl;
        }
        for (int i = 0; i < nlocal_foam_*3; i++)
        {
            Pout << couplingStep_ << "st id_foamVec_[" << i << "]=" << id_foamVec_[i] << "  -  "<<endl;
        }
        Pout << couplingStep_ << "st nlocal_lammps_=" << nlocal_lammps_ << endl;
        Pout << couplingStep_ << "st nlocal_foam_=" << nlocal_foam_ << endl;*/

    // correct mapping
    particleCloud_.clockM().start(11,"setup_Comm");

    // update communication
    delete lmp2foam_;
    delete lmp2foam_vec_;
    delete foam2lmp_vec_;
    lmp2foam_ = new Many2Many(MPI_COMM_WORLD);
    lmp2foam_vec_ = new Many2Many(MPI_COMM_WORLD);
    foam2lmp_vec_ = new Many2Many(MPI_COMM_WORLD);
 
    if(firstRun_)
    {
        lmp2foam_->setup(nlocal_lammps_,id_lammps_,nlocal_foam_,id_foam_);
        lmp2foam_vec_->setup(nlocal_lammps_*3,id_lammpsVec_,nlocal_foam_*3,id_foamVec_);
        foam2lmp_vec_->setup(nlocal_foam_*3,id_foamVec_,nlocal_lammps_*3,id_lammpsVec_);
        id_lammps_=NULL;    // free pointer from LIG
        pos_lammps_ = NULL; // free pointer from LIG
    }else
    {
        lmp2foam_->setup(nlocal_lammps_,id_lammpsSync,nlocal_foam_,id_foam_);
        lmp2foam_vec_->setup(nlocal_lammps_*3,id_lammpsVec_,nlocal_foam_*3,id_foamVec_);
        foam2lmp_vec_->setup(nlocal_foam_*3,id_foamVec_,nlocal_lammps_*3,id_lammpsVec_);
    }
    particleCloud_.clockM().stop("setup_Comm");
}

void Foam::twoWayM2M::locateParticle() const
{
    int nop = particleCloud_.numberOfParticles();

    // realloc array of lost particles     // these arrays will be too long, but we do not know their length a priori???
    allocateArray(id_foamLost_,0,nop);
    allocateArray(lost_pos_,0.,nop*3);
    allocateArray(id_foam_,0,nop);
    allocateArray(id_foamVec_,0,nop*3);
    allocateArray(cellID_foam_,0,nop);
    if(firstRun_)
		allocateArray(pos_foam_,0,nop*3);

    // stage 1 - look on proc or send or prepare for all-to-all
    particleCloud_.clockM().start(7,"locate_Stage1");
    int iterate;
    if(firstRun_) iterate=nlocal_lammps_;
    else iterate=nlocal_foam_;

    nlocal_foam_ = 0;
    nlocal_foam_lost_ = 0;
    vector pos;
    label cellID = 0;
    label searchCellID=-1;
    List< DynamicList<int> > particleTransferID(neighbourProcs_.size());
    List< DynamicList<vector> > particleTransferPos(neighbourProcs_.size());

    for (int i = 0; i < iterate; i++)
    {
        pos = vector(pos_lammps_[i][0],pos_lammps_[i][1],pos_lammps_[i][2]);
        cellID = particleCloud_.locateM().findSingleCell(pos,searchCellID);
        point oldPos(pos_foam_[nlocal_foam_*3+0],pos_foam_[nlocal_foam_*3+1],pos_foam_[nlocal_foam_*3+2]);

        // found particle on cfd proc
        if (cellID >= 0)
        {
            // IDs for scalars
            id_foam_[nlocal_foam_] = id_lammps_[i];

            // IDs for vectors
            for (int j=0;j<3;j++)
            {
                id_foamVec_[nlocal_foam_*3+j] = id_lammps_[i]*3+j;                
                pos_foam_[nlocal_foam_*3+j] = pos[j];
            }
            cellID_foam_[nlocal_foam_] = cellID;

            nlocal_foam_ += 1;
            //Pout << couplingStep_ << "st stage1 found particle at pos=" << pos << " ,id_lammps_[i]=" << id_lammps_[i] << endl;
        }
        else
        {
            // find out where particle has migrated (must have passed a CFD proc border)
            bool commPart=false;
            point newPos=pos;
            label nearestFace = particleCloud_.locateM().intersection(oldPos,newPos);
            //Pout << couplingStep_ << "st nearestFace="<< nearestFace << " oldPos="<< oldPos <<" at pos=" << pos << endl;

            if (nearestFace >= particleCloud_.mesh().nInternalFaces())
            {
                label patchI = pbm_.whichPatch(nearestFace);

                if (procPatchIndices_[patchI] != -1)
                {
                    label n = neighbourProcIndices_
                    [
                        refCast<const processorPolyPatch>
                        (
                            pbm_[patchI]
                        ).neighbProcNo()
                    ];
                    particleTransferID[n].append(id_lammps_[i]);
                    particleTransferPos[n].append(pos);
                    commPart=true;
					//Pout << couplingStep_ << "st communicating particle " << id_lammps_[i] << ", to proc# " << n << endl;
                }
            }
            if (!commPart)
            {
                // prepare for all to all comm
                id_foamLost_[nlocal_foam_lost_] = id_lammps_[i];
          
                for (int j=0; j<3; j++)
                    lost_pos_[nlocal_foam_lost_*3+j] = pos[j];

                nlocal_foam_lost_ += 1;
                //Pout << couplingStep_ << "st cellID="<< cellID << " lost particle id="<< id_lammps_[i] <<" at pos=" << pos << endl;
            }
        }
    }
    particleCloud_.clockM().stop("locate_Stage1");

    // stage 2 - recv particle, locate or all-to-all
    particleCloud_.clockM().start(8,"locate_Stage2");
 
    // Allocate transfer buffers
    PstreamBuffers pBufs(Pstream::nonBlocking);

    // Stream into send buffers
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

    // Set up transfers when in non-blocking mode. Returns sizes (in bytes)
    // to be sent/received.
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

    // Retrieve from receive buffers
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
                //Pout << couplingStep_ << " stage 2 received id="<<recvParticleTransferID[i] << "received pos="<<recvParticleTransferPos[i] << endl;
                pos = recvParticleTransferPos[i];
                cellID = particleCloud_.locateM().findSingleCell(pos,searchCellID);
        
                // found particle on cfd proc
                if (cellID >= 0)
                {
                    // IDs for scalars
                    id_foam_[nlocal_foam_] = recvParticleTransferID[i];

                    // IDs for vectors
                    for (int j=0;j<3;j++)
                    {
                        id_foamVec_[nlocal_foam_*3+j] = recvParticleTransferID[i]*3+j;
                        pos_foam_[nlocal_foam_*3+j] = pos[j];
                    }
                    cellID_foam_[nlocal_foam_] = cellID;

                    // mark that ID was finally found
                    //id_foamLostAll[i]=-1;

                    nlocal_foam_ += 1;
                    //Pout << couplingStep_ << "st stage2 found particle at pos=" << pos << " ,id_foam_[i]=" << id_foam_[i] << endl;
                }
                else // might have been comm to wrong proc
                {
                    // prepare for all to all comm
                    id_foamLost_[nlocal_foam_lost_] = recvParticleTransferID[i];
          
                    for (int j=0; j<3; j++)
                        lost_pos_[nlocal_foam_lost_*3+j] = pos[j];

                    nlocal_foam_lost_ += 1;
                }

            }
        }
    }
    particleCloud_.clockM().stop("locate_Stage2");

    // stage 3 - all-to-all
    particleCloud_.clockM().start(9,"locate_Stage3");

    // check if all-to-all is necessary
    int nlocal_foam_lostAll(-1);
	if (firstRun_ || safeRun_)
    {
        particleCloud_.clockM().start(10,"locate_Stage3_1");
	    MPI_Allreduce(&nlocal_foam_lost_, &nlocal_foam_lostAll, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        particleCloud_.clockM().stop("locate_Stage3_1");
    }

    if (nlocal_foam_lostAll > 0)
    {
        Info << "all-to-all necessary: nlocal_foam_lostAll=" << nlocal_foam_lostAll << endl;
        if(lost_posAll) 
        { 
           delete[] lost_posAll;
           lost_posAll = NULL;
        }
        if(id_foamLostAll) 
        { 
           delete[] id_foamLostAll;
           id_foamLostAll = NULL;
        }
        int nlocal_foam_lostAll = LAMMPS_NS::MPI_Allgather_Vector(lost_pos_, nlocal_foam_lost_*3, lost_posAll, MPI_COMM_WORLD)/3; // new[] für lost_posAll!!!
        LAMMPS_NS::MPI_Allgather_Vector(id_foamLost_, nlocal_foam_lost_, id_foamLostAll, MPI_COMM_WORLD);
        //Info << couplingStep_ << "st nlocal_foam_lostAll=" << nlocal_foam_lostAll << endl;

        // locate lost particles
        for (int i = 0; i < nlocal_foam_lostAll; i++)
        {
            pos = vector(lost_posAll[i*3+0],lost_posAll[i*3+1],lost_posAll[i*3+2]);
            //Pout << "stage3 look for particle at pos=" << pos << endl;
            cellID = particleCloud_.locateM().findSingleCell(pos,searchCellID);
        
            // found particle on cfd proc
            if (cellID >= 0)
            {
                // IDs for scalars
                id_foam_[nlocal_foam_] = id_foamLostAll[i];

                // IDs for vectors
                for (int j=0;j<3;j++)
                {
                    id_foamVec_[nlocal_foam_*3+j] = id_foamLostAll[i]*3+j;
                    pos_foam_[nlocal_foam_*3+j] = pos[j];
                }
                cellID_foam_[nlocal_foam_] = cellID;

                // mark that ID was finally found
                //id_foamLostAll[i]=-1;

                nlocal_foam_ += 1;
                //Pout << "stage3 found particle at pos=" << pos << " ,id="<< id_foamLostAll[i] << endl;
            }
        }
    }
    particleCloud_.clockM().stop("locate_Stage3");

  /*  // check if really all particles were found
    particleCloud_.clockM().start(10,"locate_Stage3");
    Foam::dataExchangeModel::allocateArray(id_foam_nowhere_all,1,nlocal_foam_lostAll);
    MPI_Allreduce(id_foamLostAll, id_foam_nowhere_all, nlocal_foam_lostAll, MPI_INT, MPI_MIN, MPI_COMM_WORLD);

    int i=0;
    while (i < nlocal_foam_lostAll)
    {
        // these particles where found nowhere
        if (id_foam_nowhere_all[i] > 0)
        {
            for (int j=0;j<nlocal_lammps_;j++)
            {
                if (id_lammpsComm_[j]==id_foam_nowhere_all[i])
                {
                    // re-arrange IDs
                    id_lammpsComm_[j] = id_lammpsComm_[nlocal_lammps_-1];

                    // re-arrange IDs for vectors
                    for (int k=0;k<3;k++)
                    {
                        id_lammpsVec_[j*3+k] = id_lammpsVec_[nlocal_lammps_*3+k];
                    }
                
                    nlocal_lammps_ -= 1;
                    break;
                }
            }
        }
        i++;
    }
    particleCloud_.clockM().stop("locate_Stage3");*/    

/*int gaga;
MPI_Allreduce(&nlocal_foam_, &gaga, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
Info << "nlocal_foam_ALL=" << gaga << endl;

int gugu;
MPI_Allreduce(&nlocal_lammps_, &gugu, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
Info << "nlocal_lammps_ALL=" << gugu << endl;*/
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
