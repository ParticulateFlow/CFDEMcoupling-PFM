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
    propsDict_(dict.subDict(typeName + "Props"))
{
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
    lmp2foam_ = new Many2Many(MPI_COMM_WORLD);
    lmp2foam_vec_ = new Many2Many(MPI_COMM_WORLD);
    foam2lmp_vec_ = new Many2Many(MPI_COMM_WORLD);
    nlocal_lammps_ = -1;
    id_lammps_ = NULL;
    id_lammpsComm_ = NULL;
    id_lammps_vec_ = NULL;
    nlocal_foam_ = -1;
    id_foam_ = NULL;
    id_foam_vec_ = NULL;
    nlocal_foam_lost_ = -1;
    id_foam_lost_ = NULL;
    id_foam_lost_all = NULL;
    id_foam_nowhere_all = NULL;
    lost_pos_ = NULL;
    lost_pos_all = NULL;
    cellID_foam_ = NULL;
    pos_foam_ = NULL;
    //Foam::dataExchangeModel::allocateArray(idHashTable_,0,100000);//liggghts_get_maxtag(lmp));    // define idHashTable
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

twoWayM2M::~twoWayM2M()
{
    free(id_lammps_);
    delete[] id_lammpsComm_;
    delete[] id_lammps_vec_;
    delete[] id_foam_vec_;
    delete[] id_foam_;
    delete[] id_foam_lost_;
    delete[] lost_pos_;
    delete[] cellID_foam_;
    delete[] pos_foam_;
    //delete[] lmp2foam_;
    //delete lmp;
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
        //if (nlocal_lammps_>0){        
            double **tata_ = (double **) lammps_extract_atom(lmp,charName);
            lmp2foam_vec_->exchange(&(tata_[0][0]), &(field[0][0]));
        /*}else{
            // is this else necessary?
            double *tata_ = (double *) lammps_extract_atom(lmp,charName);
            lmp2foam_vec_->exchange(tata_, &(field[0][0]));
        }*/
        //for (int i = 0; i < nlocal_foam_; i++)
        //    Pout << name <<"=" << tata_[i][0]<<","<<tata_[i][1]<<","<<tata_[i][2] <<endl;
    }else if (name != "x"){
        tmp_ = (double *) lammps_extract_atom(lmp,charName);
        lmp2foam_->exchange(tmp_, &(field[0][0]));
        //for (int i = 0; i < nlocal_foam_; i++)
        //    Pout << name <<"[0][i]=" << field[0][i] <<endl;
    }
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
        double **tata_=NULL;
        LAMMPS_NS::Fix *fix = NULL;
        fix = lmp->modify->find_fix_property(charName,"property/atom","vector",0,0,"cfd coupling",false);
        if(fix)
            tata_ = (double **) static_cast<LAMMPS_NS::FixPropertyAtom*>(fix)->array_atom;
        else
            Warning << "coupling fix not found!"<<endl;

        foam2lmp_vec_->exchange(&(field[0][0]),&(tata_[0][0]));

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

void Foam::twoWayM2M::allocateArray
(
    double**& array,
    double initVal,
    int width,
    int length
) const
{
    //if(length==-1) then LIGGGHTS uses own length data
    allocate_external_double(array, width,length,initVal,lmp);
}

void Foam::twoWayM2M::allocateArray
(
    double**& array,
    double initVal,
    int width,
    const char* length
) const
{
    allocate_external_double(array, width,max(particleCloud_.numberOfParticles(),1),initVal,lmp);
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
    //if(length==-1) then LIGGGHTS uses own length data
    allocate_external_int(array, width,length,initVal,lmp);
}

void Foam::twoWayM2M::allocateArray
(
    int**& array,
    int initVal,
    int width,
    const char* length
) const
{
    allocate_external_int(array, width,max(particleCloud_.numberOfParticles(),1),initVal,lmp);
}
//============

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
        particleCloud_.clockM().stop("CoupleSyncIDs()");

        // give nr of particles to cloud
        setNumberOfParticles(nlocal_foam_);

        // re-allocate arrays of cloud      
        particleCloud_.reAllocArrays();

        // give existing position and cellID data to cloud
        setPositions(nlocal_foam_,pos_foam_);
        setCellIDs(nlocal_foam_,cellID_foam_);
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
    // get data from lammps
    nlocal_lammps_ = *((int *) lammps_extract_global(lmp,"nlocal"));
    id_lammps_ = (int *) lammps_extract_atom(lmp,"id");

    // make a copy of liggghts ID array
    delete [] id_lammpsComm_;
    Foam::dataExchangeModel::allocateArray(id_lammpsComm_,0,nlocal_lammps_);
    for(int i=0;i<nlocal_lammps_;i++) id_lammpsComm_[i]=id_lammps_[i];
    
    Foam::dataExchangeModel::allocateArray(pos_lammps_,-1.,3,nlocal_lammps_);  // is this the right allocate???
    pos_lammps_ = (double **) lammps_extract_atom(lmp,"x");
    particleCloud_.clockM().stop("recv_DEM_ids");

    particleCloud_.clockM().start(6,"locateParticle()");
    locateParticle();
    particleCloud_.clockM().stop("locateParticle()");

        // output
        /*Info << "LAMMPS " << endl;
        for (int i = 0; i < nlocal_lammps_; i++)
        {
            Pout << "id_lammps_[" << i << "]=" << id_lammps_[i] << "  -  "<<endl;
            //Pout << "pos_lammps_ radius [" << i << "]=" << pos_lammps_[i] << endl;
        }*/
        /*for (int i = 0; i < nlocal_lammps_*3; i++)
        {
            Pout << "id_lammps_vec_[" << i << "]=" << id_lammps_vec_[i] << "  -  "<<endl;
        }*/

        /*Info << "FOAM "<< endl;
        for (int i = 0; i < nlocal_foam_; i++)
        {
            Pout << "id_foam_[" << i << "]=" << id_foam_[i] << "  -  "<<endl;
            //Pout << "pos_foam_ radius [" << i << "]=" << pos_foam_[i] << endl;
        }*/
        /*for (int i = 0; i < nlocal_foam_*3; i++)
        {
            Pout << "id_foam_vec_[" << i << "]=" << id_foam_vec_[i] << "  -  "<<endl;
        }*/
        Pout << "nlocal_lammps_=" << nlocal_lammps_ << endl;
        Pout << "nlocal_foam_=" << nlocal_foam_ << endl;

    // communicate lmp->foam
    particleCloud_.clockM().start(11,"setup_Comm");
    lmp2foam_->setup(nlocal_lammps_,id_lammpsComm_,nlocal_foam_,id_foam_);
    lmp2foam_vec_->setup(nlocal_lammps_*3,id_lammps_vec_,nlocal_foam_*3,id_foam_vec_);
    foam2lmp_vec_->setup(nlocal_foam_*3,id_foam_vec_,nlocal_lammps_*3,id_lammps_vec_);
    particleCloud_.clockM().stop("setup_Comm");
}

void Foam::twoWayM2M::locateParticle() const
{
    //?=======
    // if I do not creat new communicators it fails at liggghts porcessor transfer!
    // this should probably not be done that often!!!  
    lmp2foam_ = new Many2Many(MPI_COMM_WORLD);
    lmp2foam_vec_ = new Many2Many(MPI_COMM_WORLD);
    foam2lmp_vec_ = new Many2Many(MPI_COMM_WORLD);
    //?=======

    int nop = particleCloud_.numberOfParticles();

    // realloc array of lost particles
    if (particleCloud_.numberOfParticlesChanged())
    {
        // these arrays will be to long, but we do not know their length a priori
        delete[] id_foam_lost_;
        Foam::dataExchangeModel::allocateArray(id_foam_lost_,0,nop);

        delete [] lost_pos_;
        Foam::dataExchangeModel::allocateArray(lost_pos_,0.,nop*3);

        delete [] id_foam_;
        Foam::dataExchangeModel::allocateArray(id_foam_,0,nop);

        delete [] id_foam_vec_;
        Foam::dataExchangeModel::allocateArray(id_foam_vec_,0,nop*3);

        delete [] id_lammps_vec_;
        Foam::dataExchangeModel::allocateArray(id_lammps_vec_,0,nop*3);

        delete [] cellID_foam_;
        Foam::dataExchangeModel::allocateArray(cellID_foam_,0,nop);
        delete [] pos_foam_;
        Foam::dataExchangeModel::allocateArray(pos_foam_,0,nop*3);
    }
    else
    {
        // reset array to zero
        for (int i=0; i<nop; i++)
        {
            id_foam_lost_[i] = 0;
            for (int j=0; j<3; j++) lost_pos_[i*3+j] = 0.;
        }
    }

    // if global nr of particles changed:
    // realloc idHashTable_

    // loop all lmp particles
    nlocal_foam_ = 0;
    nlocal_foam_lost_ = 0;
    vector pos;
    label cellID = 0;
    label searchCellID;

    particleCloud_.clockM().start(7,"locate_Stage1");
    for (int i = 0; i < nlocal_lammps_; i++)
    {
        pos = vector(pos_lammps_[i][0],pos_lammps_[i][1],pos_lammps_[i][2]);
        searchCellID = -1;//idHashTable_[id_lammps_[i]];
        cellID = particleCloud_.locateM().findSingleCell(pos,searchCellID);
    
        // IDs for vectors
        for (int j=0;j<3;j++)
            id_lammps_vec_[i*3+j] = id_lammps_[i]*3+j;

        // found particle on cfd proc
        if (cellID >= 0)
        {
            // IDs for scalars
            id_foam_[nlocal_foam_] = id_lammps_[i];

            // IDs for vectors
            for (int j=0;j<3;j++)
            {
                id_foam_vec_[nlocal_foam_*3+j] = id_lammps_[i]*3+j;                
                pos_foam_[nlocal_foam_*3+j] = pos[j];
            }
            cellID_foam_[nlocal_foam_] = cellID;

            // cell id hash table
            //idHashTable_[id_lammps_[i]]=cellID;

            nlocal_foam_ += 1;
            //Pout << "found particle at pos=" << pos << endl;
        }
        else
        {
            id_foam_lost_[nlocal_foam_lost_] = id_lammps_[i];
          
            for (int j=0; j<3; j++)
                lost_pos_[nlocal_foam_lost_*3+j] = pos[j];

            nlocal_foam_lost_ += 1;
            //Pout << "cellID="<< cellID << " lost particle id="<< id_lammps_[i] <<" at pos=" << pos << endl;
        }
    }
    particleCloud_.clockM().stop("locate_Stage1");

    // using allgather to allreduce lost particles
    particleCloud_.clockM().start(8,"locate_Stage2");
    int nlocal_foam_lost_all = LAMMPS_NS::MPI_Allgather_Vector(lost_pos_, nlocal_foam_lost_*3, lost_pos_all, MPI_COMM_WORLD)/3;
    LAMMPS_NS::MPI_Allgather_Vector(id_foam_lost_, nlocal_foam_lost_, id_foam_lost_all, MPI_COMM_WORLD);

    // locate lost particles
    for (int i = 0; i < nlocal_foam_lost_all; i++)
    {
        pos = vector(lost_pos_all[i*3+0],lost_pos_all[i*3+1],lost_pos_all[i*3+2]);
        searchCellID = -1;//idHashTable_[id_foam_lost_all[i]];
        particleCloud_.clockM().start(9,"findSingleCell");   
        cellID = particleCloud_.locateM().findSingleCell(pos,searchCellID);
        particleCloud_.clockM().stop("findSingleCell");
        
        // found particle on cfd proc
        if (cellID >= 0)
        {
            // IDs for scalars
            id_foam_[nlocal_foam_] = id_foam_lost_all[i];

            // IDs for vectors
            for (int j=0;j<3;j++)
            {
                id_foam_vec_[nlocal_foam_*3+j] = id_foam_lost_all[i]*3+j;
                pos_foam_[nlocal_foam_*3+j] = pos[j];
            }
            cellID_foam_[nlocal_foam_] = cellID;

            // mark that ID was finally found
            id_foam_lost_all[i]=-1;

            // cell id hash table
            //idHashTable_[id_foam_lost_all[i]]=cellID;

            nlocal_foam_ += 1;
        }
    }
    particleCloud_.clockM().stop("locate_Stage2");

  /*  // check if really all particles were found
    particleCloud_.clockM().start(10,"locate_Stage3");
    Foam::dataExchangeModel::allocateArray(id_foam_nowhere_all,1,nlocal_foam_lost_all);
    MPI_Allreduce(id_foam_lost_all, id_foam_nowhere_all, nlocal_foam_lost_all, MPI_INT, MPI_MIN, MPI_COMM_WORLD);

    int i=0;
    while (i < nlocal_foam_lost_all)
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
                        id_lammps_vec_[j*3+k] = id_lammps_vec_[(nlocal_lammps_-1)*3+k];
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

    delete[] id_foam_nowhere_all;
    delete[] id_foam_lost_all;
    delete[] lost_pos_all;
}

void Foam::twoWayM2M::exchange(double* demDat, double* cfdDat) const
{
    lmp2foam_->exchange(demDat,cfdDat);//(pos_lammps,pos_foam);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
