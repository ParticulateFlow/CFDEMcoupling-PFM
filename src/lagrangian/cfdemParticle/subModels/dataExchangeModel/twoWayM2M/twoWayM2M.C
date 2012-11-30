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
    propsDict_(dict.subDict(typeName + "Props"))/*,
    lmp2foam_(*new Many2Many(MPI_COMM_WORLD)), // init of many2many &
    lmp2foam_vec_(*new Many2Many(MPI_COMM_WORLD)),
    foam2lmp_vec_(*new Many2Many(MPI_COMM_WORLD))*/
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
    firstRun_=true;
    lmp2foam_ = NULL;
    lmp2foam_vec_ = NULL;
    foam2lmp_vec_ = NULL;
    nlocal_lammps_ = -1;
    id_lammps_ = NULL;
    //id_lammpsComm_ = NULL;
    id_lammps_vec_ = NULL;
    nlocal_foam_ = -1;
    id_foam_ = NULL;
    id_foam_vec_ = NULL;
    pos_lammps_=NULL;
    nlocal_foam_lost_ = -1;
    id_foam_lost_ = NULL;
    id_foam_lost_all = NULL;
    id_foam_nowhere_all = NULL;
    lost_pos_ = NULL;
    lost_pos_all = NULL;
    cellID_foam_ = NULL;
    pos_foam_ = NULL;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

twoWayM2M::~twoWayM2M()
{
    //delete[] id_lammps_; // does not work?
    //delete[] id_lammpsComm_;
    delete[] id_lammps_vec_;
    delete[] id_foam_vec_;
    delete[] id_foam_;
    delete[] id_foam_lost_;
    delete[] lost_pos_;
    delete[] cellID_foam_;
    delete[] pos_foam_;
    //delete& lmp2foam_;  // suitable for m2m&
    //delete& lmp2foam_vec_;
    //delete& foam2lmp_vec_;
    delete[] lmp2foam_;
    delete[] lmp2foam_vec_;
    delete[] foam2lmp_vec_;
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
        double **tata_ = (double **) lammps_extract_atom(lmp,charName);
        lmp2foam_vec_->exchange(&(tata_[0][0]), &(field[0][0]));
        //for (int i = 0; i < nlocal_foam_; i++)
        //    Pout << "hihi getData: " << name <<"=" << field[i][0]<<","<<field[i][1]<<","<<field[i][2] <<endl;
        //    Pout << name <<"=" << tata_[i][0]<<","<<tata_[i][1]<<","<<tata_[i][2] <<endl;
    }else if (name != "x"){
        //if(nlocal_lammps_>0){
            tmp_ = (double *) lammps_extract_atom(lmp,charName);
        /*}else{
            // might use the fct from dataExchangeModel mother class
            tmp_=new double[1];
            tmp_[0]=0;
        }*/

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
    allocate_external_double(array, width,max(length,1),initVal,lmp);
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
    allocate_external_int(array, width,max(length,1),initVal,lmp);
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
        firstRun_=false;
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
    // if I do not creat new communicators it fails at liggghts porcessor transfer!
    // this should probably not be done that often!!!  
    delete lmp2foam_;
    delete lmp2foam_vec_;
    delete foam2lmp_vec_;
    lmp2foam_ = new Many2Many(MPI_COMM_WORLD);
    lmp2foam_vec_ = new Many2Many(MPI_COMM_WORLD);
    foam2lmp_vec_ = new Many2Many(MPI_COMM_WORLD);
    //?=======

//MPI_Barrier(MPI_COMM_WORLD);
//Pout << couplingStep_ << "st == syncIDs  " << endl;
//if(couplingStep_==30){
//FatalError<<"stop!!!"<< abort(FatalError);
//}

    particleCloud_.clockM().start(5,"recv_DEM_ids");
    // get data from lammps
    nlocal_lammps_ = *((int *) lammps_extract_global(lmp,"nlocal"));
    int*  id_lammps_sync;
    double** pos_lammps_sync;
    if(firstRun_)  // do not forget iterator !
    {
        // IDs for vectors
        if(nlocal_lammps_>0){
            id_lammps_ = (int *) lammps_extract_atom(lmp,"id");
        }else{
            // might use the fct from dataExchangeModel mother class
            id_lammps_=new int[1];
            id_lammps_[0]=0;
        }

        Pout << couplingStep_ << "st id_lammps_[0]=" << id_lammps_[0]<< endl;


        //delete [] id_lammps_vec_;
        Foam::dataExchangeModel::allocateArray(id_lammps_vec_,0,nlocal_lammps_*3);
        for (int i = 0; i < nlocal_lammps_; i++)
            for (int j=0;j<3;j++)
                id_lammps_vec_[i*3+j] = id_lammps_[i]*3+j;
        
        Foam::dataExchangeModel::allocateArray(pos_lammps_,-1.,3,nlocal_lammps_);  // do I need this???
        if(nlocal_lammps_>0){
            pos_lammps_ = (double **) lammps_extract_atom(lmp,"x"); 
        }else{
            // might use the fct from dataExchangeModel mother class
            pos_lammps_ = new double*[1];
            pos_lammps_[0] = new double [1];
            pos_lammps_[0][0] = 0;
        }              
    }
    else
    {
       // re-arrange data using map
        //Foam::dataExchangeModel::allocateArray(id_lammps_sync,-1.,nlocal_lammps_); // probably not necessary
        if(nlocal_lammps_>0){
            id_lammps_sync = (int *) lammps_extract_atom(lmp,"id");
        }else{
            // might use the fct from dataExchangeModel mother class
            id_lammps_sync=new int[1];
            id_lammps_sync[0]=10;
        }
        //extract_save(id_lammps_sync,"id"); // in future it should look like this!!!
        Foam::dataExchangeModel::allocateArray(id_lammps_vec_,0,nlocal_lammps_*3);
        for (int i = 0; i < nlocal_lammps_; i++)
            for (int j=0;j<3;j++)
                id_lammps_vec_[i*3+j] = id_lammps_sync[i]*3+j;

        // make setup of m2m
        lmp2foam_->setup(nlocal_lammps_,id_lammps_sync,nlocal_foam_,id_foam_);
        lmp2foam_vec_->setup(nlocal_lammps_*3,id_lammps_vec_,nlocal_foam_*3,id_foam_vec_);
        foam2lmp_vec_->setup(nlocal_foam_*3,id_foam_vec_,nlocal_lammps_*3,id_lammps_vec_);

        // map data according to last TS
        Foam::dataExchangeModel::allocateArray(id_lammps_,0,nlocal_foam_);
        lmp2foam_->exchange(id_lammps_sync, id_lammps_);

        if(nlocal_lammps_>0){
            pos_lammps_sync = (double **) lammps_extract_atom(lmp,"x");
        }else{
            // might use the fct from dataExchangeModel mother class
            pos_lammps_sync = new double*[1];
            pos_lammps_sync[0] = new double [3];
            pos_lammps_sync[0][0] = 0;
            pos_lammps_sync[0][1] = 0;
            pos_lammps_sync[0][2] = 0;
        }
        // find better solution here!!!
        //Foam::dataExchangeModel::allocateArray(pos_lammps_,-1.,3,nlocal_foam_);
        //lmp2foam_vec_->exchange(&(pos_lammps_sync[0][0]), &(pos_lammps_[0][0]));
        double** gugu;
        Foam::dataExchangeModel::allocateArray(gugu,-1.,3*nlocal_foam_,1);
        Foam::dataExchangeModel::allocateArray(pos_lammps_,-1.,3,nlocal_foam_);
        lmp2foam_vec_->exchange(&(pos_lammps_sync[0][0]), &(gugu[0][0]));

        // conversion of array (should not be necessary if above problem is solved)
        for (int i = 0; i < nlocal_foam_; i++)
            for (int j = 0; j < 3; j++)
                pos_lammps_[i][j]=gugu[0][i*3+j];
    }
    particleCloud_.clockM().stop("recv_DEM_ids");

    particleCloud_.clockM().start(6,"locateParticle()");
    locateParticle();
    particleCloud_.clockM().stop("locateParticle()");

        //for (int i = 0; i < nlocal_lammps_; i++)
        //    Pout << "getData3: " << "v" <<"=" << pos_lammps_[i][0]<<","<<pos_lammps_[i][1]<<","<<pos_lammps_[i][2] <<endl;

        // output
        /*Info << "LAMMPS " << endl;
        for (int i = 0; i < nlocal_lammps_; i++)
        {
            if(firstRun_)
            {
                Pout << couplingStep_ << "st id_lammps_[" << i << "]=" << id_lammps_[i] << "  -  "<<endl;
            }else{
                Pout << couplingStep_ << "st id_lammps_sync[" << i << "]=" << id_lammps_sync[i] << "  -  "<<endl;
            }
        }*/
        /*for (int i = 0; i < nlocal_lammps_*3; i++)
        {
            Pout << couplingStep_ << "st id_lammps_vec_[" << i << "]=" << id_lammps_vec_[i] << "  -  "<<endl;
        }*/
        /*Info << "FOAM "<< endl;
        for (int i = 0; i < nlocal_foam_; i++)
        {
            Pout << couplingStep_ << "st id_foam_[" << i << "]=" << id_foam_[i] << "  -  "<<endl;
        }*/
        /*for (int i = 0; i < nlocal_foam_*3; i++)
        {
            Pout << couplingStep_ << "st id_foam_vec_[" << i << "]=" << id_foam_vec_[i] << "  -  "<<endl;
        }*/
        Pout << couplingStep_ << "st nlocal_lammps_=" << nlocal_lammps_ << endl;
        Pout << couplingStep_ << "st nlocal_foam_=" << nlocal_foam_ << endl;

    delete lmp2foam_;
    delete lmp2foam_vec_;
    delete foam2lmp_vec_;
    lmp2foam_ = new Many2Many(MPI_COMM_WORLD);
    lmp2foam_vec_ = new Many2Many(MPI_COMM_WORLD);
    foam2lmp_vec_ = new Many2Many(MPI_COMM_WORLD);

    // correct mapping
    particleCloud_.clockM().start(11,"setup_Comm");
    if(firstRun_)
    {
        //lmp2foam_->setup(nlocal_lammps_,id_lammpsComm_,nlocal_foam_,id_foam_);
        lmp2foam_->setup(nlocal_lammps_,id_lammps_,nlocal_foam_,id_foam_);
        lmp2foam_vec_->setup(nlocal_lammps_*3,id_lammps_vec_,nlocal_foam_*3,id_foam_vec_);
        foam2lmp_vec_->setup(nlocal_foam_*3,id_foam_vec_,nlocal_lammps_*3,id_lammps_vec_);
    }else
    {
        //lmp2foam_->setup(nlocal_lammps_,id_lammpsComm_,nlocal_foam_,id_foam_);
        lmp2foam_->setup(nlocal_lammps_,id_lammps_sync,nlocal_foam_,id_foam_);
        lmp2foam_vec_->setup(nlocal_lammps_*3,id_lammps_vec_,nlocal_foam_*3,id_foam_vec_);
        foam2lmp_vec_->setup(nlocal_foam_*3,id_foam_vec_,nlocal_lammps_*3,id_lammps_vec_);
    }
    particleCloud_.clockM().stop("setup_Comm");
}

void Foam::twoWayM2M::locateParticle() const
{
    int nop = particleCloud_.numberOfParticles();

    // realloc array of lost particles
    if (particleCloud_.numberOfParticlesChanged())
    {
        // these arrays will be to long, but we do not know their length a priori
        //delete[] id_foam_lost_;
        Foam::dataExchangeModel::allocateArray(id_foam_lost_,0,nop);

        //delete [] lost_pos_;
        Foam::dataExchangeModel::allocateArray(lost_pos_,0.,nop*3);

        //delete [] id_foam_;
        Foam::dataExchangeModel::allocateArray(id_foam_,0,nop);

        //delete [] id_foam_vec_;
        Foam::dataExchangeModel::allocateArray(id_foam_vec_,0,nop*3);

        //delete [] cellID_foam_;
        Foam::dataExchangeModel::allocateArray(cellID_foam_,0,nop);
        //delete [] pos_foam_;
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

    // loop all lmp particles
    int iterate;
    if(firstRun_) iterate=nlocal_lammps_;
    else iterate=nlocal_foam_;

    nlocal_foam_ = 0;
    nlocal_foam_lost_ = 0;
    vector pos;
    label cellID = 0;
    label searchCellID;

    particleCloud_.clockM().start(7,"locate_Stage1");

/*//===============
    // move this to top level
            const polyBoundaryMesh& pbm = particleCloud_.mesh().boundaryMesh(); 
            const globalMeshData& pData = particleCloud_.mesh().globalData(); // polyMesh???

            // Which patches are processor patches
            const labelList& procPatches = pData.processorPatches();

            // Indexing of patches into the procPatches list
            const labelList& procPatchIndices = pData.processorPatchIndices();

            // Which processors this processor is connected to
            const labelList& neighbourProcs = pData[Pstream::myProcNo()];

            // Indexing from the processor number into the neighbourProcs list
            labelList neighbourProcIndices(Pstream::nProcs(), -1);

            forAll(neighbourProcs, i)
            {
                neighbourProcIndices[neighbourProcs[i]] = i;
            }

            List< DynamicList<int> > particleTransferLists(neighbourProcs.size());
//===============*/

    for (int i = 0; i < iterate; i++)
    {
        pos = vector(pos_lammps_[i][0],pos_lammps_[i][1],pos_lammps_[i][2]);
        searchCellID = -1;
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
                id_foam_vec_[nlocal_foam_*3+j] = id_lammps_[i]*3+j;                
                pos_foam_[nlocal_foam_*3+j] = pos[j];
            }
            cellID_foam_[nlocal_foam_] = cellID;

            nlocal_foam_ += 1;
            //Pout << couplingStep_ << "st stage1 found particle at pos=" << pos << " ,id_lammps_[i]=" << id_lammps_[i] << endl;
        }
        else
        {
            //-----------------
            id_foam_lost_[nlocal_foam_lost_] = id_lammps_[i];
          
            for (int j=0; j<3; j++)
                lost_pos_[nlocal_foam_lost_*3+j] = pos[j];

            nlocal_foam_lost_ += 1;
            //Pout << couplingStep_ << "st cellID="<< cellID << " lost particle id="<< id_lammps_[i] <<" at pos=" << pos << endl;
            //-----------------
            /*// find out where particle has migrated (must have passed a CFD proc border)
            point newPos=pos;
            label nearestFace = particleCloud_.locateM().intersection(oldPos,newPos);
            Pout << "nearest face=" << nearestFace << endl;


            // If we hit a boundary face
            if (nearestFace >= particleCloud_.mesh().nInternalFaces())
            {
                Pout << " face=" << nearestFace << " , is a boundary face!" << endl;
                label patchI = pbm.whichPatch(nearestFace);

                // ... and the face is on a processor patch
                // prepare it for transfer
                if (procPatchIndices[patchI] != -1)
                {
                    Pout << " face=" << nearestFace << " , is a proc face!" << endl;
                    label n = neighbourProcIndices
                    [
                        refCast<const processorPolyPatch>
                        (
                            pbm[patchI]
                        ).neighbProcNo()
                    ];
                    Pout << " communicate to n=" << n << endl;
                    particleTransferLists[n].append(id_lammps_[i]);
                }
            }*/
            //-----------------

        }
    }

    /*forAll(particleTransferLists, i)
    {
        forAll(particleTransferLists[i],j)
        {
            Pout << "communicate particle   i="<< particleTransferLists[i][j]<<" to proc "<< i << endl;
        }
    }*/

    particleCloud_.clockM().stop("locate_Stage1");

    // using allgather to allreduce lost particles
    particleCloud_.clockM().start(8,"locate_Stage2");

    int nlocal_foam_lost_all = LAMMPS_NS::MPI_Allgather_Vector(lost_pos_, nlocal_foam_lost_*3, lost_pos_all, MPI_COMM_WORLD)/3;
    LAMMPS_NS::MPI_Allgather_Vector(id_foam_lost_, nlocal_foam_lost_, id_foam_lost_all, MPI_COMM_WORLD);
    Info << couplingStep_ << "st nlocal_foam_lost_all=" << nlocal_foam_lost_all << endl;

    // locate lost particles
    for (int i = 0; i < nlocal_foam_lost_all; i++)
    {
        pos = vector(lost_pos_all[i*3+0],lost_pos_all[i*3+1],lost_pos_all[i*3+2]);
        //Pout << "stage2 looking for particle pos="<< pos << endl;
        searchCellID = -1;
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

            nlocal_foam_ += 1;
            //Pout << "stage2 found particle at pos=" << pos << endl;
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
                        id_lammps_vec_[j*3+k] = id_lammps_vec_[nlocal_lammps_*3+k];
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

    //delete[] id_foam_nowhere_all;
    //delete[] id_foam_lost_all;
    //delete[] lost_pos_all;
}

/*void Foam::twoWayM2M::exchange(double* demDat, double* cfdDat) const
{
    lmp2foam_->exchange(demDat,cfdDat);//(pos_lammps,pos_foam);
}*/

/*template <typename T>
T*& Foam::twoWayM2M::extract_save(T *& a,char name)
{
    a = (T *) lammps_extract_atom(lmp,name);
    if(!a)
        a=new T[1];
    return a;//static_const<T>(a); 
}*/

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
