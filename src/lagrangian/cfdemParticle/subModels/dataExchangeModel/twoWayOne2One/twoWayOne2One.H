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
    enhanced two way DEM-CFD coupling via MPI.

    Compared to twoWayMPI, no Allreduces are used for communication.
    Instead, a geomatric map between FOAM and LIG domains is created and
    subsequently used for communication. 

Class
    twoWayOne2One

SourceFiles
    twoWayOne2One.C

Contributing authors
    Paul Kieckhefen (TUHH)  2018

\*---------------------------------------------------------------------------*/

#ifndef twoWayOne2One_H
#define twoWayOne2One_H

#include "dataExchangeModel.H"
#include "liggghtsCommandModel.H"
#include "OFstream.H"
#include <sys/stat.h>
#include "pair.h"
#include "force.h"
#include "forceModel.H"
#include "one2one.H"

#include "meshSearch.H"

//=================================//
//LAMMPS/LIGGGHTS
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#include <lammps.h>         // these are LAMMPS include files
#include <input.h>
#include <atom.h>
#include <memory.h>
#include <library.h>
#include <library_cfd_coupling_one2one.h>
#include <update.h>
#include <comm.h>
#include <fix.h>
#include <fix_property_atom.h>
//=================================//

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class noDrag Declaration
\*---------------------------------------------------------------------------*/

class twoWayOne2One
:
    public dataExchangeModel
{
private:

  // private data
    dictionary        propsDict_;

    mutable MPI_Comm  comm_liggghts_;


    // LIG ranks from which to retrieve particle data
    mutable labelList thisLigPartner_;
    mutable labelList thisFoamPartner_;

    mutable One2One*  lig2foam_;
    mutable One2One*  foam2lig_;

    mutable bool*     lig2foam_mask_;

    mutable int*      lig2foam_ids_;
    mutable int*      foam2lig_ids_;

    mutable double*   lig2foam_vec_tmp_;
    mutable double*   lig2foam_scl_tmp_;

    mutable double*   foam2lig_vec_tmp_;
    mutable double*   foam2lig_scl_tmp_;

    mutable Switch    staticProcMap_;
    mutable Switch    cellIdComm_;
    mutable double*   prev_cell_ids_;
    mutable double**  dbl_cell_ids_;
    mutable LAMMPS_NS::FixPropertyAtom* my_prev_cell_ids_fix_;

    const scalar      boundBoxMargin_;

    treeBoundBox      thisFoamBox_;

    Switch            verbose_;


  // private member functions

    //- creates a geometric mapping between FOAM and LIG domains
    void createProcMap() const;

    //- create a One2One communicator which transfers from LIG to FOAM
    void setupLig2FoamCommunication() const;

    //- locates particles received from Lig
    void locateParticles() const;

    //- create a One2One communicator which transfers from FOAM to LIG
    void setupFoam2LigCommunication() const;

protected:
    LAMMPS_NS::LAMMPS *lmp;

public:

    //- Runtime type information
    TypeName("twoWayOne2One");


    // Constructors

        //- Construct from components
        twoWayOne2One
        (
            const dictionary& dict,
            cfdemCloud& sm
        );

    // Destructor

        ~twoWayOne2One();


    // Member Functions
        void getData
        (
            word name,
            word type,
            double ** const& field,
            label step
        ) const;

        void getData
        (
            word name,
            word type,
            int ** const& field,
            label step
        ) const;

        void giveData
        (
            word name,
            word type,
            double ** const& field,
            const char* datatype
        ) const;

        //============
        // double **
        void allocateArray(double**&, double, int, int) const;
        void allocateArray(double**&, double, int,const char* = "nparticles") const;
        void inline destroy(double**,int=0) const;
        //============
        // int **
        void allocateArray(int**&, int, int, int) const;
        void allocateArray(int**&, int, int,const char* = "nparticles") const;
        void inline destroy(int**,int=0) const;
        //==============

        //==============
        // double *
        void allocateArray(double*&, double, int) const;
        void inline destroy(double*) const;
        //==============
        // int *
        void allocateArray(int*&, int, int) const;
        void inline destroy(int*) const;
        //==============

        bool couple(int) const;

        //- extractCollected takes the collected data from Lig
        //  present in this Foam domain and applies the mask.
        //  the width parameter can be used for reshaping. 
        template <typename T>
        void extractCollected(T**&, T**&, int width=1) const;
        template <typename T>
        void extractCollected(T*&, T*&, int width=1) const;
        template <typename T>
        void extractCollected(T*&, T**&, int width=1) const;

        int getNumberOfParticles() const;

        word myType() const { return typeName; }

        void setCG() { particleCloud_.setCG(lmp->force->cg()); }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
