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

    two way DEM-CFD coupling via many 2 many scheme

Class
    twoWayMany2Many

SourceFiles
    twoWayMany2Many.C

\*---------------------------------------------------------------------------*/

#ifndef twoWayMany2Many_H
#define twoWayMany2Many_H

#include "dataExchangeModel.H"
#include "liggghtsCommandModel.H"
#include "OFstream.H"
#include <sys/stat.h>
#include "pair.h"
#include "force.h"
#include "forceModel.H"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#include <lammps.h>
#include <input.h>
#include <atom.h>
#include <library.h>
#include <error.h>
#include <library_cfd_coupling.h>
#include <update.h>

#include <many2many.h>
#include <mpi_liggghts.h>
#include "locateModel.H"

#include <fix.h>
#include <modify.h>
#include <fix_property_atom.h>
#include "meshSearch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class noDrag Declaration
\*---------------------------------------------------------------------------*/

class twoWayMany2Many
:
    public dataExchangeModel
{
private:

  // private data
    dictionary propsDict_;

    // new vars
    mutable bool firstRun_;
    Switch allowDiagComm_;
    mutable bool particleLost_;
    mutable double Npart_;
    mutable Many2Many * lmp2foam_;
    mutable Many2Many * lmp2foam_vec_;
    mutable Many2Many * foam2lmp_vec_;
    mutable Many2Many * foam2lmp_;
    mutable int nlocal_lammps_;
    mutable int *id_lammps_;
    mutable int *id_lammpsVec_;
    mutable int nlocal_foam_;
    mutable int *id_foam_;
    mutable int *id_foamVec_;
    mutable double *tmp_;
    mutable int *tmpI_;
    mutable double **pos_lammps_;
    mutable int nlocal_foam_lost_;
    mutable int *id_foamLost_;
    mutable int *id_foamLostAll;
    mutable double *lost_pos_;
    mutable double *lost_posAll;
    mutable int *cellID_foam_;
    mutable double *pos_foam_;

    const polyBoundaryMesh& pbm_;
    const globalMeshData& pData_;
    const labelList& procPatches_;
    const labelList& procPatchIndices_;
    const labelList& neighbourProcs_;
    labelList neighbourProcIndices_;

    // variables
    int me;

    int nprocs;

    int liggghts;

    MPI_Comm comm_liggghts;

    LAMMPS_NS::LAMMPS *lmp;

  // private member functions

public:

    //- Runtime type information
    TypeName("twoWayMany2Many");


    // Constructors

        //- Construct from components
        twoWayMany2Many
        (
            const dictionary& dict,
            cfdemCloud& sm
        );

    // Destructor

        ~twoWayMany2Many();


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
        void inline destroy(double**,int) const;
        //============
        // int **
        void allocateArray(int**&, int, int, int) const;
        void allocateArray(int**&, int, int,const char* = "nparticles") const;
        void inline destroy(int**,int) const;
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

        bool couple(int);

        int getNumberOfParticles() const;
        int getNumberOfClumps() const;

        void syncIDs() const;
        void locateParticle(int*, bool) const;
        scalar getCG() const { return lmp->force->cg(); }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
