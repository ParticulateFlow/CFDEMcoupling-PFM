/*---------------------------------------------------------------------------*\
    CFDEMcoupling academic - Open Source CFD-DEM coupling
    
    Contributing authors:
    Thomas Lichtenegger
    Copyright (C) 2015- Johannes Kepler University, Linz
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
    Foam::cfdemCloudRec

SourceFiles
    cfdemCloudRec.C

\*---------------------------------------------------------------------------*/

#ifndef cfdemCloudRec_H
#define cfdemCloudRec_H

#include "cfdemCloud.H"
#include "recModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
class recModel;
/*---------------------------------------------------------------------------*\
                           Class cfdemCloudRec Declaration
\*---------------------------------------------------------------------------*/

class cfdemCloudRec : public cfdemCloud
{
protected:

    autoPtr<recModel> recModel_;

    void getDEMdata();

    void giveDEMdata();
    
    void setForces();

public:

    friend class dataExchangeModel;
    friend class voidFractionModel;
    friend class forceModel;
    friend class forceSubModel;

// Constructors

    //- Construct from mesh and a list of particles
    cfdemCloudRec
    (
            const fvMesh& mesh
    );

    //- Destructor
    virtual ~cfdemCloudRec();

// public Member Functions

    virtual bool evolve(volScalarField&,volVectorField&,volVectorField&);
    
    inline const recModel& recM() const
    {
        return recModel_;
    }

    void initRecFields(volScalarField&,volScalarField&,volVectorField&,volVectorField&);
    void updateRecFields(volScalarField&,volScalarField&,volVectorField&,volVectorField&);
    
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //



#endif

// ************************************************************************* //
