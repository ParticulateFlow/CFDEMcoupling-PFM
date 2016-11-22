/*---------------------------------------------------------------------------*\
    CFDEMcoupling academic - Open Source CFD-DEM coupling
    
    Contributing authors:
    Thomas Lichtenegger
    Copyright (C) 2015- Johannes Kepler University, Linz
-------------------------------------------------------------------------------
License
    This file is part of CFDEMcoupling academic.

    CFDEMcoupling academic is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    CFDEMcoupling academic is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with CFDEMcoupling academic.  If not, see <http://www.gnu.org/licenses/>.
\*---------------------------------------------------------------------------*/

#include "recBase.H"
#include "recModel.H"
#include "recNorm.H"
#include "recPath.H"

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
recBase::recBase
(
    const fvMesh& mesh
)
:
    mesh_(mesh),
    recProperties_
    (
        IOobject
        (
            "recProperties",
            mesh_.time().constant(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    recModel_
    (
        recModel::New
        (
            recProperties_,
            *this
        )
    ),
    recNorm_
    (
        recNorm::New
        (
            recProperties_,
            *this
        )
    ),
    recPath_
    (
        recPath::New
        (
            recProperties_,
            *this
        )
    ),
    couplingSubStep_(recProperties_.lookupOrDefault<label>("couplingSubStep",0))
{
  recModel_ ->  readFieldSeries();
  recNorm_  ->  computeRecMatrix();
  recPath_  ->  getRecPath();
 
  recModel_ ->  init();
  
  recModel_ ->  writeRecMatrix();    
  recModel_ ->  writeRecPath();
}

// * * * * * * * * * * * * * * * * Destructors  * * * * * * * * * * * * * * //
recBase::~recBase()
{}

// * * * * * * * * * * * * * * * public Member Functions  * * * * * * * * * * * * * //

const fvMesh& recBase::mesh() const
{
    return mesh_;
}

recModel& recBase::recM()
{
    return recModel_();
}

void recBase::updateRecFields()
{
    recModel_->updateRecFields(); 
}

}

label recBase::couplingSubStep() const
{
    return couplingSubStep_;
}