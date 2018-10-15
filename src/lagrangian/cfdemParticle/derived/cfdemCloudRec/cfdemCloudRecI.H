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

Description
    This code is designed to realize coupled rCFD-DEM simulations using LIGGGHTS
    and OpenFOAM(R). Note: this code is not part of OpenFOAM(R) (see DISCLAIMER).
\*---------------------------------------------------------------------------*/

//#include "cfdemCloudRec.H"
#include "forceModel.H"
namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
template <class baseCloud>
cfdemCloudRec<baseCloud>::cfdemCloudRec
(
    const fvMesh& mesh
)
:
    baseCloud(mesh),
    coupleRecFluc_(true),
    coupleRecForce_(true)
{
    if (baseCloud::couplingProperties().found("coupleRecFluc"))
    {
        coupleRecFluc_ = readBool(baseCloud::couplingProperties().lookup("coupleRecFluc"));
    }

    if (baseCloud::couplingProperties().found("coupleRecForce"))
    {
        coupleRecForce_ = readBool(baseCloud::couplingProperties().lookup("coupleRecForce"));
    }
}

// * * * * * * * * * * * * * * * * Destructors  * * * * * * * * * * * * * * //
template <class baseCloud>
cfdemCloudRec<baseCloud>::~cfdemCloudRec()
{}

// * * * * * * * * * * * * * * * protected Member Functions  * * * * * * * * * * * * * //
template <class baseCloud>
void cfdemCloudRec<baseCloud>::giveDEMdata()
{
    baseCloud::dataExchangeM().giveData("vrec","vector-atom",baseCloud::particleConvVel_);
    if(coupleRecFluc_)
    {
        baseCloud::dataExchangeM().giveData("vfluc","vector-atom",baseCloud::particleFlucVel_);
    }
    if(coupleRecForce_)
    {
        baseCloud::dataExchangeM().giveData("dragforce","vector-atom",baseCloud::DEMForces_);
    }
    if(baseCloud::verbose_) Info << "giveDEMdata done." << endl;
}

template <class baseCloud>
void cfdemCloudRec<baseCloud>::setForces()
{
    baseCloud::resetArray(baseCloud::particleConvVel_,baseCloud::numberOfParticles(),3);
    baseCloud::resetArray(baseCloud::particleFlucVel_,baseCloud::numberOfParticles(),3);
    baseCloud::resetArray(baseCloud::DEMForces_,baseCloud::numberOfParticles(),3);
    for (int i=0;i<baseCloud::nrForceModels();i++) baseCloud::forceM(i).setForce();
}

template <class baseCloud>
void cfdemCloudRec<baseCloud>::setParticleForceField()
{}

template <class baseCloud>
void cfdemCloudRec<baseCloud>::setVectorAverages()
{}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
