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

#include "momCoupleModel.H"
#include "smoothingModel.H"
#include "meshMotionModel.H"
#include "averagingModel.H"
#include "clockModel.H"
#include "IOModel.H"
#include "voidFractionModel.H"
#include "locateModel.H"
#include "probeModel.H"

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
inline void cfdemCloud::setCG(double cg)
{
    cg_ = cg;
    Info << "cg is set to: " << cg_ << endl;
}

inline bool cfdemCloud::impDEMdrag() const
{
    return impDEMdrag_;
}

inline bool cfdemCloud::impDEMdragAcc() const
{
    return impDEMdragAcc_;
}

inline scalar cfdemCloud::imExSplitFactor() const
{
    return imExSplitFactor_;
}

inline bool cfdemCloud::treatVoidCellsAsExplicitForce() const
{
    return treatVoidCellsAsExplicitForce_;
}

inline scalar cfdemCloud::cg() const
{
    return cg_;
}

inline bool cfdemCloud::ignore() const
{
    return ignore_;
}

inline const fvMesh& cfdemCloud::mesh() const
{
    return mesh_;
}

inline bool cfdemCloud::solveFlow() const
{
    return solveFlow_;
}

inline bool cfdemCloud::verbose() const
{
    return verbose_;
}

inline const IOdictionary& cfdemCloud::couplingProperties() const
{
    return couplingProperties_;
}

inline double ** cfdemCloud::positions() const
{
    return positions_;
}

inline double ** cfdemCloud::velocities() const
{
    return velocities_;
}

inline double ** cfdemCloud::fluidVels() const
{
    return fluidVel_;
}

inline double ** cfdemCloud::fAccs() const
{
    return fAcc_;
}

inline double ** cfdemCloud::impForces() const
{
    return impForces_;
}

inline double ** cfdemCloud::expForces() const
{
    return expForces_;
}

inline double ** cfdemCloud::DEMForces() const
{
    return DEMForces_;
}

inline double ** cfdemCloud::Cds() const
{
    return Cds_;
}

inline double ** cfdemCloud::radii() const
{
    return radii_;
}

inline double ** cfdemCloud::voidfractions() const
{
    return voidfractions_;
}

inline void cfdemCloud::get_radii(double **& values) const
{
//    Info << "set_radii level=" << numberOfParticles_ << endl;
    // make a copy of the array entries
//    for (int i=0;i<numberOfParticles_;i++)
//        radii_[0][i]=values[0][i];
    values=radii_;
}

inline int ** cfdemCloud::cellIDs() const
{
    return cellIDs_;
}

inline void cfdemCloud::get_cellIDs(int **& values) const
{
//    // make a copy of the array entries
//    for (int i=0;i<numberOfParticles_;i++)
//        cellIDs_[0][i]=values[0][i];
    values=cellIDs_;
}

inline double ** cfdemCloud::particleWeights() const
{
    return particleWeights_;
}

inline label cfdemCloud::body(int index) const
{
    return index;
}

inline double cfdemCloud::particleVolume(int index) const
{
    return particleV_[index][0];
}

inline scalar cfdemCloud::radius(int index) const
{
    return radii_[index][0];
}

inline double cfdemCloud::d(int index) const
{
    return 2.*radii_[index][0];
}

inline double cfdemCloud::d32(bool recalc)
{
    if(d32_<0 || recalc)
    {
        scalar Ntot(0);
        scalar Dtot(0);
        scalar r(0);
        for(int index = 0;index <  numberOfParticles(); ++index)
        {
            r=radii_[index][0];
            Ntot+=2*r*r*r;
            Dtot+=r*r;
        }
        d32_=Ntot/Dtot;
    }

    return d32_;
}

inline bool cfdemCloud::getParticleDensities() const
{
    return getParticleDensities_;
}

inline double ** cfdemCloud::particleDensities() const
{
    return particleDensities_;
}

inline double cfdemCloud::particleDensity(label index) const
{
    if(!getParticleDensities_) return -1.0;
    else
    {
        return particleDensities_[index][0];
    }
}

inline bool cfdemCloud::getParticleEffVolFactors() const
{
    return getParticleEffVolFactors_;
}

inline double cfdemCloud::particleEffVolFactor(label index) const
{
    if(!getParticleEffVolFactors_) return -1.0;
    else
    {
        return particleEffVolFactors_[index][0];
    }
}

inline bool cfdemCloud::getParticleTypes() const
{
    return getParticleTypes_;
}

inline int ** cfdemCloud::particleTypes() const
{
    return particleTypes_;
}

inline label cfdemCloud::particleType(label index) const
{
    if(!getParticleDensities_) return -1;
    else
    {
        return particleTypes_[index][0];
    }
}

inline int cfdemCloud::numberOfParticles() const
{
    return numberOfParticles_;
}

inline bool cfdemCloud::numberOfParticlesChanged() const
{
    return numberOfParticlesChanged_;
}

inline int cfdemCloud::numberOfClumps() const
{
    Warning << "cfdemCloud::numberOfClumps() is not called correctly!" << endl;
    return -1;
}

inline bool cfdemCloud::arraysReallocated() const
{
    return arraysReallocated_;
}

inline const wordList& cfdemCloud::forceModels() const
{
    return forceModels_;
}

inline const locateModel& cfdemCloud::locateM() const
{
    return locateModel_;
}

inline momCoupleModel& cfdemCloud::momCoupleM(int i)
{
    return momCoupleModel_[i];
}

inline const momCoupleModel& cfdemCloud::momCoupleM(int i) const
{
    return momCoupleModel_[i];
}

inline dataExchangeModel& cfdemCloud::dataExchangeM()
{
    return dataExchangeModel_();
}

inline const IOModel& cfdemCloud::IOM() const
{
    return IOModel_;
}

inline probeModel& cfdemCloud::probeM()
{
    return probeModel_();
}

inline voidFractionModel& cfdemCloud::voidFractionM()
{
    return voidFractionModel_();
}

inline const voidFractionModel& cfdemCloud::voidFractionM() const
{
    return voidFractionModel_();
}

inline averagingModel& cfdemCloud::averagingM()
{
    return averagingModel_();
}

inline const averagingModel& cfdemCloud::averagingM() const
{
    return averagingModel_();
}

inline const clockModel& cfdemCloud::clockM() const
{
    return clockModel_;
}

inline const smoothingModel& cfdemCloud::smoothingM() const
{
    return smoothingModel_;
}

inline const meshMotionModel& cfdemCloud::meshMotionM() const
{
    return meshMotionModel_;
}

inline const wordList& cfdemCloud::liggghtsCommandModelList() const
{
    return liggghtsCommandModelList_;
}

inline liggghtsCommandModel& cfdemCloud::liggghtsCommand(const label modeli)
{
    return liggghtsCommand_[modeli];
}

inline const turbulenceModel& cfdemCloud::turbulence() const
{
    return turbulence_;
}

}
// ************************************************************************* //
