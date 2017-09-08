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

#include "cfdemCloudMS.H"

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline label cfdemCloudMS::body(int index) const
{
    return bodies_[0][index]-1;
}

inline double cfdemCloudMS::particleVolume(int index) const
{
    int ind = body(index); // particle to clump ID

    double Vp(SMALL);
    if(ind >= 0 && ind <numberOfClumps())
        Vp=clumpVol_[ind][0];
    else
        Warning << "ind=" << ind << endl;

    int nR(nrigid(ind));
    if(nR>0) Vp/=nR;
    return Vp;
}

inline vector cfdemCloudMS::positionCM(int index)
{
    vector pos;
    for(int i=0;i<3;i++) pos[i] = positionsCM()[index][i];
    return pos;
}

inline vector cfdemCloudMS::velocityCM(int index)
{
    vector vel;
    for(int i=0;i<3;i++) vel[i] = velocitiesCM()[index][i];
    return vel;
}

inline label cfdemCloudMS::cellIDCM(int index)
{
    return cellIDsCM_[index][0];
}

inline label cfdemCloudMS::nrigid(int index) const
{
    return nrigids_[0][index];
}

inline label cfdemCloudMS::nrForceModels() const
{
    return forceModels_.size();
}

inline double **& cfdemCloudMS::positionsCM() const
{
    return positionsCM_;
}

inline double **& cfdemCloudMS::velocitiesCM() const
{
    return velocitiesCM_;
}

inline int cfdemCloudMS::nClumpTypes() const
{
    return nClumpTypes_;
}

inline int **& cfdemCloudMS::clumpType() const
{
    return clumpType_;
}

inline double **& cfdemCloudMS::clumpVol() const
{
    return clumpVol_;
}

inline double **& cfdemCloudMS::clumpDH() const
{
    return clumpDH_;
}

inline double **& cfdemCloudMS::impForcesCM() const
{
    return impForcesCM_;
}

inline double **& cfdemCloudMS::expForcesCM() const
{
    return expForcesCM_;
}

inline double **& cfdemCloudMS::DEMForcesCM() const
{
    return DEMForcesCM_;
}

inline int cfdemCloudMS::numberOfClumps() const
{
    return numberOfClumps_;
}

inline bool cfdemCloudMS::numberOfClumpsChanged() const
{
    return numberOfClumpsChanged_;
}

inline void cfdemCloudMS::setManDHdev(Switch v)
{
    manDHdev_ = v;
}

inline void cfdemCloudMS::setDHbyV(scalarList v)
{
    dHbyV_ = v;
}

inline scalarList& cfdemCloudMS::dHbyV()
{
    return dHbyV_;
}

inline bool cfdemCloudMS::useforcePerClump() const
{
    return useforcePerClump_;
}

}
// ************************************************************************* //
