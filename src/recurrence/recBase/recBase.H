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
#ifndef recBase_H
#define recBase_H

// choose version
#include <vector>

#include "fvCFD.H"
#include "IFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// forward declarations
class recModel;
class recNorm;
class recPath;
class recStatAnalysis;

/*---------------------------------------------------------------------------*\
                           Class recBase Declaration
\*---------------------------------------------------------------------------*/

class recBase
{

protected:

    const fvMesh& mesh_;

    IOdictionary recProperties_;

    autoPtr<recModel> recModel_;

    autoPtr<recNorm> recNorm_;

    autoPtr<recPath> recPath_;

    autoPtr<recStatAnalysis> recStatAnalysis_;

    // in case of substepping: couple passive quantity at this substep
    label couplingSubStep_;

public:

      // Constructors

    //- Construct from mesh and a list of particles
    recBase
    (
            const fvMesh& mesh
    );
    //- Construct from mesh and a list of recProperties
    recBase
    (
            const fvMesh& mesh, const word recDictName_
    );

    //- Destructor
    virtual ~recBase();

    // public Member Functions

    const fvMesh& mesh() const;

    recModel& recM();

    recPath& recP();

    recStatAnalysis& recStatA();

    void extendPath();

    void updateRecFields();

    label couplingSubStep() const;

};

}

#endif
