/*---------------------------------------------------------------------------*\
    CFDEMcoupling academic - Open Source CFD-DEM coupling
    
    Contributing authors:
    Thomas Lichtenegger, Sanaz Abbasi
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

#ifndef terminalVelocity_H
#define terminalVelocity_H

#include "forceModel.H"
#include "interpolationCellPoint.H"
#include "cellSet.H"
#include "recBase.H"
#include "recModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class freeStreaming Declaration
\*---------------------------------------------------------------------------*/

class terminalVelocity
:
    public forceModel
{
private:
    dictionary propsDict_;
    
    bool interpolate_;

    const fvMesh& mesh_;

    // ignore streaming velocity in region
    word ignoreCellsName_;

    autoPtr<cellSet> ignoreCells_;

    bool existIgnoreCells_;

    word turbDissipationRateFieldName_;

    word turbKineticEnergyFieldName_;

    word recurrenceBaseName_;

    volScalarField* turbKineticEnergy_;

    volScalarField* turbDissipationRate_;

    mutable volScalarField delta;

    mutable volScalarField wallIndicatorField_;

    bool existturbDissipationRateInObjReg_;

    scalar liquidViscosity_;

    // drag reduction factor in turbulent flow
    scalar dragReductionFactor_;

    // terminal rising velocity of the bubble in the resting liquid
    scalar terminalVel_;

    int risingVelDir_;

    bool ignoreCell(label) const;

    void updateEpsilon() const;

    recBase* recurrenceBase_;

public:

    //- Runtime type information
    TypeName("terminalVelocity");


    // Constructors

        //- Construct from components
        terminalVelocity
        (
            const dictionary& dict,
            cfdemCloud& sm
        );

    // Destructor

        ~terminalVelocity();


    // Member Functions
        void setForce() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
