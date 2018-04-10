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

#include "error.H"

#include "evaluateFluctuations.H"
#include "addToRunTimeSelectionTable.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(evaluateFluctuations, 0);

addToRunTimeSelectionTable
(
    forceModel,
    evaluateFluctuations,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
evaluateFluctuations::evaluateFluctuations
(
    const dictionary& dict,
    cfdemCloud& sm
)
:
    forceModel(dict,sm),
    propsDict_(dict.subDict(typeName + "Props")),
    interpolate_(propsDict_.lookupOrDefault<bool>("interpolation", false)),
    activeCellsName_(propsDict_.lookupOrDefault<word>("activeCellsName","all")),
    activeCells_(),
    voidfractionFieldName_(propsDict_.lookupOrDefault<word>("voidfractionFieldName","voidfraction")),
    voidfraction_(sm.mesh().lookupObject<volScalarField> (voidfractionFieldName_)),
    voidfractionTest_(voidfraction_),
    voidfractionRecFieldName_(propsDict_.lookupOrDefault<word>("voidfractionRecFieldName","voidfractionRec")),
    voidfractionRec_(sm.mesh().lookupObject<volScalarField> (voidfractionRecFieldName_)),
    critVoidfraction_(propsDict_.lookupOrDefault<scalar>("critVoidfraction", 1.0)),
    dtDEM_(particleCloud_.dataExchangeM().DEMts())
{
    if(activeCellsName_ != "all")
    {
       activeCells_.set(new cellSet(particleCloud_.mesh(),activeCellsName_));
       Info << "evaluateFluctuations: evaluating fluctuations in cellSet " << activeCells_().name() <<
        " with " << activeCells_().size() << " cells." << endl;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

evaluateFluctuations::~evaluateFluctuations()
{
}

// * * * * * * * * * * * * * * * private Member Functions  * * * * * * * * * * * * * //

bool evaluateFluctuations::activeCell(label cell) const
{
    if (activeCellsName_ == "all") return true;
    else return activeCells_()[cell];
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void evaluateFluctuations::setForce() const
{
    voidfractionTest_ = voidfraction_;
    vector position(0,0,0);
    scalar voidfraction(0.0);
    scalar voidfractionRec(0.0);
    scalar deltaVoidfrac(0.0);

    vector flucU(0,0,0);
    label cellI=0;
    bool undoStep = false;
   
    interpolationCellPoint<scalar> voidfractionInterpolator_(voidfractionTest_);
    interpolationCellPoint<scalar> voidfractionRecInterpolator_(voidfractionRec_);
    
    for(int index = 0;index <  particleCloud_.numberOfParticles(); ++index)
    {
            cellI = particleCloud_.cellIDs()[index][0];
            flucU = vector(0,0,0);
            voidfraction = 0.0;
            voidfractionRec = 0.0;
            deltaVoidfrac = 0.0;
            undoStep = false;

            if (cellI > -1 && activeCell(cellI)) // particle found
            {
                // particles in empty regions follow trajectories subject to gravity
                if(voidfractionRec_[cellI] < critVoidfraction_)
                {
                    if( interpolate_ )
                    {
                      position = particleCloud_.position(index);
                      voidfraction = voidfractionInterpolator_.interpolate(position,cellI);
                      voidfractionRec = voidfractionRecInterpolator_.interpolate(position,cellI);
                    }
                    else
                    {
                        voidfraction = voidfractionTest_[cellI];
                        voidfractionRec = voidfractionRec_[cellI];
                    }
                    

                    for(int j=0;j<3;j++)
                    {
                        flucU[j] = particleCloud_.particleFlucVels()[index][j];
                    }

                    vector newpos = particleCloud_.position(index) + flucU * dtDEM_;
                    label newcell = particleCloud_.mesh().findCell(newpos);
                    if(newcell > -1 && activeCell(newcell))
                    {
                        // avoid in-cell mixing
                        if(newcell == cellI) undoStep = true;
                        else
                        {
                            // is the move making anything better?
                            scalar volP = particleCloud_.particleVolume(index);
                            scalar cellVolNew = particleCloud_.mesh().V()[newcell];
                            scalar cellVolOld = particleCloud_.mesh().V()[cellI];
                            deltaVoidfrac=voidfractionRec-voidfraction;
                            scalar oldErr = deltaVoidfrac - volP/cellVolOld;
                            scalar newErr = voidfractionRec_[newcell] - voidfractionTest_[newcell] + volP/cellVolNew;
                            if(newErr > oldErr) undoStep = true;
                            else
                            {
                                voidfractionTest_[newcell] -= volP/cellVolNew;
                                voidfractionTest_[cellI] += volP/cellVolOld;
                            }
                        }
                    }

                    if(undoStep)
                    {
                        for(int j=0;j<3;j++)
                        {
                            particleCloud_.particleFlucVels()[index][j] = 0.0;
                        }
                    }
                }
            }
    }

}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
