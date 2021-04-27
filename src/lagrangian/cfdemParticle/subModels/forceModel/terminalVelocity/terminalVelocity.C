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

#include "error.H"

#include "terminalVelocity.H"
#include "addToRunTimeSelectionTable.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(terminalVelocity, 0);

addToRunTimeSelectionTable
(
    forceModel,
    terminalVelocity,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
terminalVelocity::terminalVelocity
(
    const dictionary& dict,
    cfdemCloud& sm
)
:
    forceModel(dict,sm),
    propsDict_(dict.subDict(typeName + "Props")),
    interpolate_(propsDict_.lookupOrDefault<bool>("interpolation", false)),
    ignoreCellsName_(propsDict_.lookupOrDefault<word>("ignoreCellsName","none")),
    ignoreCells_(),
    existIgnoreCells_(true),
    mesh_(sm.mesh()),
    turbKineticEnergyFieldName_(propsDict_.lookupOrDefault<word>("turbKineticEnergyFieldName","")),
    turbDissipationRateFieldName_(propsDict_.lookupOrDefault<word>("turbDissipationRateFieldName","")),
    turbKineticEnergy_(NULL),
    turbDissipationRate_(NULL),
//    recurrenceBaseName_(propsDict_.lookupOrDefault<word>("recurrenceBaseName","")),
//    recurrenceBase_(NULL),
    existturbDissipationRateInObjReg_(false),
    delta(  IOobject
                (
                      "delta",
                      mesh_.time().timeName(),
                      mesh_,
                      IOobject::NO_READ,
                      IOobject::NO_WRITE
                ),
                mesh_,
                dimensionedScalar("delta", dimLength, 1.0)),
    wallIndicatorField_
    (   IOobject
        (
            "wallIndicator",
            sm.mesh().time().timeName(),
            sm.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        sm.mesh(),
        dimensionedScalar("zero", dimensionSet(0,0,0,0,0,0,0), 0.0)
    ),
    liquidViscosity_(propsDict_.lookupOrDefault<scalar>("liquidViscosity", 1.0)),
    dragReductionFactor_(propsDict_.lookupOrDefault<scalar>("dragReductionFactor", 0.0)),
    terminalVel_(propsDict_.lookupOrDefault<scalar>("terminalVelocity", 0.0)),
    risingVelDir_(propsDict_.lookupOrDefault<scalar>("risingVelDir", 2))
{
    if(ignoreCellsName_ != "none")
    {
       ignoreCells_.set(new cellSet(particleCloud_.mesh(),ignoreCellsName_));
       Info << "terminalVelocity: ignoring rising velocity in cellSet " << ignoreCells_().name() <<
        " with " << ignoreCells_().size() << " cells." << endl;
    }
    else existIgnoreCells_ = false;

/*
    if (recurrenceBaseName_ != "")
    {
        Info << "creating the recBase pointer" << endl;
        Foam::recBase& rBase_ = const_cast<recBase&>(sm.mesh().lookupObject<Foam::recBase>(recurrenceBaseName_));
        recurrenceBase_ = &rBase_;
        existturbDissipationRateInObjReg_ = true;
    }
*/
       turbDissipationRate_ = new volScalarField
        (
            IOobject
            (
                "turbDissipationRate",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("zero", dimensionSet(0,2,-3,0,0), 0)
        );

       if (turbKineticEnergyFieldName_ != "")
       {
           turbKineticEnergy_ = new volScalarField
                   (
                       IOobject
                       (
                           "turbKineticEnergy",
                           mesh_.time().timeName(),
                           mesh_,
                           IOobject::NO_READ,
                           IOobject::NO_WRITE
                       ),
                       mesh_,
                       dimensionedScalar("zero", dimensionSet(0,2,-2,0,0), 0)
                   );
           delta.primitiveFieldRef()=Foam::pow(mesh_.V(),1.0/3.0);
       }

       // define a field to indicate if a cell is next to boundary
        label cellI = -1;
        forAll (mesh_.boundary(),patchI)
        {
            word patchName = mesh_.boundary()[patchI].name();
            if (patchName.rfind("procB",0) == 0) continue;

            forAll(mesh_.boundary()[patchI], faceI)
            {
                cellI = mesh_.boundary()[patchI].faceCells()[faceI];
                wallIndicatorField_[cellI] = 1.0;
            }
        }

}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

terminalVelocity::~terminalVelocity()
{
        if (!existturbDissipationRateInObjReg_) delete turbDissipationRate_;
}

// * * * * * * * * * * * * * * * private Member Functions  * * * * * * * * * * * * * //

bool terminalVelocity::ignoreCell(label cell) const
{
    if (!existIgnoreCells_) return false;
    else return ignoreCells_()[cell];
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void terminalVelocity::setForce() const
{
    updateEpsilon();

    vector position(0,0,0);
    scalar Urising=0.0;
    label cellI=-1;
    scalar radius=0.0;
    scalar epsilon=0.0;
    vector Uparticle(0,0,0);

    label patchID = -1;
    label faceIGlobal = -1;
    scalar velProjection = 0.0;
    vector faceINormal = vector::zero;
    word patchName("");

    interpolationCellPoint<scalar> turbDissipationRateInterpolator_(*turbDissipationRate_);
    for(int index = 0;index <  particleCloud_.numberOfParticles(); ++index)
    {
            cellI = particleCloud_.cellIDs()[index][0];
            Urising =0.0;
            Uparticle = vector::zero;
            if (cellI > -1 && !ignoreCell(cellI)) // particle found
            {
                if (interpolate_)
                {
                    position = particleCloud_.position(index);
                    epsilon = turbDissipationRateInterpolator_.interpolate(position,cellI);
                }
                else
                {
                    epsilon = (*turbDissipationRate_)[cellI];
                }

                    position = particleCloud_.position(index);
                    radius = particleCloud_.radius(index);

                    // d * kolmogorov length scale
                    scalar dLambda = 2*radius*pow(epsilon,0.25)/pow(liquidViscosity_,0.75);

                    Urising = terminalVel_ / Foam::sqrt(1+(dragReductionFactor_*pow(dLambda,3)) );
                    particleCloud_.particleConvVels()[index][risingVelDir_] += Urising;

                    // read the new particle velocity
                    for(int j=0;j<3;j++)
                    {
                        Uparticle[j]= particleCloud_.particleConvVels()[index][j];
                    }

                    // prevent particles being pushed through walls
                    // check if cell is adjacent to wall and remove the normal velocity to the wall
                    if (wallIndicatorField_[cellI] > 0.5)
                    {
                        const cell& faces = mesh_.cells()[cellI];
                        forAll (faces, faceI)
                        {
                            faceIGlobal = faces[faceI];
                            patchID = mesh_.boundaryMesh().whichPatch(faceIGlobal);
                            if (patchID < 0) continue;
                            patchName = mesh_.boundary()[patchID].name();

                            if (patchName.rfind("procB",0) == 0) continue;

                            faceINormal = mesh_.Sf()[faceIGlobal];
                            faceINormal /= mag(faceINormal);
                            velProjection = faceINormal&Uparticle;
                            if (velProjection > 0.0)
                            {
                                // removes the value normal to the face
                                for(int j=0;j<3;j++)
                                {
                                    particleCloud_.particleConvVels()[index][j] -= velProjection*faceINormal[j];
                                }
                            }
                        }
                    }


	    }

            if (index >0 && index <2)
            {
                Pout << "cellI = " << cellI << endl;
                Pout << "index = " << index << endl;
                Pout << "epsilon = " << epsilon << endl;
                Pout << "rising velocity = " << Urising << endl;
            }
    }
}

void terminalVelocity::updateEpsilon() const
{
    if (!existturbDissipationRateInObjReg_)
    {
        Info << "epsilon is calculated from the turbulence model. " << endl;
        *turbDissipationRate_ = particleCloud_.turbulence().epsilon()();
    }
/*
    else if (recurrenceBaseName_ != "")
    {
         Info << "updating epsilon from the database. " << endl;
         label currTimeIndex = recurrenceBase_ ->recM().currentTimeIndex();
         Info << "current time index = " << currTimeIndex << endl;

         if(turbDissipationRateFieldName_ != "")
         {
             *turbDissipationRate_ = recurrenceBase_ -> recM().exportVolScalarField(turbDissipationRateFieldName_,currTimeIndex);
         }

         else if(turbKineticEnergyFieldName_ !="")
         {
             *turbKineticEnergy_ = recurrenceBase_ -> recM().exportVolScalarField(turbKineticEnergyFieldName_,currTimeIndex);
             *turbDissipationRate_ = 1.048*Foam::pow(*turbKineticEnergy_,1.5)/delta;
         }
     }
    else
    {
        FatalError <<"turbulent dissipation enery must be calculated either from the turbulence model or recurrence dataBase!\n" << abort(FatalError);
    }
*/
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
