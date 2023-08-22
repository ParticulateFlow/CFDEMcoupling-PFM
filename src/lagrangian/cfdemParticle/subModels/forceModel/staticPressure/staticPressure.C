/*---------------------------------------------------------------------------*\
    CFDEMcoupling - Open Source CFD-DEM coupling

    Copyright (C) 2023  Behrad Esgandari, JKU Linz, Austria
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

\*---------------------------------------------------------------------------*/

#include "error.H"
#include "cellSet.H"
#include "staticPressure.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(staticPressure, 0);

addToRunTimeSelectionTable
(
    forceModel,
    staticPressure,
    dictionary
);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
// Construct from components
staticPressure::staticPressure
(
    const dictionary& dict,
    cfdemCloud& sm
)
:
    forceModel(dict,sm),
    propsDict_(dict.subDict(typeName + "Props")),
    rhoGas_(readScalar(propsDict_.lookup ("rhoGas"))),
    g_(propsDict_.lookup ("g_DEM")),
    voidfractionFieldName_(propsDict_.lookup("voidfractionFieldName")),
    voidfraction_(sm.mesh().lookupObject<volScalarField> (voidfractionFieldName_)),
    rhoPart_(readScalar(propsDict_.lookup ("rhoPart"))),
    solidFraction_(readScalar(propsDict_.lookup ("DomainSolidVolumeFraction")))
{
    // init force sub model
    setForceSubModels(propsDict_);
    // define switches which can be read from dict
    forceSubM(0).setSwitchesList(SW_TREAT_FORCE_EXPLICIT,true); // activate treatExplicit switch
    forceSubM(0).setSwitchesList(SW_TREAT_FORCE_DEM,true); // activate treatForceDEM switch
    forceSubM(0).setSwitchesList(SW_INTERPOLATION,true); // activate search for interpolate switch
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //
staticPressure::~staticPressure()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
void staticPressure::setForce() const
{
    #include "setupProbeModel.H"
    label cellI;
    scalar rhoMix_ = solidFraction_*rhoPart_ + (1.0-solidFraction_)*rhoGas_;
    vector force;

    for(int index = 0; index < particleCloud_.numberOfParticles(); ++index)
    {
        force = vector::zero;
        cellI = particleCloud_.cellIDs()[index][0];
        if (cellI > -1) // particle Found
        {
            scalar Vs = particleCloud_.particleVolume(index);
            // set force on particle
            force = -Vs * rhoMix_ * g_;
        }
        // write particle based data to global array
        forceSubM(0).partToArray(index,force,vector::zero);
     }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
