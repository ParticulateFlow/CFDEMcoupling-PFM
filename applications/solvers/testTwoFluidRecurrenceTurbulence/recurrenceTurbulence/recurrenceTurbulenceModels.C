/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2014-2015 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "phaseCompressibleTurbulenceModel.H"
#include "addToRunTimeSelectionTable.H"
#include "makeTurbulenceModel.H"

#include "laminar.H"
#include "RASModel.H"
#include "LESModel.H"

#include "recurrenceKEpsilon.H"
#include "recurrenceKOmega.H"
#include "recurrenceSmagorinsky.H"



// Instructions for OpenFOAM-5.0 
/*makeTurbulenceModelTypes
(
    volScalarField,
    volScalarField,
    compressibleTurbulenceModel,
    PhaseCompressibleTurbulenceModel,
    ThermalDiffusivity,
    phaseModel
);


makeTurbulenceModel
(phaseModelPhaseCompressibleTurbulenceModel, RAS, recurrenceKEpsilon);

makeTurbulenceModel
(phaseModelPhaseCompressibleTurbulenceModel, RAS, recurrenceKOmega);

makeTurbulenceModel
(phaseModelPhaseCompressibleTurbulenceModel, LES, recurrenceSmagorinsky);
*/



// Instructions for OpenFOAM-4.0 
makeBaseTurbulenceModel
(
    volScalarField,
    volScalarField,
    compressibleTurbulenceModel,
    PhaseCompressibleTurbulenceModel,
    ThermalDiffusivity,
    phaseModel
);

#define makeRASModel(Type)                                            \
    makeTurbulenceModel                                               \
    (phaseModelPhaseCompressibleTurbulenceModel, RAS, Type)

#define makeLESModel(Type)                                            \
    makeTurbulenceModel                                               \
    (phaseModelPhaseCompressibleTurbulenceModel, LES, Type)

makeRASModel(recurrenceKEpsilon);

makeRASModel(recurrenceKOmega);

makeLESModel(recurrenceSmagorinsky);

// ************************************************************************* //
