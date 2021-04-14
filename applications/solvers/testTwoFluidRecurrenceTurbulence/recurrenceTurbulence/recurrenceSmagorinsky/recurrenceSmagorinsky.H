/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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

Class
    Foam::LESModels::recurrenceSmagorinsky

Description
    The recurrenceSmagorinsky SGS model.

    This model provides Smagorinsky's turbulence quantities.
    To be used by recurrence solvers.

SourceFiles
    recurrenceSmagorinsky.C

\*---------------------------------------------------------------------------*/

#ifndef recurrenceSmagorinsky_H
#define recurrenceSmagorinsky_H

#include "LESModel.H"
#include "LESeddyViscosity.H"
#include "phaseCompressibleTurbulenceModel.H"
#include "EddyDiffusivity.H"

#include "recurrenceTurbulenceModel.H"

#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace LESModels
{

/*---------------------------------------------------------------------------*\
                           Class recurrenceSmagorinsky Declaration
\*---------------------------------------------------------------------------*/

class recurrenceSmagorinsky
:
    public LESeddyViscosity<EddyDiffusivity<phaseCompressibleTurbulenceModel>>,
    public recurrenceTurbulenceModel
{
    // Private Member Functions

        // Disallow default bitwise copy construct and assignment
        recurrenceSmagorinsky(const recurrenceSmagorinsky&);
        void operator=(const recurrenceSmagorinsky&);


protected:

    // Protected data

        dimensionedScalar Ck_;


    // Protected Member Functions

        //- Return SGS kinetic energy
        //  calculated from the given velocity gradient
        tmp<volScalarField> k(const tmp<volTensorField>& gradU) const;

        //- Update the SGS eddy viscosity
        virtual void correctNut();


public:


    //- Runtime type information
    TypeName("recurrenceSmagorinsky");


    // Constructors

        //- Construct from components
        recurrenceSmagorinsky
        (
            const volScalarField& alpha,
            const volScalarField& rho,
            const volVectorField& U,
            const surfaceScalarField& alphaRhoPhi,
            const surfaceScalarField& phi,
            const phaseModel& transport,
            const word& propertiesName = turbulenceModel::propertiesName,
            const word& type = typeName
        );


    //- Destructor
    virtual ~recurrenceSmagorinsky()
    {}


    // Member Functions

        //- Read model coefficients if they have changed
        virtual bool read();

        //- Return SGS kinetic energy
        virtual tmp<volScalarField> k() const
        {
            return k(fvc::grad(this->U_));
        }

        //- Return sub-grid disipation rate
        virtual tmp<volScalarField> epsilon() const;

        //- Correct Eddy-Viscosity and related properties
        virtual void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
