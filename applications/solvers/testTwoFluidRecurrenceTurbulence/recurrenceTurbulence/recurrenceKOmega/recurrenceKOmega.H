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
    Foam::recurrenceKOmega

Description
    recurrence-based kOmega turbulence model

    This model provides kOmega's turbulence quantities that were computed 
    elsewhere, i.e. taken from the recurrence dataBase
    To be used by recurrence solvers

SourceFiles
    recurrenceKOmega.C

\*---------------------------------------------------------------------------*/

#ifndef recurrenceKOmega_H
#define recurrenceKOmega_H

#include "RASModel.H"
#include "eddyViscosity.H"
#include "phaseCompressibleTurbulenceModel.H"
#include "EddyDiffusivity.H"

#include "recurrenceTurbulenceModel.H"

#include "autoPtr.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
                           Class recurrenceKOmega Declaration
\*---------------------------------------------------------------------------*/

class recurrenceKOmega
:
    public eddyViscosity<RASModel<EddyDiffusivity<phaseCompressibleTurbulenceModel>>>,
    public recurrenceTurbulenceModel
{
    // Private data

    // Private Member Functions


        //- Disallow default bitwise copy construct
        recurrenceKOmega(const recurrenceKOmega&);

        //- Disallow default bitwise assignment
        void operator=(const recurrenceKOmega&);
        

protected:

    // Protected data

        // Model coefficients

            dimensionedScalar Cmu_;

        // Fields

            volScalarField k_;
            volScalarField omega_;
    
    // Protected Member Functions

        virtual void correctNut();


public:


    //- Runtime type information
    TypeName("recurrenceKOmega");


    // Constructors

        //- Construct from components
        recurrenceKOmega
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
    virtual ~recurrenceKOmega();


    // Member Functions

        //- Re-read model coefficients if they have changed
        virtual bool read();
        
        
        //- Return the turbulence kinetic energy
        virtual tmp<volScalarField> k() const
        {
            return k_;
        }

        //- Return the turbulence specific dissipation rate
        virtual tmp<volScalarField> omega() const
        {
            return omega_;
        }

        //- Return the turbulence kinetic energy dissipation rate
        virtual tmp<volScalarField> epsilon() const
        {
            return tmp<volScalarField>
            (
                new volScalarField
                (
                    IOobject
                    (
                        "epsilon",
                        this->mesh_.time().timeName(),
                        this->mesh_
                    ),
                    Cmu_*k_*omega_,
                    omega_.boundaryField().types()
                )
            );
        }

        //- Update the turbulent fields
        virtual void correct();
        
        //- Check model settings
        virtual void validate();
        
        
    // Setters
        
        void setRecurrenceBasePtr(recBase* recurrenceBasePtr)
        {
            recurrenceBasePtr_ = recurrenceBasePtr;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
