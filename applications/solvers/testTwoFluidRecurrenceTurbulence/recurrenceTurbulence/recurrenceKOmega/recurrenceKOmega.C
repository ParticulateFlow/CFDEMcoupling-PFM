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

\*---------------------------------------------------------------------------*/

#include "recurrenceKOmega.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


void Foam::RASModels::recurrenceKOmega::correctNut()
{
    this->nut_ = k_/omega_;
    this->nut_.correctBoundaryConditions();
}



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::RASModels::recurrenceKOmega::recurrenceKOmega
(
    const volScalarField& alpha,
    const volScalarField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const transportModel& phase,
    const word& propertiesName,
    const word& type

)
:
    eddyViscosity<RASModel<EddyDiffusivity<phaseCompressibleTurbulenceModel>>>
    (
        type,
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        phase,
        propertiesName
    ),
    recurrenceTurbulenceModel(U.group()),
    Cmu_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cmu",
            this->coeffDict_,
            0.09
        )
    ),
    k_
    (
        IOobject
        (
            IOobject::groupName("k", U.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_,
        dimensionedScalar("k0", dimensionSet(0,2,-2,0,0), 0.0)
    ),
    omega_
    (
        IOobject
        (
            IOobject::groupName("omega", U.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_,
        dimensionedScalar("om0", dimensionSet(0,0,-1,0,0), 0.0)
    )
{
    if (type == typeName)
    {
        printCoeffs(type);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //


Foam::RASModels::recurrenceKOmega::~recurrenceKOmega()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::RASModels::recurrenceKOmega::read()
{
    if
    (
        eddyViscosity<RASModel<EddyDiffusivity<phaseCompressibleTurbulenceModel>>>::read()
    )
    {
        Cmu_.readIfPresent(this->coeffDict());
        
        return true;
    }
    else
    {
        return false;
    }
}


void Foam::RASModels::recurrenceKOmega::correct()
{
    // update turbulence fields
    recurrenceBasePtr_->recM().exportVolScalarField("k."+group_, this->k_);
    recurrenceBasePtr_->recM().exportVolScalarField("omega."+group_, this->omega_);
    
    // update nut
    correctNut();
}


void Foam::RASModels::recurrenceKOmega::validate()
{
    /*
        Check whether k and omega are included in the dataBase.
        The check only makes sure that these fields are included in the
            volScalarFields list of recProperties.
        Whether the fields are actually contained in the dataBase is 
            done by the recurrenceModel itself.
    */
    bool foundK(false);
    bool foundOmega(false);
    
    wordList fieldNames(recurrenceBasePtr_->recM().volScalarFieldNames());
    
    forAll(fieldNames, i)
    {
        word curFieldName = fieldNames[i];
        
        if (curFieldName == k_.name())
        {
            Info << "Found " << k_.name()<< endl;
            foundK = true;
        }
        
        if (curFieldName == omega_.name())
        {
            Info << "Found " << omega_.name()<< endl;
            foundOmega = true;
        }
    }
    
    if (not (foundK and foundOmega))
    {
        FatalError
            << "Fields " << k_.name() << " and " << omega_.name()
            << " not specified in the volScalarFields list of recProperties!" << nl
            << "volScalarFields : " << fieldNames << nl
            << "Add these fields and make sure they are contained in the dataBase." << nl
            << exit(FatalError);
    }
}


// ************************************************************************* //
