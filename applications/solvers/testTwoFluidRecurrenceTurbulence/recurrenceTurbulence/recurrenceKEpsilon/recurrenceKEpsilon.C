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

#include "recurrenceKEpsilon.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


void Foam::RASModels::recurrenceKEpsilon::correctNut()
{
    this->nut_ = Cmu_*sqr(k_)/epsilon_;
    this->nut_.correctBoundaryConditions();
}



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::RASModels::recurrenceKEpsilon::recurrenceKEpsilon
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
    epsilon_
    (
        IOobject
        (
            IOobject::groupName("epsilon", U.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_,
        dimensionedScalar("eps0", dimensionSet(0,2,-3,0,0), 0.0)
    )
{
    if (type == typeName)
    {
        printCoeffs(type);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //


Foam::RASModels::recurrenceKEpsilon::~recurrenceKEpsilon()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::RASModels::recurrenceKEpsilon::read()
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


void Foam::RASModels::recurrenceKEpsilon::correct()
{
    // update turbulence fields
    recurrenceBasePtr_->recM().exportVolScalarField("k."+group_, this->k_);
    recurrenceBasePtr_->recM().exportVolScalarField("epsilon."+group_, this->epsilon_);
    
    // update nut
    correctNut();
}


void Foam::RASModels::recurrenceKEpsilon::validate()
{
    /*
        Check whether k and epsilon are included in the dataBase.
        The check only makes sure that these fields are included in the
            volScalarFields list of recProperties.
        Whether the fields are actually contained in the dataBase is 
            done by the recurrenceModel itself.
    */
    bool foundK(false);
    bool foundEpsilon(false);
    
    wordList fieldNames(recurrenceBasePtr_->recM().volScalarFieldNames());
    
    forAll(fieldNames, i)
    {
        word curFieldName = fieldNames[i];
        
        if (curFieldName == k_.name())
        {
            foundK = true;
        }
        
        if (curFieldName == epsilon_.name())
        {
            foundEpsilon = true;
        }
        
        
    }
    
    if (not (foundK and foundEpsilon))
    {
        FatalError
            << "Fields " << k_.name() << " and " << epsilon_.name()
            << " not specified in the volScalarFields list of recProperties!" << nl
            << "volScalarFields : " << fieldNames << nl
            << "Add these fields and make sure they are contained in the dataBase." << nl
            << exit(FatalError);
    }
}


// ************************************************************************* //
