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
    Foam::recurrenceTurbulenceModel

Description
    recurrence-based turbulence model base class

    This turbulence model class provides a framework for derived turbulence
    models to provide turbulence quantities that were computed elsewhere, 
    i.e. taken from the recurrence dataBase.
    
    A concrete recurrence-based turbulence model implementation, e.g. recurrenceKEpsilon, 
    needs to be derived from the base class of the original turbulence model 
    and this class. This class provides the link to the recurrence model. The original 
    base class provides the proper turbulence modelling interfaces of OpenFOAM's 
    turbulence modelling framework.
    
    To be used by recurrence solvers.

SourceFiles
    recurrenceTurbulenceModel.C

\*---------------------------------------------------------------------------*/

#ifndef recurrenceTurbulenceModel_H
#define recurrenceTurbulenceModel_H

#include "recBase.H"
#include "recModel.H"

#include "autoPtr.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class recurrenceTurbulenceModel Declaration
\*---------------------------------------------------------------------------*/

class recurrenceTurbulenceModel
{
    // Private Member Functions


        //- Disallow default bitwise copy construct
        recurrenceTurbulenceModel(const recurrenceTurbulenceModel&);

        //- Disallow default bitwise assignment
        void operator=(const recurrenceTurbulenceModel&);
        

protected:

    // Protected data

        recBase* recurrenceBasePtr_;
        const word group_;
        
    
    // Protected Member Functions


public:

    //- Runtime type information
    TypeName("recurrenceTurbulenceModel");


    // Constructors

        //- Construct from components
        recurrenceTurbulenceModel
        (
            const word group
        );


    //- Destructor
    virtual ~recurrenceTurbulenceModel();


    // Member Functions
        
        void setRecurrenceBasePtr(recBase* recurrenceBasePtr)
        {
            recurrenceBasePtr_ = recurrenceBasePtr;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
