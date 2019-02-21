/*---------------------------------------------------------------------------*\
License

    This is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This code is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with this code.  If not, see <http://www.gnu.org/licenses/>.

    Copyright (C) 2018- Mathias Vångö, JKU Linz, Austria

Class
    temporalSmoothing

Description
    Smoothens a field using a temporal relaxation approach with a specified
    reference field.

SourceFiles
    temporalSmoothing.C
\*---------------------------------------------------------------------------*/

#ifndef temporalSmoothing_H
#define temporalSmoothing_H

#include "smoothingModel.H"
//#include "wallFvPatch.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class temporalSmoothing Declaration
\*---------------------------------------------------------------------------*/

class temporalSmoothing
:
    public smoothingModel
{

private:

    dictionary propsDict_;
    scalar lowerLimit_;
    scalar upperLimit_;
    bool verbose_;
    word refFieldName_;
    scalar gamma_;

public:

    //- Runtime type information
    TypeName("temporalSmoothing");


    // Constructors

        //- Construct from components
        temporalSmoothing
        (
            const dictionary& dict,
            cfdemCloud& sm
        );

    // Destructor

        ~temporalSmoothing();


    // Member Functions
        bool doSmoothing() const;

        //void dSmoothing(volScalarField&) const;

        void smoothen(volScalarField&) const;

        void smoothen(volVectorField&) const;

        void smoothenReferenceField(volVectorField&) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
