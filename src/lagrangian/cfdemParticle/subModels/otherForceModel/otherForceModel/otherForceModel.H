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

    Copyright (C) 2015- Thomas Lichtenegger, JKU Linz, Austria


Class
    otherForceModel

SourceFiles
    otherForceModel.C

\*---------------------------------------------------------------------------*/

#ifndef otherForceModel_H
#define otherForceModel_H

#include "fvCFD.H"
#include "cfdemCloud.H"
#include "probeModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class otherForceModel Declaration
\*---------------------------------------------------------------------------*/

class otherForceModel
{

protected:

    // Protected data
        const dictionary& dict_;

        cfdemCloud& particleCloud_;

public:

    //- Runtime type information
    TypeName("otherForceModel");

    // Declare runtime constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            otherForceModel,
            dictionary,
            (
                const dictionary& dict,
                cfdemCloud& sm
            ),
            (dict,sm)
        );


    // Constructors

        //- Construct from components
        otherForceModel
        (
            const dictionary& dict,
            cfdemCloud& sm
        );


    // Destructor

        virtual ~otherForceModel();


    // Selector

        static autoPtr<otherForceModel> New
        (
            const dictionary& dict,
            cfdemCloud& sm,
            word otherForceType
        );


    // Member Functions

        virtual tmp<volVectorField> exportForceField() = 0;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
