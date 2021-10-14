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


Description
    model for species diffusion coefficient in fluid phase in presence of particles

Class
    diffCoeffModel

SourceFiles
    diffCoeffModel.C

\*---------------------------------------------------------------------------*/

#ifndef diffCoeffModel_H
#define diffCoeffModel_H

#include "fvCFD.H"
#include "cfdemCloud.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class diffCoeffModel Declaration
\*---------------------------------------------------------------------------*/

class diffCoeffModel
{

protected:

    // Protected data
        const dictionary& dict_;

        cfdemCloud& particleCloud_;

        IOdictionary transportProperties_;

        mutable volScalarField D0Field_;

public:

    //- Runtime type information
    TypeName("diffCoeffModel");

    // Declare runtime constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            diffCoeffModel,
            dictionary,
            (
                const dictionary& dict,
                cfdemCloud& sm
            ),
            (dict,sm)
        );


    // Constructors

        //- Construct from components
        diffCoeffModel
        (
            const dictionary& dict,
            cfdemCloud& sm
        );


    // Destructor

        virtual ~diffCoeffModel();


    // Selector

        static autoPtr<diffCoeffModel> New
        (
            const dictionary& dict,
            cfdemCloud& sm
        );


    // Member Functions

        virtual tmp<volScalarField> diffCoeff() const = 0;

        const volScalarField& D0Field() const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
