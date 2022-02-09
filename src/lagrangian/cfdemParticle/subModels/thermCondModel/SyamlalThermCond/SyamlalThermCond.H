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
    thermal conductivity of fluid phase in presence of particles according to
    Syamlal, M. and Gidaspow, M. AIChE Journal 31.1 (1985)

Class
    SyamlalThermCond

SourceFiles
    SyamlalThermCond.C

\*---------------------------------------------------------------------------*/


#ifndef SyamlalThermCond_H
#define SyamlalThermCond_H

#include "thermCondModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class SyamlalThermCond Declaration
\*---------------------------------------------------------------------------*/

class SyamlalThermCond
:
    public thermCondModel
{

private:

    dictionary propsDict_;

    word voidfractionFieldName_;

    const volScalarField& voidfraction_;

public:

    //- Runtime type information
    TypeName("SyamlalThermCond");


    // Constructors

        //- Construct from components
        SyamlalThermCond
        (
            const dictionary& dict,
            cfdemCloud& sm
        );

    // Destructor

        ~SyamlalThermCond();


    // Member Functions

        void calcThermCond();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
